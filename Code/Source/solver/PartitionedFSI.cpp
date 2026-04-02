// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "PartitionedFSI.h"
#include "fsi_coupling.h"
#include "set_bc.h"
#include "all_fun.h"

#include <cmath>
#include <iostream>
#include <iomanip>

//----------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------
PartitionedFSI::PartitionedFSI(Simulation* simulation, Integrator* integrator,
                               const PartitionedFSIConfig& config)
  : simulation_(simulation), integrator_(integrator), config_(config),
    omega_(config.initial_relaxation)
{
  resolve_faces();
}

//----------------------------------------------------------------------
// resolve_faces
//----------------------------------------------------------------------
void PartitionedFSI::resolve_faces()
{
  auto& com_mod = simulation_->com_mod;

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    for (int iFa = 0; iFa < com_mod.msh[iM].nFa; iFa++) {
      auto& fa = com_mod.msh[iM].fa[iFa];
      if (fa.name == config_.fluid_interface_face) {
        fluid_face_ = &fa;
        fluid_mesh_ = &com_mod.msh[iM];
      }
      if (fa.name == config_.solid_interface_face) {
        solid_face_ = &fa;
      }
    }
  }

  if (!fluid_face_) {
    throw std::runtime_error("[PartitionedFSI] Fluid interface face '"
        + config_.fluid_interface_face + "' not found.");
  }
  if (!solid_face_) {
    throw std::runtime_error("[PartitionedFSI] Solid interface face '"
        + config_.solid_interface_face + "' not found.");
  }
  if (!fluid_mesh_) {
    throw std::runtime_error("[PartitionedFSI] Fluid mesh not found.");
  }

  // Auto-detect equation indices from equation types
  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    auto phys = com_mod.eq[iEq].phys;
    if (phys == consts::EquationType::phys_fluid) {
      config_.fluid_eq_index = iEq;
    } else if (phys == consts::EquationType::phys_struct ||
               phys == consts::EquationType::phys_ustruct) {
      config_.solid_eq_index = iEq;
    } else if (phys == consts::EquationType::phys_mesh) {
      config_.mesh_eq_index = iEq;
    }
  }

  if (config_.fluid_eq_index < 0) {
    throw std::runtime_error("[PartitionedFSI] No fluid equation found.");
  }
  if (config_.solid_eq_index < 0) {
    throw std::runtime_error("[PartitionedFSI] No struct/ustruct equation found.");
  }
  if (config_.mesh_eq_index < 0) {
    throw std::runtime_error("[PartitionedFSI] No mesh equation found.");
  }
}

//----------------------------------------------------------------------
// compute_aitken_omega
//----------------------------------------------------------------------
double PartitionedFSI::compute_aitken_omega(
    const Array<double>& residual,
    const Array<double>& residual_prev,
    double omega_prev)
{
  // Aitken delta-squared: omega_{k+1} = -omega_k * (r_{k-1} . (r_k - r_{k-1})) / ||r_k - r_{k-1}||^2
  double num = 0.0;
  double den = 0.0;
  for (int a = 0; a < residual.ncols(); a++) {
    for (int i = 0; i < residual.nrows(); i++) {
      double dr = residual(i, a) - residual_prev(i, a);
      num += residual_prev(i, a) * dr;
      den += dr * dr;
    }
  }

  if (den < 1e-30) return omega_prev;

  double omega = -omega_prev * num / den;

  // Clamp to reasonable range
  omega = std::max(0.01, std::min(2.0, std::abs(omega)));

  return omega;
}

//----------------------------------------------------------------------
// step
//----------------------------------------------------------------------
bool PartitionedFSI::step()
{
  auto& com_mod = simulation_->com_mod;
  auto& cm_mod = simulation_->cm_mod;
  auto& cm = com_mod.cm;
  const int nsd = com_mod.nsd;

  auto& solutions = integrator_->get_solutions();
  auto& solid_eq = com_mod.eq[config_.solid_eq_index];
  auto& fluid_eq = com_mod.eq[config_.fluid_eq_index];

  // Initialize omega for this time step
  omega_ = config_.initial_relaxation;

  // Get initial displacement from predictor
  auto disp_current = fsi_coupling::extract_solid_displacement(
      com_mod, solid_eq, *solid_face_, solutions);
  disp_prev_ = disp_current;

  bool converged = false;

  for (int outer = 0; outer < config_.max_coupling_iterations; outer++) {

    // 1. Transfer solid displacement to fluid mesh interface
    auto mesh_disp = fsi_coupling::transfer_face_data(
        com_mod, *solid_face_, *fluid_face_, disp_prev_);

    // 2. Apply displacement as Dirichlet BC on mesh interface
    fsi_coupling::apply_displacement_on_mesh(
        com_mod, com_mod.eq[config_.mesh_eq_index],
        *fluid_face_, mesh_disp, solutions);

    // Apply strong Dirichlet BCs (needed to enforce interface displacement)
    set_bc::set_bc_dir(com_mod, solutions);

    // 3. Solve mesh equation with interface displacement enforced as Dirichlet BC
    integrator_->step_equation(config_.mesh_eq_index,
        [&]() {
          // Enforce Dirichlet BC at interface: zero R and diagonalize Val
          // so the Newton correction is zero at prescribed displacement nodes
          fsi_coupling::enforce_dirichlet_on_face(com_mod, *fluid_face_, nsd);
        });

    // 4. Solve fluid equation to convergence
    integrator_->step_equation(config_.fluid_eq_index);

    // 5. Extract fluid traction at interface
    auto fluid_traction = fsi_coupling::extract_fluid_traction(
        com_mod, cm_mod, *fluid_mesh_, *fluid_face_, fluid_eq,
        integrator_->get_Yg(), integrator_->get_Dg(), solutions);

    // Transfer traction from fluid face to solid face
    auto solid_traction = fsi_coupling::transfer_face_data(
        com_mod, *fluid_face_, *solid_face_, fluid_traction);

    // 6. Solve solid equation with traction as Neumann BC
    integrator_->step_equation(config_.solid_eq_index,
        [&]() {
          fsi_coupling::apply_traction_on_solid(
              com_mod, solid_eq, *solid_face_, solid_traction);
        });

    // 7. Extract new solid displacement at interface
    disp_current = fsi_coupling::extract_solid_displacement(
        com_mod, solid_eq, *solid_face_, solutions);

    // 8. Compute displacement residual
    Array<double> disp_residual(nsd, solid_face_->nNo);
    for (int a = 0; a < solid_face_->nNo; a++) {
      for (int i = 0; i < nsd; i++) {
        disp_residual(i, a) = disp_current(i, a) - disp_prev_(i, a);
      }
    }

    // 9. Aitken relaxation
    if (outer > 0 && config_.use_aitken) {
      omega_ = compute_aitken_omega(disp_residual, disp_residual_prev_, omega_);
    }

    // Apply relaxation: d_prev = d_prev + omega * (d_current - d_prev)
    for (int a = 0; a < solid_face_->nNo; a++) {
      for (int i = 0; i < nsd; i++) {
        disp_prev_(i, a) += omega_ * disp_residual(i, a);
      }
    }

    // Store residual for next Aitken iteration
    disp_residual_prev_ = disp_residual;

    // 10. Check convergence
    double res_norm = 0.0;
    double disp_norm = 0.0;
    for (int a = 0; a < solid_face_->nNo; a++) {
      for (int i = 0; i < nsd; i++) {
        res_norm += disp_residual(i, a) * disp_residual(i, a);
        disp_norm += disp_current(i, a) * disp_current(i, a);
      }
    }
    res_norm = sqrt(res_norm);
    disp_norm = sqrt(disp_norm);
    double rel_change = (disp_norm > 1e-30) ? res_norm / disp_norm : res_norm;

    // Print coupling iteration info (matching solver output format)
    if (cm.mas(cm_mod)) {
      int dB = (rel_change > 0 && disp_norm > 0)
               ? static_cast<int>(10.0 * log10(rel_change)) : 0;
      char buf[128];
      snprintf(buf, sizeof(buf),
               " CP %d-%d  %4.3e  [%d %.3e %.3e]",
               com_mod.cTS, outer + 1,
               com_mod.timer.get_elapsed_time(),
               dB, rel_change, omega_);
      std::cout << buf << std::endl;
    }

    if (rel_change < config_.coupling_tolerance) {
      converged = true;
      break;
    }
  }

  return converged;
}
