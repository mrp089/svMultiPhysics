// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "PartitionedFSI.h"
#include "fsi_coupling.h"
#include "set_bc.h"
#include "all_fun.h"
#include "distribute.h"
#include "initialize.h"
#include "output.h"
#include "vtk_xml.h"
#include "txt.h"
#include "read_files.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>

// Forward declaration of add_eq_linear_algebra (defined in main.cpp)
void add_eq_linear_algebra(ComMod& com_mod, eqType& lEq);

/// Check if any value in the solution arrays is NaN
static bool has_nan(const SolutionStates& sol) {
  const Array<double>* arrays[] = {
    &sol.current.get_velocity(),
    &sol.current.get_acceleration(),
    &sol.current.get_displacement()
  };
  for (auto* arr : arrays)
    for (int a = 0; a < arr->ncols(); a++)
      for (int i = 0; i < arr->nrows(); i++)
        if (std::isnan((*arr)(i, a))) return true;
  return false;
}

//----------------------------------------------------------------------
// Helper: initialize one sub-simulation through the standard pipeline
//----------------------------------------------------------------------
static void init_sub_sim(Simulation* sim, const std::string& xml_path)
{
  read_files_ns::read_files(sim, xml_path);
  distribute(sim);
  Vector<double> init_time(3);
  initialize(sim, init_time);
  for (int iEq = 0; iEq < sim->com_mod.nEq; iEq++) {
    add_eq_linear_algebra(sim->com_mod, sim->com_mod.eq[iEq]);
  }
}

//----------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------
PartitionedFSI::PartitionedFSI(Simulation* main_simulation,
                               const PartitionedFSIConfig& config,
                               const std::string& xml_file_path)
  : main_sim_(main_simulation), config_(config),
    xml_file_path_(xml_file_path), omega_(config.initial_relaxation)
{
  auto& cm = main_sim_->com_mod.cm;
  auto& cm_mod = main_sim_->cm_mod;

  // Resolve XML paths relative to the main XML directory
  std::string dir;
  auto slash = xml_file_path_.find_last_of('/');
  if (slash != std::string::npos) {
    dir = xml_file_path_.substr(0, slash + 1);
  }

  std::string fluid_xml = dir + config_.fluid_xml;
  std::string solid_xml = dir + config_.solid_xml;
  std::string mesh_xml  = dir + config_.mesh_xml;

  if (cm.mas(cm_mod)) {
    std::cout << "[PartitionedFSI] Initializing fluid sub-sim: " << fluid_xml << std::endl;
  }
  fluid_sim_ = std::make_unique<Simulation>();
  init_sub_sim(fluid_sim_.get(), fluid_xml);

  if (cm.mas(cm_mod)) {
    std::cout << "[PartitionedFSI] Initializing solid sub-sim: " << solid_xml << std::endl;
  }
  solid_sim_ = std::make_unique<Simulation>();
  init_sub_sim(solid_sim_.get(), solid_xml);

  if (cm.mas(cm_mod)) {
    std::cout << "[PartitionedFSI] Initializing mesh sub-sim:  " << mesh_xml << std::endl;
  }
  mesh_sim_ = std::make_unique<Simulation>();
  init_sub_sim(mesh_sim_.get(), mesh_xml);

  if (cm.mas(cm_mod)) {
    std::cout << "[PartitionedFSI] Sub-sims ready:"
              << "  fluid=" << fluid_sim_->com_mod.tnNo << "n/" << fluid_sim_->com_mod.eq[0].dof << "dof"
              << "  solid=" << solid_sim_->com_mod.tnNo << "n/" << solid_sim_->com_mod.eq[0].dof << "dof"
              << "  mesh="  << mesh_sim_->com_mod.tnNo  << "n/" << mesh_sim_->com_mod.eq[0].dof  << "dof"
              << std::endl;
  }

  resolve_faces();
  build_node_maps();
}

PartitionedFSI::~PartitionedFSI() {}

//----------------------------------------------------------------------
// resolve_faces
//----------------------------------------------------------------------
void PartitionedFSI::resolve_faces()
{
  auto find_face = [](Simulation* sim, const std::string& face_name,
                      const faceType*& face_out, const mshType*& mesh_out) {
    for (int iM = 0; iM < sim->com_mod.nMsh; iM++) {
      auto& msh = sim->com_mod.msh[iM];
      for (int iFa = 0; iFa < msh.nFa; iFa++) {
        if (msh.fa[iFa].name == face_name) {
          face_out = &msh.fa[iFa];
          mesh_out = &msh;
          return;
        }
      }
    }
    throw std::runtime_error("[PartitionedFSI] Face '" + face_name + "' not found.");
  };

  find_face(fluid_sim_.get(), config_.fluid_interface_face, fluid_face_, fluid_mesh_);
  find_face(solid_sim_.get(), config_.solid_interface_face, solid_face_, solid_mesh_);
  find_face(mesh_sim_.get(),  config_.fluid_interface_face, mesh_face_,  mesh_mesh_);
}

//----------------------------------------------------------------------
// build_face_node_map
//----------------------------------------------------------------------
void PartitionedFSI::build_face_node_map(
    const faceType& face_a, const ComMod& com_a,
    const faceType& face_b, const ComMod& com_b,
    std::vector<int>& a_to_b)
{
  const int nsd = com_a.nsd;
  const double tol = 1e-8;
  a_to_b.assign(face_a.nNo, -1);

  for (int a = 0; a < face_a.nNo; a++) {
    int Ac = face_a.gN(a);
    double best = 1e30;
    int best_b = -1;
    for (int b = 0; b < face_b.nNo; b++) {
      int Bc = face_b.gN(b);
      double d2 = 0.0;
      for (int i = 0; i < nsd; i++) {
        double d = com_a.x(i, Ac) - com_b.x(i, Bc);
        d2 += d * d;
      }
      if (d2 < best) { best = d2; best_b = b; }
    }
    if (best < tol * tol) a_to_b[a] = best_b;
  }
}

//----------------------------------------------------------------------
// build_node_maps
//----------------------------------------------------------------------
void PartitionedFSI::build_node_maps()
{
  build_face_node_map(*solid_face_, solid_sim_->com_mod,
                      *fluid_face_, fluid_sim_->com_mod, solid_to_fluid_map_);
  build_face_node_map(*fluid_face_, fluid_sim_->com_mod,
                      *solid_face_, solid_sim_->com_mod, fluid_to_solid_map_);
  build_face_node_map(*solid_face_, solid_sim_->com_mod,
                      *mesh_face_,  mesh_sim_->com_mod,  solid_to_mesh_map_);

  auto& cm = main_sim_->com_mod.cm;
  auto& cm_mod = main_sim_->cm_mod;
  if (cm.mas(cm_mod)) {
    int matched = 0;
    for (int v : solid_to_fluid_map_) if (v >= 0) matched++;
    std::cout << "[PartitionedFSI] Interface node maps: " << matched
              << "/" << solid_face_->nNo << " matched" << std::endl;
  }
}

//----------------------------------------------------------------------
// transfer_data
//----------------------------------------------------------------------
Array<double> PartitionedFSI::transfer_data(
    const std::vector<int>& src_to_tgt_map,
    const Array<double>& src_data, int tgt_nNo)
{
  int nrows = src_data.nrows();
  Array<double> result(nrows, tgt_nNo);
  for (int a = 0; a < static_cast<int>(src_to_tgt_map.size()); a++) {
    int b = src_to_tgt_map[a];
    if (b >= 0) {
      for (int i = 0; i < nrows; i++) result(i, b) = src_data(i, a);
    }
  }
  return result;
}

//----------------------------------------------------------------------
// compute_aitken_omega
//----------------------------------------------------------------------
double PartitionedFSI::compute_aitken_omega(
    const Array<double>& residual,
    const Array<double>& residual_prev,
    double omega_prev)
{
  double num = 0.0, den = 0.0;
  for (int a = 0; a < residual.ncols(); a++) {
    for (int i = 0; i < residual.nrows(); i++) {
      double dr = residual(i, a) - residual_prev(i, a);
      num += residual_prev(i, a) * dr;
      den += dr * dr;
    }
  }
  if (den < 1e-30) return omega_prev;
  double omega = -omega_prev * num / den;
  return std::max(0.01, std::min(1.0, std::abs(omega)));
}

//======================================================================
// run — full time-stepping loop with Dirichlet-Neumann coupling
//======================================================================
void PartitionedFSI::run()
{
  auto& main_com = main_sim_->com_mod;
  auto& cm_mod = main_sim_->cm_mod;
  auto& cm = main_com.cm;

  int nTS = main_com.nTS;
  int& cTS = main_com.cTS;
  double& dt = main_com.dt;
  double& time = main_com.time;
  int nITs = main_com.nITs;

  if (cTS <= nITs) dt = dt / 10.0;

  Simulation* sims[3] = {fluid_sim_.get(), solid_sim_.get(), mesh_sim_.get()};

  while (true) {
    if (cTS == nITs) dt = 10.0 * dt;
    cTS = cTS + 1;
    time = time + dt;

    // Sync time to sub-sims
    for (auto* sim : sims) {
      sim->com_mod.cTS = cTS;
      sim->com_mod.time = time;
      sim->com_mod.dt = dt;
      for (auto& eq : sim->com_mod.eq) { eq.itr = 0; eq.ok = false; }
    }

    if (cm.mas(cm_mod)) {
      std::cout << std::string(70, '=') << std::endl;
      std::cout << "  TIME STEP " << cTS << "  t=" << time << "  dt=" << dt << std::endl;
      std::cout << std::string(70, '=') << std::endl;
    }

    // Predictor + Dirichlet BCs for each sub-sim
    for (auto* sim : sims) {
      sim->get_integrator().predictor();
      set_bc::set_bc_dir(sim->com_mod, sim->get_integrator().get_solutions());
    }

    // Coupling loop
    bool converged = step();

    if (cm.mas(cm_mod)) {
      if (converged) {
        std::cout << "  TIME STEP " << cTS << " CONVERGED" << std::endl;
      } else {
        std::cout << "  TIME STEP " << cTS << " FAILED (NaN or no convergence)" << std::endl;
      }
    }

    // Stop on failure
    if (!converged) break;

    // Save results
    save_results();

    // Copy current -> old
    for (auto* sim : sims) {
      auto& sol = sim->get_integrator().get_solutions();
      sol.old.get_acceleration() = sol.current.get_acceleration();
      sol.old.get_velocity()     = sol.current.get_velocity();
      if (sim->com_mod.dFlag)
        sol.old.get_displacement() = sol.current.get_displacement();
    }

    // Stop condition
    int stopTS = nTS;
    if (cm.mas(cm_mod)) {
      if (FILE* fp = fopen(main_com.stopTrigName.c_str(), "r")) {
        int count = fscanf(fp, "%d", &stopTS);
        if (count == 0) stopTS = cTS;
        fclose(fp);
      }
    }
    cm.bcast(cm_mod, &stopTS);
    if (cTS >= stopTS) break;
  }
}

//======================================================================
// step — one time step of Dirichlet-Neumann coupling with Aitken
//======================================================================
bool PartitionedFSI::step()
{
  auto& fluid_com = fluid_sim_->com_mod;
  auto& solid_com = solid_sim_->com_mod;
  auto& mesh_com  = mesh_sim_->com_mod;
  auto& cm_mod = main_sim_->cm_mod;
  auto& cm = main_sim_->com_mod.cm;
  const int nsd = main_sim_->com_mod.nsd;
  const int cTS = main_sim_->com_mod.cTS;

  auto& fluid_int = fluid_sim_->get_integrator();
  auto& solid_int = solid_sim_->get_integrator();
  auto& mesh_int  = mesh_sim_->get_integrator();
  auto& fluid_sol = fluid_int.get_solutions();
  auto& solid_sol = solid_int.get_solutions();
  auto& mesh_sol  = mesh_int.get_solutions();

  omega_ = config_.initial_relaxation;

  // Save predictor state for each sub-sim. Each coupling iteration must
  // start from the same initial guess (post-predictor), otherwise Newton
  // corrections accumulate and the coupling diverges.
  struct SavedState {
    Array<double> An, Yn, Dn;
  };
  auto save_state = [](SolutionStates& sol) -> SavedState {
    return {sol.current.get_acceleration(),
            sol.current.get_velocity(),
            sol.current.get_displacement()};
  };
  auto restore_state = [](SolutionStates& sol, const SavedState& s) {
    sol.current.get_acceleration() = s.An;
    sol.current.get_velocity() = s.Yn;
    sol.current.get_displacement() = s.Dn;
  };

  SavedState fluid_pred = save_state(fluid_sol);
  SavedState solid_pred = save_state(solid_sol);
  SavedState mesh_pred  = save_state(mesh_sol);

  // Initial displacement and velocity from predictor
  auto disp_current = fsi_coupling::extract_solid_displacement(
      solid_com, solid_com.eq[0], *solid_face_, solid_sol);
  auto vel_current = fsi_coupling::extract_solid_velocity(
      solid_com, solid_com.eq[0], *solid_face_, solid_sol);
  disp_prev_ = disp_current;
  vel_prev_ = vel_current;

  bool converged = false;

  for (int cp = 0; cp < config_.max_coupling_iterations; cp++) {

    // Restore all sub-sims to predictor state
    restore_state(fluid_sol, fluid_pred);
    restore_state(solid_sol, solid_pred);
    restore_state(mesh_sol, mesh_pred);

    if (cm.mas(cm_mod)) {
      std::cout << std::string(50, '-') << std::endl;
      std::cout << "  COUPLING ITERATION " << cp + 1
                << " (time step " << cTS << ")" << std::endl;
      std::cout << std::string(50, '-') << std::endl;
    }

    // ---- Use relaxed displacement and velocity at interface ----
    auto mesh_disp = transfer_data(solid_to_mesh_map_, disp_prev_, mesh_face_->nNo);
    auto fluid_vel = transfer_data(solid_to_fluid_map_, vel_prev_, fluid_face_->nNo);

    // Debug: print max interface values
    if (cm.mas(cm_mod)) {
      double max_disp = 0, max_vel = 0;
      for (int a = 0; a < mesh_disp.ncols(); a++)
        for (int i = 0; i < nsd; i++) {
          max_disp = std::max(max_disp, std::abs(mesh_disp(i, a)));
          max_vel  = std::max(max_vel,  std::abs(fluid_vel(i, a)));
        }
      std::cout << "    interface: max|disp|=" << max_disp
                << "  max|vel|=" << max_vel << std::endl;
    }

    // ---- MESH SOLVE ----
    if (cm.mas(cm_mod)) {
      std::cout << "    [mesh] solving..." << std::endl;
    }
    // Apply XML BCs first, then overwrite interface with coupling values
    set_bc::set_bc_dir(mesh_com, mesh_sol);
    fsi_coupling::apply_displacement_on_mesh(
        mesh_com, mesh_com.eq[0], *mesh_face_, mesh_disp, mesh_sol);
    mesh_int.step_equation(0, [&]() {
      fsi_coupling::enforce_dirichlet_on_face(mesh_com, *mesh_face_, nsd);
    });
    if (has_nan(mesh_sol)) {
      if (cm.mas(cm_mod)) std::cout << "  ABORT: NaN in mesh solve" << std::endl;
      return false;
    }

    // ---- ALE: deform fluid mesh ----
    auto& mesh_Dg = mesh_int.get_Dg();
    for (int a = 0; a < fluid_com.tnNo; a++)
      for (int i = 0; i < nsd; i++)
        fluid_com.x(i, a) += mesh_Dg(i, a);

    // ---- FLUID SOLVE ----
    if (cm.mas(cm_mod)) {
      std::cout << "    [fluid] solving..." << std::endl;
    }
    // Apply XML BCs first, then overwrite interface with coupling values
    set_bc::set_bc_dir(fluid_com, fluid_sol);
    fsi_coupling::apply_velocity_on_fluid(
        fluid_com, fluid_com.eq[0], *fluid_face_, fluid_vel, fluid_sol);
    fluid_int.step_equation(0, [&]() {
      fsi_coupling::enforce_dirichlet_dofs_on_face(fluid_com, *fluid_face_, 0, nsd);
    });
    if (has_nan(fluid_sol)) {
      if (cm.mas(cm_mod)) std::cout << "  ABORT: NaN in fluid solve" << std::endl;
      for (int a = 0; a < fluid_com.tnNo; a++)
        for (int i = 0; i < nsd; i++)
          fluid_com.x(i, a) -= mesh_Dg(i, a);
      return false;
    }

    // ---- Extract traction ON DEFORMED MESH ----
    // Must extract before restoring coordinates: the velocity/pressure
    // solution was computed on the deformed mesh, so the element geometry
    // used for velocity gradients must match.
    auto fluid_traction = fsi_coupling::extract_fluid_traction(
        fluid_com, fluid_sim_->cm_mod,
        *fluid_mesh_, *fluid_face_, fluid_com.eq[0],
        fluid_int.get_Yg(), fluid_int.get_Dg(), fluid_sol);
    auto solid_traction = transfer_data(fluid_to_solid_map_,
                                        fluid_traction, solid_face_->nNo);

    // ---- Restore fluid mesh coordinates ----
    for (int a = 0; a < fluid_com.tnNo; a++)
      for (int i = 0; i < nsd; i++)
        fluid_com.x(i, a) -= mesh_Dg(i, a);

    // ---- SOLID SOLVE ----
    if (cm.mas(cm_mod)) {
      double trac_max = 0.0;
      for (int a = 0; a < solid_face_->nNo; a++)
        for (int i = 0; i < nsd; i++)
          trac_max = std::max(trac_max, std::abs(solid_traction(i, a)));
      std::cout << "    [solid] solving (max traction=" << trac_max << ")..." << std::endl;
    }
    // Re-apply XML BCs (inlet/outlet Dir) before solid solve
    set_bc::set_bc_dir(solid_com, solid_sol);
    solid_int.step_equation(0, [&]() {
      fsi_coupling::apply_traction_on_solid(
          solid_com, solid_com.eq[0], *solid_face_, solid_traction);
    });
    if (has_nan(solid_sol)) {
      if (cm.mas(cm_mod)) std::cout << "  ABORT: NaN in solid solve" << std::endl;
      return false;
    }

    // ---- Extract new solid displacement and velocity ----
    disp_current = fsi_coupling::extract_solid_displacement(
        solid_com, solid_com.eq[0], *solid_face_, solid_sol);
    vel_current = fsi_coupling::extract_solid_velocity(
        solid_com, solid_com.eq[0], *solid_face_, solid_sol);

    // ---- Residual + Aitken relaxation (applied to both disp and vel) ----
    Array<double> disp_residual(nsd, solid_face_->nNo);
    for (int a = 0; a < solid_face_->nNo; a++)
      for (int i = 0; i < nsd; i++)
        disp_residual(i, a) = disp_current(i, a) - disp_prev_(i, a);

    if (cp > 0 && config_.use_aitken)
      omega_ = compute_aitken_omega(disp_residual, disp_residual_prev_, omega_);

    // Relax displacement and velocity with the same omega
    for (int a = 0; a < solid_face_->nNo; a++)
      for (int i = 0; i < nsd; i++) {
        disp_prev_(i, a) += omega_ * (disp_current(i, a) - disp_prev_(i, a));
        vel_prev_(i, a)  += omega_ * (vel_current(i, a)  - vel_prev_(i, a));
      }
    disp_residual_prev_ = disp_residual;

    // ---- Convergence check ----
    double res_norm = 0.0, disp_norm = 0.0;
    for (int a = 0; a < solid_face_->nNo; a++)
      for (int i = 0; i < nsd; i++) {
        res_norm  += disp_residual(i, a) * disp_residual(i, a);
        disp_norm += disp_current(i, a) * disp_current(i, a);
      }
    res_norm  = sqrt(res_norm);
    disp_norm = sqrt(disp_norm);
    double rel = (disp_norm > 1e-30) ? res_norm / disp_norm : res_norm;

    if (cm.mas(cm_mod)) {
      std::cout << "  COUPLING " << cTS << "-" << cp + 1
                << "  rel_change=" << rel
                << "  omega=" << omega_
                << "  |disp|=" << disp_norm
                << std::endl;
    }

    // Check for NaN — abort early
    if (std::isnan(rel) || std::isnan(disp_norm)) {
      if (cm.mas(cm_mod)) {
        std::cout << "  COUPLING " << cTS << "-" << cp + 1
                  << "  ABORT: NaN detected" << std::endl;
      }
      return false;
    }

    if (rel < config_.coupling_tolerance) { converged = true; break; }
  }
  return converged;
}

//----------------------------------------------------------------------
// save_results
//----------------------------------------------------------------------
void PartitionedFSI::save_results()
{
  int cTS = main_sim_->com_mod.cTS;
  Simulation* sims[3] = {fluid_sim_.get(), solid_sim_.get(), mesh_sim_.get()};

  for (auto* sim : sims) {
    auto& com = sim->com_mod;
    auto& sol = sim->get_integrator().get_solutions();
    if (com.saveVTK) {
      bool l2 = ((cTS % com.saveIncr) == 0);
      bool l3 = (cTS >= com.saveATS);
      if (l2 && l3) {
        output::output_result(sim, com.timeP, 3, 0);
        vtk_xml::write_vtus(sim, sol, false);
      }
    }
  }
}

//----------------------------------------------------------------------
// Stubs for unused methods declared in the header
//----------------------------------------------------------------------
void PartitionedFSI::create_sub_simulations() {}
std::string PartitionedFSI::generate_sub_xml(const std::string&) { return ""; }
void PartitionedFSI::init_sub_simulation(Simulation*, const std::string&) {}
