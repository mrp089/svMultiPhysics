// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "PartitionedFSI.h"
#include "Integrator.h"
#include "fsi_coupling.h"
#include "post.h"
#include "set_bc.h"
#include "distribute.h"
#include "initialize.h"
#include "output.h"
#include "vtk_xml.h"
#include "read_files.h"

#include <cmath>
#include <iostream>
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
  sim->logger.set_cout_write(false);
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

  // 3 separate sub-sims: fluid, solid, mesh
  fluid_sim_ = std::make_unique<Simulation>();
  init_sub_sim(fluid_sim_.get(), fluid_xml);

  solid_sim_ = std::make_unique<Simulation>();
  init_sub_sim(solid_sim_.get(), solid_xml);

  mesh_sim_ = std::make_unique<Simulation>();
  init_sub_sim(mesh_sim_.get(), mesh_xml);


  if (cm.mas(cm_mod)) {
    // Open log files
    std::string log_dir = fluid_sim_->get_chnl_mod().appPath;
    coupling_log_.open(log_dir + "coupling.dat");
    char hdr[256];
    snprintf(hdr, sizeof(hdr), "# %4s %3s %10s %5s %10s %10s %10s %10s",
             "cTS", "cp", "time", "dB", "Ri/R1", "Ri/R0", "omega", "|disp|");
    coupling_log_ << hdr << std::endl;

    histor_log_.open(log_dir + "histor.dat");
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
// relax_interface — updates disp_prev_ and vel_prev_
//----------------------------------------------------------------------
void PartitionedFSI::relax_interface(int cp, int nsd,
                                     const Array<double>& disp_current)
{
  switch (config_.coupling_method) {
    case CouplingMethod::constant:
      relax_constant(cp, nsd, disp_current);
      break;
    case CouplingMethod::aitken:
      relax_aitken(cp, nsd, disp_current);
      break;
  }
}

//----------------------------------------------------------------------
// relax_constant — fixed relaxation
//----------------------------------------------------------------------
void PartitionedFSI::relax_constant(int cp, int nsd,
                                     const Array<double>& disp_current)
{
  omega_ = config_.initial_relaxation;
  for (int a = 0; a < solid_face_->nNo; a++)
    for (int i = 0; i < nsd; i++)
      disp_prev_(i, a) += omega_ * (disp_current(i, a) - disp_prev_(i, a));
}

//----------------------------------------------------------------------
// relax_aitken — Aitken Delta^2 (Küttler & Wall 2008, Eq. 44)
//----------------------------------------------------------------------
void PartitionedFSI::relax_aitken(int cp, int nsd,
                                   const Array<double>& disp_current)
{
  const int u = nsd * solid_face_->nNo;

  // Build residual r = x_tilde - x
  std::vector<double> r(u);
  for (int a = 0; a < solid_face_->nNo; a++)
    for (int i = 0; i < nsd; i++)
      r[a * nsd + i] = disp_current(i, a) - disp_prev_(i, a);

  // Aitken update: omega = -omega * r^T (r_new - r_old) / |r_new - r_old|^2
  // Negative omega allowed (corrects overshoot)
  if (cp > 0 && !r_prev_.empty()) {
    double num = 0, den = 0;
    for (int j = 0; j < u; j++) {
      double dr = r[j] - r_prev_[j];
      num += r_prev_[j] * dr;
      den += dr * dr;
    }
    if (den > 1e-30) {
      omega_ = -omega_ * num / den;
      if (std::abs(omega_) > config_.omega_max)
        omega_ = (omega_ > 0) ? config_.omega_max : -config_.omega_max;
    }
  }
  r_prev_ = r;

  // Apply: x_{k+1} = x_k + omega * r
  for (int a = 0; a < solid_face_->nNo; a++)
    for (int i = 0; i < nsd; i++)
      disp_prev_(i, a) += omega_ * (disp_current(i, a) - disp_prev_(i, a));
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
      if (histor_log_.is_open()) {
        histor_log_ << std::string(70, '=') << std::endl;
        histor_log_ << "  TIME STEP " << cTS << "  t=" << time << "  dt=" << dt << std::endl;
        histor_log_ << std::string(70, '=') << std::endl;
      }
    }

    // Predictor + Dirichlet BCs for each sub-sim
    for (auto* sim : sims) {
      sim->get_integrator().predictor();
      set_bc::set_bc_dir(sim->com_mod, sim->get_integrator().get_solutions());
    }

    // Coupling loop
    bool converged = step();

    if (!converged && cm.mas(cm_mod)) {
      std::cout << "  TIME STEP " << cTS << " FAILED (NaN or no convergence)" << std::endl;
      if (histor_log_.is_open())
        histor_log_ << "  TIME STEP " << cTS << " FAILED (NaN or no convergence)" << std::endl;
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

//----------------------------------------------------------------------
// compute_interface_velocity — Newmark-consistent velocity from disp_prev_
//----------------------------------------------------------------------
void PartitionedFSI::compute_interface_velocity()
{
  auto& solid_com = solid_sim_->com_mod;
  auto& solid_sol = solid_sim_->get_integrator().get_solutions();
  const auto& eq = solid_com.eq[0];
  const int s = eq.s;
  const int nsd = main_sim_->com_mod.nsd;
  const double dt = solid_com.dt;
  const auto& Do = solid_sol.old.get_displacement();
  const auto& Yo = solid_sol.old.get_velocity();
  const auto& Ao = solid_sol.old.get_acceleration();

  vel_prev_.resize(nsd, solid_face_->nNo);
  for (int a = 0; a < solid_face_->nNo; a++) {
    int Ac = solid_face_->gN(a);
    for (int i = 0; i < nsd; i++) {
      double a_new, v_new;
      newmark::state_from_displacement(
          disp_prev_(i, a), Do(i + s, Ac), Yo(i + s, Ac), Ao(i + s, Ac),
          dt, eq.beta, eq.gam, a_new, v_new);
      vel_prev_(i, a) = v_new;
    }
  }
}

//----------------------------------------------------------------------
// solve_fluid — fluid equation with interface velocity and ALE
//----------------------------------------------------------------------
bool PartitionedFSI::solve_fluid(
    const Array<double>& mesh_vel_Yo, const Array<double>& mesh_vel_Yn)
{
  auto& fluid_com = fluid_sim_->com_mod;
  auto& fluid_int = fluid_sim_->get_integrator();
  auto& fluid_sol = fluid_int.get_solutions();
  const int nsd = main_sim_->com_mod.nsd;

  auto fluid_vel = transfer_data(solid_to_fluid_map_, vel_prev_, fluid_face_->nNo);
  set_bc::set_bc_dir(fluid_com, fluid_sol);
  fsi_coupling::apply_velocity_on_fluid(
      fluid_com, fluid_com.eq[0], *fluid_face_, fluid_vel, fluid_sol);

  // ALE mesh velocity at generalized-alpha intermediate time
  double af = fluid_com.eq[0].af;
  fluid_com.ale_mesh_velocity.resize(nsd, fluid_com.tnNo);
  for (int a = 0; a < fluid_com.tnNo; a++)
    for (int i = 0; i < nsd; i++)
      fluid_com.ale_mesh_velocity(i, a) = (1.0 - af) * mesh_vel_Yo(i, a)
                                         + af * mesh_vel_Yn(i, a);

  fluid_int.step_equation(0, [&]() {
    set_bc::enforce_dirichlet_dofs_on_face(fluid_com, *fluid_face_, 0, nsd);
  });
  return !has_nan(fluid_sol);
}

//----------------------------------------------------------------------
// solve_solid — extract traction from fluid, solve solid
//----------------------------------------------------------------------
bool PartitionedFSI::solve_solid()
{
  auto& fluid_com = fluid_sim_->com_mod;
  auto& solid_com = solid_sim_->com_mod;
  auto& fluid_int = fluid_sim_->get_integrator();
  auto& solid_int = solid_sim_->get_integrator();
  auto& solid_sol = solid_int.get_solutions();

  auto fluid_traction = post::compute_face_traction(
      fluid_com, fluid_sim_->cm_mod,
      *fluid_mesh_, *fluid_face_, fluid_com.eq[0],
      fluid_int.get_Yg(), fluid_int.get_Dg(), fluid_int.get_solutions());
  auto solid_traction = transfer_data(fluid_to_solid_map_,
                                      fluid_traction, solid_face_->nNo);

  set_bc::set_bc_dir(solid_com, solid_sol);
  solid_int.step_equation(0, [&]() {
    fsi_coupling::apply_traction_on_solid(
        solid_com, solid_com.eq[0], *solid_face_, solid_traction);
  });
  return !has_nan(solid_sol);
}

//----------------------------------------------------------------------
// solve_mesh — mesh equation with relaxed displacement, deform fluid mesh
//----------------------------------------------------------------------
bool PartitionedFSI::solve_mesh(const Array<double>& x_ref, int mesh_s)
{
  auto& fluid_com = fluid_sim_->com_mod;
  auto& mesh_com  = mesh_sim_->com_mod;
  auto& mesh_int  = mesh_sim_->get_integrator();
  auto& mesh_sol  = mesh_int.get_solutions();
  const int nsd = main_sim_->com_mod.nsd;

  auto mesh_disp = transfer_data(solid_to_mesh_map_, disp_prev_, mesh_face_->nNo);
  set_bc::set_bc_dir(mesh_com, mesh_sol);
  fsi_coupling::apply_displacement_on_mesh(
      mesh_com, mesh_com.eq[0], *mesh_face_, mesh_disp, mesh_sol);
  mesh_int.step_equation(0, [&]() {
    set_bc::enforce_dirichlet_on_face(mesh_com, *mesh_face_, nsd);
  });
  if (has_nan(mesh_sol)) return false;

  // Deform fluid mesh: apply only the INCREMENT (Dn - Do) to x_ref
  // so that fluid_com.x = x_original + Dn
  auto& mesh_Dn = mesh_sol.current.get_displacement();
  auto& mesh_Do = mesh_sol.old.get_displacement();
  for (int a = 0; a < fluid_com.tnNo; a++)
    for (int i = 0; i < nsd; i++)
      fluid_com.x(i, a) = x_ref(i, a)
                         + mesh_Dn(i + mesh_s, a) - mesh_Do(i + mesh_s, a);
  return true;
}

//======================================================================
// step — one coupling iteration loop for one time step
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

  auto& fluid_sol = fluid_sim_->get_integrator().get_solutions();
  auto& solid_sol = solid_sim_->get_integrator().get_solutions();
  auto& mesh_sol  = mesh_sim_->get_integrator().get_solutions();

  omega_ = config_.initial_relaxation;
  r_prev_.clear();

  // Save predictor state
  struct SavedState { Array<double> An, Yn, Dn; };
  auto save_state = [](SolutionStates& s) -> SavedState {
    return {s.current.get_acceleration(), s.current.get_velocity(), s.current.get_displacement()};
  };
  auto restore_state = [](SolutionStates& s, const SavedState& st) {
    s.current.get_acceleration() = st.An;
    s.current.get_velocity() = st.Yn;
    s.current.get_displacement() = st.Dn;
  };
  SavedState fluid_pred = save_state(fluid_sol);
  SavedState solid_pred = save_state(solid_sol);
  SavedState mesh_pred  = save_state(mesh_sol);

  // Save mesh coordinates at start of time step = x_original + Do
  Array<double> x_ref(fluid_com.x);

  // ALE mesh velocity from predictor (updated after each mesh solve)
  const int mesh_s = mesh_com.eq[0].s;
  Array<double> mesh_vel_Yn(nsd, mesh_com.tnNo);
  Array<double> mesh_vel_Yo(nsd, mesh_com.tnNo);
  {
    auto& mYn = mesh_sol.current.get_velocity();
    auto& mYo = mesh_sol.old.get_velocity();
    for (int a = 0; a < mesh_com.tnNo; a++)
      for (int i = 0; i < nsd; i++) {
        mesh_vel_Yn(i, a) = mYn(mesh_s + i, a);
        mesh_vel_Yo(i, a) = mYo(mesh_s + i, a);
      }
  }

  // Initial interface state from predictor
  auto disp_current = fsi_coupling::extract_solid_displacement(
      solid_com, solid_com.eq[0], *solid_face_, solid_sol);
  disp_prev_ = disp_current;
  compute_interface_velocity();

  bool converged = false;

  for (int cp = 0; cp < config_.max_coupling_iterations; cp++) {

    // Restore all sub-sims to predictor state
    restore_state(fluid_sol, fluid_pred);
    restore_state(solid_sol, solid_pred);
    restore_state(mesh_sol, mesh_pred);
    fluid_com.x = x_ref;

    // ---- 1. Fluid solve ----
    if (!solve_fluid(mesh_vel_Yo, mesh_vel_Yn)) {
      if (cm.mas(cm_mod)) std::cout << "  ABORT: NaN in fluid solve" << std::endl;
      return false;
    }

    // ---- 2. Solid solve ----
    if (!solve_solid()) {
      if (cm.mas(cm_mod)) std::cout << "  ABORT: NaN in solid solve" << std::endl;
      return false;
    }

    // ---- 3. Extract displacement, check convergence ----
    disp_current = fsi_coupling::extract_solid_displacement(
        solid_com, solid_com.eq[0], *solid_face_, solid_sol);

    double res_norm = 0.0, disp_norm = 0.0;
    for (int a = 0; a < solid_face_->nNo; a++)
      for (int i = 0; i < nsd; i++) {
        double res = disp_current(i, a) - disp_prev_(i, a);
        res_norm  += res * res;
        disp_norm += disp_current(i, a) * disp_current(i, a);
      }
    res_norm  = sqrt(res_norm);
    disp_norm = sqrt(disp_norm);
    double rel = (disp_norm > 1e-30) ? res_norm / disp_norm : res_norm;

    // ---- 4. Relaxation ----
    relax_interface(cp, nsd, disp_current);
    compute_interface_velocity();

    // Check for NaN/divergence
    {
      bool bad = false;
      double max_disp = 0;
      for (int a = 0; a < solid_face_->nNo && !bad; a++)
        for (int i = 0; i < nsd; i++) {
          if (std::isnan(disp_prev_(i, a)) || std::isinf(disp_prev_(i, a)))
            { bad = true; break; }
          max_disp = std::max(max_disp, std::abs(disp_prev_(i, a)));
        }
      if (bad || max_disp > 1e10) {
        if (cm.mas(cm_mod)) std::cout << "  ABORT: NaN/divergence after relaxation" << std::endl;
        return false;
      }
    }

    // ---- 5. Mesh solve + deform fluid mesh ----
    if (!solve_mesh(x_ref, mesh_s)) {
      if (cm.mas(cm_mod)) std::cout << "  ABORT: NaN in mesh solve" << std::endl;
      return false;
    }

    // Update ALE mesh velocity for next coupling iteration
    {
      auto& mYn = mesh_sol.current.get_velocity();
      for (int a = 0; a < mesh_com.tnNo; a++)
        for (int i = 0; i < nsd; i++)
          mesh_vel_Yn(i, a) = mYn(mesh_s + i, a);
    }

    // ---- 6. Output ----
    if (cp == 0) first_res_norm_ = res_norm;
    int dB_val = 0;
    double ri_r1 = 1.0;
    if (first_res_norm_ > 1e-30 && res_norm > 0) {
      ri_r1 = res_norm / first_res_norm_;
      dB_val = static_cast<int>(20.0 * log10(ri_r1));
    }

    if (cm.mas(cm_mod)) {
      bool conv = rel < config_.coupling_tolerance;
      bool saved = conv
                && (cTS % fluid_sim_->com_mod.saveIncr == 0)
                && (cTS >= fluid_sim_->com_mod.saveATS);
      char buf[256];
      snprintf(buf, sizeof(buf), " CP %d-%d%s %10.3e %5d %10.3e %10.3e %10.3e %10.3e",
               cTS, cp + 1, saved ? "s" : " ",
               main_sim_->com_mod.timer.get_elapsed_time(),
               dB_val, ri_r1, rel, omega_, disp_norm);
      std::cout << buf << std::endl;
      if (coupling_log_.is_open()) coupling_log_ << buf << std::endl;
      if (histor_log_.is_open()) histor_log_ << buf << std::endl;
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

