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

  // 3 separate sub-sims: fluid, solid, mesh
  if (cm.mas(cm_mod)) std::cout << "[PartitionedFSI] Initializing fluid: " << fluid_xml << std::endl;
  fluid_sim_ = std::make_unique<Simulation>();
  init_sub_sim(fluid_sim_.get(), fluid_xml);

  if (cm.mas(cm_mod)) std::cout << "[PartitionedFSI] Initializing solid: " << solid_xml << std::endl;
  solid_sim_ = std::make_unique<Simulation>();
  init_sub_sim(solid_sim_.get(), solid_xml);

  if (cm.mas(cm_mod)) std::cout << "[PartitionedFSI] Initializing mesh:  " << mesh_xml << std::endl;
  mesh_sim_ = std::make_unique<Simulation>();
  init_sub_sim(mesh_sim_.get(), mesh_xml);

  if (cm.mas(cm_mod)) {
    std::cout << "[PartitionedFSI] Sub-sims ready:"
              << "  fluid=" << fluid_sim_->com_mod.tnNo << "n/" << fluid_sim_->com_mod.tDof << "tDof"
              << "  solid=" << solid_sim_->com_mod.tnNo << "n/" << solid_sim_->com_mod.eq[0].dof << "dof"
              << "  mesh=" << mesh_sim_->com_mod.tnNo << "n/" << mesh_sim_->com_mod.eq[0].dof << "dof"
              << std::endl;

    // Open coupling log file
    std::string log_dir = fluid_sim_->get_chnl_mod().appPath;
    coupling_log_.open(log_dir + "coupling.dat");
    char hdr[256];
    snprintf(hdr, sizeof(hdr), "# %4s %3s %10s %5s %10s %10s %10s %10s",
             "cTS", "cp", "time", "dB", "Ri/R1", "Ri/R0", "omega", "|disp|");
    coupling_log_ << hdr << std::endl;
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

  // ---- Sanity checks ----
  verify_node_maps();
}

//----------------------------------------------------------------------
// verify_node_maps — sanity checks for interface coupling
//----------------------------------------------------------------------
void PartitionedFSI::verify_node_maps()
{
  auto& cm = main_sim_->com_mod.cm;
  auto& cm_mod = main_sim_->cm_mod;
  if (!cm.mas(cm_mod)) return;

  const int nsd = main_sim_->com_mod.nsd;
  auto& solid_com = solid_sim_->com_mod;
  auto& fluid_com = fluid_sim_->com_mod;
  auto& mesh_com  = mesh_sim_->com_mod;

  std::cout << "[PartitionedFSI] Running interface sanity checks..." << std::endl;

  // Check 1: Coordinate round-trip (solid → fluid → solid)
  // Extract solid face coordinates, transfer to fluid face, compare with
  // actual fluid face coordinates
  {
    Array<double> solid_coords(nsd, solid_face_->nNo);
    for (int a = 0; a < solid_face_->nNo; a++) {
      int Ac = solid_face_->gN(a);
      for (int i = 0; i < nsd; i++)
        solid_coords(i, a) = solid_com.x(i, Ac);
    }

    auto fluid_coords_transferred = transfer_data(solid_to_fluid_map_,
                                                   solid_coords, fluid_face_->nNo);

    double max_err = 0.0;
    int n_checked = 0;
    for (int b = 0; b < fluid_face_->nNo; b++) {
      int Bc = fluid_face_->gN(b);
      // Check if this node was mapped to
      bool mapped = false;
      for (int a = 0; a < solid_face_->nNo; a++) {
        if (solid_to_fluid_map_[a] == b) { mapped = true; break; }
      }
      if (!mapped) continue;
      n_checked++;
      for (int i = 0; i < nsd; i++) {
        double err = std::abs(fluid_coords_transferred(i, b) - fluid_com.x(i, Bc));
        max_err = std::max(max_err, err);
      }
    }
    std::cout << "  Check 1 (coord transfer solid→fluid): max_err=" << max_err
              << " (" << n_checked << " nodes)" << std::endl;
  }

  // Check 2: solid→mesh uses same map as solid→fluid (mesh face = fluid face)
  std::cout << "  Check 2 (solid→mesh = solid→fluid, same sub-sim)" << std::endl;

  // Check 3: Round-trip (solid → fluid → solid)
  {
    Array<double> ones(nsd, solid_face_->nNo);
    for (int a = 0; a < solid_face_->nNo; a++)
      for (int i = 0; i < nsd; i++)
        ones(i, a) = 1.0 + 0.001 * a;  // Unique per node to detect scrambling

    auto on_fluid = transfer_data(solid_to_fluid_map_, ones, fluid_face_->nNo);
    auto back_on_solid = transfer_data(fluid_to_solid_map_, on_fluid, solid_face_->nNo);

    double max_err = 0.0;
    for (int a = 0; a < solid_face_->nNo; a++)
      for (int i = 0; i < nsd; i++)
        max_err = std::max(max_err, std::abs(back_on_solid(i, a) - ones(i, a)));
    std::cout << "  Check 3 (round-trip solid→fluid→solid): max_err=" << max_err << std::endl;
  }

  // Check 4: Print a few sample node coordinates to eyeball
  {
    std::cout << "  Check 4 (sample nodes, first 3):" << std::endl;
    for (int a = 0; a < std::min(3, solid_face_->nNo); a++) {
      int Ac_s = solid_face_->gN(a);
      int b = solid_to_fluid_map_[a];
      int Ac_f = (b >= 0) ? fluid_face_->gN(b) : -1;
      std::cout << "    solid[" << a << "] gN=" << Ac_s << " ("
                << solid_com.x(0, Ac_s) << ", "
                << solid_com.x(1, Ac_s) << ", "
                << solid_com.x(2, Ac_s) << ") → fluid[" << b << "] gN=" << Ac_f;
      if (Ac_f >= 0) {
        std::cout << " (" << fluid_com.x(0, Ac_f) << ", "
                  << fluid_com.x(1, Ac_f) << ", "
                  << fluid_com.x(2, Ac_f) << ")";
      }
      std::cout << std::endl;
    }
  }

  std::cout << "[PartitionedFSI] Sanity checks done." << std::endl;
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

    if (!converged && cm.mas(cm_mod)) {
      std::cout << "  TIME STEP " << cTS << " FAILED (NaN or no convergence)" << std::endl;
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

  // Save predictor state
  struct SavedState { Array<double> An, Yn, Dn; };
  auto save = [](SolutionStates& s) -> SavedState {
    return {s.current.get_acceleration(), s.current.get_velocity(), s.current.get_displacement()};
  };
  auto restore = [](SolutionStates& s, const SavedState& st) {
    s.current.get_acceleration() = st.An;
    s.current.get_velocity() = st.Yn;
    s.current.get_displacement() = st.Dn;
  };
  SavedState fluid_pred = save(fluid_sol);
  SavedState solid_pred = save(solid_sol);
  SavedState mesh_pred  = save(mesh_sol);

  // Initial displacement/velocity from predictor
  auto disp_current = fsi_coupling::extract_solid_displacement(
      solid_com, solid_com.eq[0], *solid_face_, solid_sol);
  disp_prev_ = disp_current;

  bool converged = false;

  for (int cp = 0; cp < config_.max_coupling_iterations; cp++) {

    // Restore all sub-sims to predictor state
    restore(fluid_sol, fluid_pred);
    restore(solid_sol, solid_pred);
    restore(mesh_sol, mesh_pred);

    if (cm.mas(cm_mod) && cp == 0) {
      std::cout << std::string(69, '-') << std::endl;
      std::cout << " Eq     N-i     T       dB  Ri/R1   Ri/R0    R/Ri     lsIt   dB  %t" << std::endl;
      std::cout << std::string(69, '-') << std::endl;
    }

    // ---- 1. FLUID SOLVE with solid velocity at wall ----
    // Prescribe solid velocity at lumen_wall, overwriting the Dir BC (zero).
    auto solid_vel = fsi_coupling::extract_solid_velocity(
        solid_com, solid_com.eq[0], *solid_face_, solid_sol);
    auto fluid_vel = transfer_data(solid_to_fluid_map_, solid_vel, fluid_face_->nNo);

    set_bc::set_bc_dir(fluid_com, fluid_sol);
    fsi_coupling::apply_velocity_on_fluid(
        fluid_com, fluid_com.eq[0], *fluid_face_, fluid_vel, fluid_sol);

    fluid_int.step_equation(0, [&]() {
      fsi_coupling::enforce_dirichlet_dofs_on_face(fluid_com, *fluid_face_, 0, nsd);
    });
    if (has_nan(fluid_sol)) {
      if (cm.mas(cm_mod)) std::cout << "  ABORT: NaN in fluid solve" << std::endl;
      return false;
    }

    // ---- 2. Extract traction ----
    auto fluid_traction = fsi_coupling::extract_fluid_traction(
        fluid_com, fluid_sim_->cm_mod,
        *fluid_mesh_, *fluid_face_, fluid_com.eq[0],
        fluid_int.get_Yg(), fluid_int.get_Dg(), fluid_sol);
    auto solid_traction = transfer_data(fluid_to_solid_map_,
                                        fluid_traction, solid_face_->nNo);

    // ---- 3. SOLID SOLVE with traction ----
    set_bc::set_bc_dir(solid_com, solid_sol);
    solid_int.step_equation(0, [&]() {
      fsi_coupling::apply_traction_on_solid(
          solid_com, solid_com.eq[0], *solid_face_, solid_traction);
    });
    if (has_nan(solid_sol)) {
      if (cm.mas(cm_mod)) std::cout << "  ABORT: NaN in solid solve" << std::endl;
      return false;
    }

    // ---- 4. Extract displacement, relax ----
    disp_current = fsi_coupling::extract_solid_displacement(
        solid_com, solid_com.eq[0], *solid_face_, solid_sol);

    for (int a = 0; a < solid_face_->nNo; a++)
      for (int i = 0; i < nsd; i++)
        disp_prev_(i, a) += omega_ * (disp_current(i, a) - disp_prev_(i, a));

    // ---- 5. MESH SOLVE with relaxed displacement ----
    auto mesh_disp = transfer_data(solid_to_mesh_map_, disp_prev_, mesh_face_->nNo);
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

    // ---- 6. Deform fluid mesh for next iteration ----
    auto& mesh_Dg = mesh_int.get_Dg();
    for (int a = 0; a < fluid_com.tnNo; a++)
      for (int i = 0; i < nsd; i++)
        fluid_com.x(i, a) += mesh_Dg(i, a);

    // Store deformed state in fluid predictor for next iteration restore
    fluid_pred = save(fluid_sol);

    // ---- Convergence check ----
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

    // Check for NaN
    if (std::isnan(rel) || std::isnan(disp_norm)) {
      if (cm.mas(cm_mod)) {
        std::cout << " CP " << cTS << "-" << cp + 1 << "  ABORT: NaN detected" << std::endl;
      }
      return false;
    }

    // Compute convergence metrics matching solver output style
    // dB = 20*log10(Ri/R0), Ri/R1 = rel_change / first_rel_change
    int dB_val = 0;
    double ri_r1 = 1.0;  // ratio to first iteration residual

    if (cp == 0) {
      first_res_norm_ = res_norm;
    }

    if (first_res_norm_ > 1e-30 && res_norm > 0) {
      ri_r1 = res_norm / first_res_norm_;
      dB_val = static_cast<int>(20.0 * log10(ri_r1));
    }

    // Format: CP  cTS-cp  time  [dB  Ri/R1  Ri/R0  omega]
    // Matches: EQ  cTS-N  time  [dB  Ri/R1  Ri/R0  R/Ri]
    if (cm.mas(cm_mod)) {
      bool conv = rel < config_.coupling_tolerance;
      char buf[256];
      snprintf(buf, sizeof(buf),
               " CP %d-%d%s %4.3e  [%d %4.3e %4.3e %4.3e]",
               cTS, cp + 1, conv ? "s" : " ",
               main_sim_->com_mod.timer.get_elapsed_time(),
               dB_val, ri_r1, rel, omega_);
      std::cout << buf << std::endl;

      // Write to coupling log file
      if (coupling_log_.is_open()) {
        char log_buf[256];
        snprintf(log_buf, sizeof(log_buf), "  %4d %3d %10.3e %5d %10.3e %10.3e %10.3e %10.3e",
                 cTS, cp + 1, main_sim_->com_mod.timer.get_elapsed_time(),
                 dB_val, ri_r1, rel, omega_, disp_norm);
        coupling_log_ << log_buf << std::endl;
      }
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
