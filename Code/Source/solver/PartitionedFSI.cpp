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
    const char* method_name = "unknown";
    switch (config_.coupling_method) {
      case CouplingMethod::constant: method_name = "constant"; break;
      case CouplingMethod::aitken:   method_name = "aitken"; break;
      case CouplingMethod::iqn_ils:  method_name = "iqn-ils"; break;
    }
    std::cout << "[PartitionedFSI] Sub-sims ready:"
              << "  fluid=" << fluid_sim_->com_mod.tnNo << "n/" << fluid_sim_->com_mod.tDof << "tDof"
              << "  solid=" << solid_sim_->com_mod.tnNo << "n/" << solid_sim_->com_mod.eq[0].dof << "dof"
              << "  mesh=" << mesh_sim_->com_mod.tnNo << "n/" << mesh_sim_->com_mod.eq[0].dof << "dof"
              << "  coupling=" << method_name;
    if (config_.coupling_method == CouplingMethod::iqn_ils)
      std::cout << " (q=" << config_.iqn_ils_q
                << ", eps=" << config_.iqn_ils_eps
                << ", warmup=" << config_.iqn_ils_warmup << ")";
    std::cout << std::endl;

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
// relax_interface — updates disp_prev_ and vel_prev_
//----------------------------------------------------------------------
void PartitionedFSI::relax_interface(int cp, int nsd,
                                     const Array<double>& disp_current,
                                     const Array<double>& vel_current)
{
  switch (config_.coupling_method) {
    case CouplingMethod::constant:
      relax_constant(cp, nsd, disp_current, vel_current);
      break;
    case CouplingMethod::aitken:
      relax_aitken(cp, nsd, disp_current, vel_current);
      break;
    case CouplingMethod::iqn_ils:
      relax_iqn_ils(cp, nsd, disp_current, vel_current);
      break;
  }
}

//----------------------------------------------------------------------
// relax_constant — fixed relaxation
//----------------------------------------------------------------------
void PartitionedFSI::relax_constant(int cp, int nsd,
                                     const Array<double>& disp_current,
                                     const Array<double>& vel_current)
{
  omega_ = config_.initial_relaxation;
  for (int a = 0; a < solid_face_->nNo; a++)
    for (int i = 0; i < nsd; i++) {
      disp_prev_(i, a) += omega_ * (disp_current(i, a) - disp_prev_(i, a));
      vel_prev_(i, a)  += omega_ * (vel_current(i, a) - vel_prev_(i, a));
    }
}

//----------------------------------------------------------------------
// relax_aitken — Aitken Delta^2 (Küttler & Wall 2008, Eq. 44)
//----------------------------------------------------------------------
void PartitionedFSI::relax_aitken(int cp, int nsd,
                                   const Array<double>& disp_current,
                                   const Array<double>& vel_current)
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
    for (int i = 0; i < nsd; i++) {
      disp_prev_(i, a) += omega_ * (disp_current(i, a) - disp_prev_(i, a));
      vel_prev_(i, a)  += omega_ * (vel_current(i, a) - vel_prev_(i, a));
    }
}

//----------------------------------------------------------------------
// relax_iqn_ils — IQN-ILS following StanfordCBCL/svFSGe implementation
//
// Columns persist across time steps (trimmed to iqn_ils_q max).
// First time step uses Aitken. QR filtering with eps threshold
// removes linearly dependent columns.
//----------------------------------------------------------------------
void PartitionedFSI::relax_iqn_ils(int cp, int nsd,
                                    const Array<double>& disp_current,
                                    const Array<double>& vel_current)
{
  const int n = nsd * solid_face_->nNo;
  const int cTS = main_sim_->com_mod.cTS;

  // Current x_tilde (unrelaxed solver output) and residual
  std::vector<double> x_tilde(n), r(n);
  for (int a = 0; a < solid_face_->nNo; a++)
    for (int i = 0; i < nsd; i++) {
      int idx = a * nsd + i;
      x_tilde[idx] = disp_current(i, a);
      r[idx] = disp_current(i, a) - disp_prev_(i, a);
    }

  // Store history for difference vectors
  x_tilde_hist_.push_back(x_tilde);
  r_hist_.push_back(r);

  // Append difference vectors (need at least 2 history entries)
  if (x_tilde_hist_.size() >= 2) {
    auto& xt1 = x_tilde_hist_[x_tilde_hist_.size() - 1];
    auto& xt0 = x_tilde_hist_[x_tilde_hist_.size() - 2];
    auto& r1 = r_hist_[r_hist_.size() - 1];
    auto& r0 = r_hist_[r_hist_.size() - 2];
    std::vector<double> dw(n), dv(n);
    for (int j = 0; j < n; j++) {
      dw[j] = xt1[j] - xt0[j];
      dv[j] = r1[j] - r0[j];
    }
    W_cols_.push_back(dw);
    V_cols_.push_back(dv);
  }

  // Use Aitken until we have enough V/W columns
  if (static_cast<int>(V_cols_.size()) < config_.iqn_ils_warmup) {
    relax_aitken(cp, nsd, disp_current, vel_current);
    return;
  }

  // Trim to max columns
  const int nq = config_.iqn_ils_q;
  while (static_cast<int>(V_cols_.size()) > nq) {
    V_cols_.erase(V_cols_.begin());
    W_cols_.erase(W_cols_.begin());
  }

  int q = static_cast<int>(V_cols_.size());

  // QR decomposition with filtering (modified Gram-Schmidt)
  // Removes columns where ||v_orth|| < eps * ||v_orig||
  const double eps = config_.iqn_ils_eps;
  // Work on copies so we can remove columns
  auto V_work = V_cols_;
  auto W_work = W_cols_;

  // Iterative QR with column removal (matching svFSGe QRfiltering_mod)
  std::vector<std::vector<double>> Q;
  std::vector<std::vector<double>> R_mat;
  bool restart;

  do {
    restart = false;
    q = static_cast<int>(V_work.size());
    Q.assign(q, std::vector<double>(n, 0.0));
    R_mat.assign(q, std::vector<double>(q, 0.0));

    // First column
    double norm0 = 0;
    for (int k = 0; k < n; k++) norm0 += V_work[0][k] * V_work[0][k];
    norm0 = sqrt(norm0);
    if (norm0 < 1e-30) { V_work.erase(V_work.begin()); W_work.erase(W_work.begin()); restart = true; continue; }
    R_mat[0][0] = norm0;
    for (int k = 0; k < n; k++) Q[0][k] = V_work[0][k] / norm0;

    for (int j = 1; j < q; j++) {
      auto vbar = V_work[j];
      double orig_norm = 0;
      for (int k = 0; k < n; k++) orig_norm += vbar[k] * vbar[k];
      orig_norm = sqrt(orig_norm);

      for (int i = 0; i < j; i++) {
        double dot = 0;
        for (int k = 0; k < n; k++) dot += Q[i][k] * vbar[k];
        R_mat[i][j] = dot;
        for (int k = 0; k < n; k++) vbar[k] -= dot * Q[i][k];
      }

      double norm = 0;
      for (int k = 0; k < n; k++) norm += vbar[k] * vbar[k];
      norm = sqrt(norm);

      if (norm < eps * orig_norm) {
        // Linearly dependent: remove and restart QR
        V_work.erase(V_work.begin() + j);
        W_work.erase(W_work.begin() + j);
        restart = true;
        break;
      }
      R_mat[j][j] = norm;
      for (int k = 0; k < n; k++) Q[j][k] = vbar[k] / norm;
    }
  } while (restart && !V_work.empty());

  // Update stored matrices after filtering
  V_cols_ = V_work;
  W_cols_ = W_work;
  q = static_cast<int>(V_cols_.size());

  if (q == 0) {
    relax_aitken(cp, nsd, disp_current, vel_current);
    return;
  }

  // Solve R * c = Q^T * (-r)
  std::vector<double> rhs(q);
  for (int j = 0; j < q; j++) {
    double dot = 0;
    for (int k = 0; k < n; k++) dot += Q[j][k] * (-r[k]);
    rhs[j] = dot;
  }

  std::vector<double> c(q, 0.0);
  for (int j = q - 1; j >= 0; j--) {
    c[j] = rhs[j];
    for (int k = j + 1; k < q; k++) c[j] -= R_mat[j][k] * c[k];
    c[j] /= R_mat[j][j];
  }

  // Check for NaN/Inf
  for (int j = 0; j < q; j++) {
    if (std::isnan(c[j]) || std::isinf(c[j])) {
      V_cols_.clear(); W_cols_.clear();
      relax_aitken(cp, nsd, disp_current, vel_current);
      return;
    }
  }

  // Update: x_{k+1} = x_tilde + W * c  (svFSGe formula)
  for (int a = 0; a < solid_face_->nNo; a++)
    for (int i = 0; i < nsd; i++) {
      int idx = a * nsd + i;
      double correction = 0;
      for (int j = 0; j < q; j++)
        correction += W_cols_[j][idx] * c[j];
      disp_prev_(i, a) = disp_current(i, a) + correction;
    }

  // Velocity: full step (follows displacement via Newmark)
  vel_prev_ = vel_current;

  // omega_ for logging
  double corr2 = 0, res2 = 0;
  for (int j = 0; j < n; j++) res2 += r[j] * r[j];
  for (int a = 0; a < solid_face_->nNo; a++)
    for (int i = 0; i < nsd; i++) {
      double d = disp_prev_(i, a) - disp_current(i, a);
      corr2 += d * d;
    }
  omega_ = (res2 > 1e-30) ? sqrt(corr2 / res2) : 1.0;
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
  r_prev_.clear();

  // Clear per-time-step history (V/W persist across time steps for IQN-ILS)
  x_tilde_hist_.clear();
  r_hist_.clear();

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

  // Save reference mesh coordinates (restored each coupling iteration)
  Array<double> x_ref(fluid_com.x);

  // ALE mesh velocity: save predictor mesh velocity for injection into fluid.
  // mesh_vel_Yn is updated after each mesh solve so the next coupling
  // iteration's fluid sees the latest mesh motion.
  const int mesh_s = mesh_com.eq[0].s;  // mesh DOF offset (should be 0)
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

  // Initial displacement and velocity from predictor
  auto disp_current = fsi_coupling::extract_solid_displacement(
      solid_com, solid_com.eq[0], *solid_face_, solid_sol);
  auto vel_current = fsi_coupling::extract_solid_velocity(
      solid_com, solid_com.eq[0], *solid_face_, solid_sol);
  disp_prev_ = disp_current;
  vel_prev_ = vel_current;

  bool converged = false;

  for (int cp = 0; cp < config_.max_coupling_iterations; cp++) {

    // Restore all sub-sims to predictor state AND mesh coordinates
    restore(fluid_sol, fluid_pred);
    restore(solid_sol, solid_pred);
    restore(mesh_sol, mesh_pred);
    fluid_com.x = x_ref;

    if (cm.mas(cm_mod) && cp == 0) {
      if (histor_log_.is_open()) {
        histor_log_ << std::string(69, '-') << std::endl;
        histor_log_ << " Eq     N-i     T       dB  Ri/R1   Ri/R0    R/Ri     lsIt   dB  %t" << std::endl;
        histor_log_ << std::string(69, '-') << std::endl;
      }
    }

    // ---- 1. FLUID SOLVE with relaxed wall velocity ----
    // vel_prev_ comes directly from the solid solver (no Newmark recomputation).
    auto fluid_vel = transfer_data(solid_to_fluid_map_, vel_prev_, fluid_face_->nNo);

    set_bc::set_bc_dir(fluid_com, fluid_sol);
    fsi_coupling::apply_velocity_on_fluid(
        fluid_com, fluid_com.eq[0], *fluid_face_, fluid_vel, fluid_sol);

    // Set ALE mesh velocity for the fluid assembly.
    // construct_fluid and b_assem_neu_bc check this array and extend the local
    // element velocity array yl with mesh velocity at DOFs nsd+1..2*nsd.
    {
      double af = fluid_com.eq[0].af;
      fluid_com.ale_mesh_velocity.resize(nsd, fluid_com.tnNo);
      for (int a = 0; a < fluid_com.tnNo; a++)
        for (int i = 0; i < nsd; i++)
          fluid_com.ale_mesh_velocity(i, a) = (1.0 - af) * mesh_vel_Yo(i, a)
                                             + af * mesh_vel_Yn(i, a);
    }

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

    // ---- 4. Extract displacement AND velocity from solid ----
    disp_current = fsi_coupling::extract_solid_displacement(
        solid_com, solid_com.eq[0], *solid_face_, solid_sol);
    auto vel_current = fsi_coupling::extract_solid_velocity(
        solid_com, solid_com.eq[0], *solid_face_, solid_sol);

    // ---- Convergence check (BEFORE relaxation) ----
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

    // ---- 5. Relaxation (updates disp_prev_ and vel_prev_) ----
    relax_interface(cp, nsd, disp_current, vel_current);

    // Check for NaN/large values in relaxed displacement
    {
      bool has_nan_relax = false;
      double max_disp = 0;
      for (int a = 0; a < solid_face_->nNo && !has_nan_relax; a++)
        for (int i = 0; i < nsd; i++) {
          if (std::isnan(disp_prev_(i, a)) || std::isinf(disp_prev_(i, a)))
            { has_nan_relax = true; break; }
          max_disp = std::max(max_disp, std::abs(disp_prev_(i, a)));
        }
      if (has_nan_relax || max_disp > 1e10) {
        if (cm.mas(cm_mod)) std::cout << "  ABORT: NaN/divergence after relaxation (max_disp=" << max_disp << ")" << std::endl;
        return false;
      }
    }

    // ---- 6. MESH SOLVE with relaxed displacement ----
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

    // Update ALE mesh velocity for next coupling iteration's fluid solve
    {
      auto& mYn = mesh_sol.current.get_velocity();
      for (int a = 0; a < mesh_com.tnNo; a++)
        for (int i = 0; i < nsd; i++)
          mesh_vel_Yn(i, a) = mYn(mesh_s + i, a);
    }

    // ---- 7. Deform fluid mesh for next iteration ----
    // The mesh equation solves for TOTAL displacement from the original
    // reference.  x_ref already contains the old displacement (x_original + Do),
    // so we must apply only the INCREMENT (Dn - Do) to avoid accumulating
    // total displacements across time steps.
    {
      auto& mesh_Dn = mesh_sol.current.get_displacement();
      auto& mesh_Do = mesh_sol.old.get_displacement();
      for (int a = 0; a < fluid_com.tnNo; a++)
        for (int i = 0; i < nsd; i++)
          fluid_com.x(i, a) = x_ref(i, a)
                             + mesh_Dn(i + mesh_s, a) - mesh_Do(i + mesh_s, a);
    }

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

    if (cm.mas(cm_mod)) {
      bool conv = rel < config_.coupling_tolerance;
      double time_elapsed = main_sim_->com_mod.timer.get_elapsed_time();

      // Screen: tabular format (same as coupling.dat) for easy logging
      char buf[256];
      snprintf(buf, sizeof(buf), "  %4d %3d %10.3e %5d %10.3e %10.3e %10.3e %10.3e",
               cTS, cp + 1, time_elapsed,
               dB_val, ri_r1, rel, omega_, disp_norm);
      std::cout << buf << std::endl;

      // coupling.dat: same tabular format
      if (coupling_log_.is_open()) {
        coupling_log_ << buf << std::endl;
      }

      // histor.dat: compact solver-style format
      if (histor_log_.is_open()) {
        char hbuf[256];
        snprintf(hbuf, sizeof(hbuf),
                 " CP %d-%d%s %4.3e  [%d %4.3e %4.3e %4.3e]",
                 cTS, cp + 1, conv ? "s" : " ",
                 time_elapsed,
                 dB_val, ri_r1, rel, omega_);
        histor_log_ << hbuf << std::endl;
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
