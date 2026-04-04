// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PARTITIONED_FSI_H
#define PARTITIONED_FSI_H

#include "Simulation.h"
#include "Integrator.h"
#include "Array.h"

#include <fstream>
#include <memory>
#include <vector>
#include <string>

/// @brief Configuration for partitioned FSI coupling, read from XML input.
struct PartitionedFSIConfig {
  int max_coupling_iterations = 50;
  double coupling_tolerance = 1e-6;
  double initial_relaxation = 1.0;
  bool use_aitken = true;

  // Face names for the FSI interface
  std::string fluid_interface_face;
  std::string solid_interface_face;

  // Paths to standalone XML files for each sub-field
  std::string fluid_xml;
  std::string solid_xml;
  std::string mesh_xml;
};

/// @brief Partitioned FSI coupling with 3 independent sub-Simulations.
///
/// Each sub-field (fluid, struct, mesh) has its own Simulation object with
/// independent mesh, solution arrays, and linear system. No shared global
/// arrays, no DOF offsets (each eq.s=0), no regularization of inactive nodes.
///
/// Implements Dirichlet-Neumann coupling with Aitken relaxation:
///   1. Transfer solid displacement to mesh interface, solve mesh equation
///   2. Deform fluid mesh using mesh displacement, solve fluid equation
///   3. Extract fluid traction, apply to solid, solve solid equation
///   4. Extract solid displacement, apply Aitken relaxation
///   5. Check coupling convergence
///
/// Related to GitHub issue #431: Implement partitioned FSI in svMultiPhysics
class PartitionedFSI {
public:
  /// @brief Construct partitioned FSI with 3 sub-simulations.
  /// Creates temp XML files, initializes each sub-sim through the standard
  /// pipeline (read_files → distribute → initialize → add_eq_linear_algebra).
  PartitionedFSI(Simulation* main_simulation, const PartitionedFSIConfig& config,
                 const std::string& xml_file_path);

  ~PartitionedFSI();

  /// @brief Run the complete partitioned FSI time-stepping loop.
  /// Replaces iterate_solution() for partitioned coupling.
  void run();

  /// @brief Execute one time step of partitioned FSI coupling.
  /// Called from within run() after predictor and set_bc_dir.
  /// @return true if coupling converged within max iterations
  bool step();

private:
  Simulation* main_sim_;
  PartitionedFSIConfig config_;
  std::string xml_file_path_;

  // Sub-simulations (owned)
  std::unique_ptr<Simulation> fluid_sim_;
  std::unique_ptr<Simulation> solid_sim_;
  std::unique_ptr<Simulation> mesh_sim_;

  // Interface face pointers within each sub-sim
  const faceType* fluid_face_ = nullptr;   // fluid_interface_face in fluid_sim
  const faceType* solid_face_ = nullptr;   // solid_interface_face in solid_sim
  const faceType* mesh_face_ = nullptr;    // fluid_interface_face in mesh_sim

  // Mesh pointers within each sub-sim
  const mshType* fluid_mesh_ = nullptr;
  const mshType* solid_mesh_ = nullptr;
  const mshType* mesh_mesh_ = nullptr;

  // Node maps between interface faces: src_local_idx → tgt_local_idx
  std::vector<int> solid_to_fluid_map_;   // wall_inner → lumen_wall (fluid)
  std::vector<int> fluid_to_solid_map_;   // lumen_wall (fluid) → wall_inner
  std::vector<int> solid_to_mesh_map_;    // wall_inner → lumen_wall (mesh)

  // Coupling state
  Array<double> disp_prev_;
  Array<double> vel_prev_;
  double omega_;
  double first_res_norm_ = 0.0;

  // IQN-ILS state (Degroote 2013, Algorithm 10)
  std::vector<std::vector<double>> V_cols_;  // residual differences
  std::vector<std::vector<double>> W_cols_;  // displacement differences
  std::vector<double> r_prev_;               // previous residual
  std::vector<double> x_tilde_prev_;         // previous S(F(x))

  // Output file for coupling convergence history
  std::ofstream coupling_log_;

  // Temp XML file paths (cleaned up in destructor)
  std::vector<std::string> temp_xml_paths_;

  /// Create and initialize the 3 sub-simulations from temp XML files
  void create_sub_simulations();

  /// Generate a temporary XML file for one sub-field ("fluid", "struct", "mesh")
  std::string generate_sub_xml(const std::string& field_type);

  /// Initialize one sub-simulation through the standard pipeline
  void init_sub_simulation(Simulation* sim, const std::string& xml_path);

  /// Resolve face/mesh pointers within the initialized sub-simulations
  void resolve_faces();

  /// Build coordinate-based node maps between interface faces of different sub-sims
  void build_node_maps();

  /// Run sanity checks on node maps and data transfer
  void verify_node_maps();

  /// Build a one-directional node map from face_a to face_b using coordinate matching
  static void build_face_node_map(const faceType& face_a, const ComMod& com_a,
                                  const faceType& face_b, const ComMod& com_b,
                                  std::vector<int>& a_to_b);

  /// Transfer data from source face to target face using pre-built node map
  static Array<double> transfer_data(const std::vector<int>& src_to_tgt_map,
                                     const Array<double>& src_data, int tgt_nNo);

  /// Compute Aitken relaxation factor
  static double compute_aitken_omega(const Array<double>& residual,
                                     const Array<double>& residual_prev,
                                     double omega_prev);

  /// Save VTK and restart output for all sub-simulations
  void save_results();
};

#endif // PARTITIONED_FSI_H
