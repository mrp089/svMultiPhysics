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

/// @brief Coupling method for interface relaxation
enum class CouplingMethod { constant, aitken };

/// @brief Configuration for partitioned FSI coupling, read from XML input.
struct PartitionedFSIConfig {
  int max_coupling_iterations = 50;
  double coupling_tolerance = 1e-6;
  double initial_relaxation = 1.0;
  double omega_max = 1.0;
  CouplingMethod coupling_method = CouplingMethod::aitken;

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
  PartitionedFSI(Simulation* main_simulation, const PartitionedFSIConfig& config,
                 const std::string& xml_file_path);

  ~PartitionedFSI();

  void run();
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
  const faceType* fluid_face_ = nullptr;
  const faceType* solid_face_ = nullptr;
  const faceType* mesh_face_ = nullptr;

  // Mesh pointers within each sub-sim
  const mshType* fluid_mesh_ = nullptr;
  const mshType* solid_mesh_ = nullptr;
  const mshType* mesh_mesh_ = nullptr;

  // Node maps between interface faces: src_local_idx → tgt_local_idx
  std::vector<int> solid_to_fluid_map_;
  std::vector<int> fluid_to_solid_map_;
  std::vector<int> solid_to_mesh_map_;

  // Coupling state
  Array<double> disp_prev_;
  Array<double> vel_prev_;
  double omega_;
  double first_res_norm_ = 0.0;

  // Aitken state
  std::vector<double> r_prev_;

  // Output files for coupling convergence history
  std::ofstream coupling_log_;
  std::ofstream histor_log_;

  // Temp XML file paths (cleaned up in destructor)
  std::vector<std::string> temp_xml_paths_;

  void resolve_faces();
  void build_node_maps();

  /// Solve fluid equation with current interface velocity and ALE mesh velocity
  bool solve_fluid(const Array<double>& mesh_vel_Yo, const Array<double>& mesh_vel_Yn);

  /// Extract fluid traction, transfer to solid, solve solid equation
  bool solve_solid();

  /// Solve mesh equation with relaxed displacement, deform fluid mesh
  bool solve_mesh(const Array<double>& x_ref, int mesh_s);

  /// Compute vel_prev_ from disp_prev_ using Newmark relationship
  void compute_interface_velocity();

  void relax_interface(int cp, int nsd, const Array<double>& disp_current);
  void relax_constant(int cp, int nsd, const Array<double>& disp_current);
  void relax_aitken(int cp, int nsd, const Array<double>& disp_current);

  static void build_face_node_map(const faceType& face_a, const ComMod& com_a,
                                  const faceType& face_b, const ComMod& com_b,
                                  std::vector<int>& a_to_b);

  static Array<double> transfer_data(const std::vector<int>& src_to_tgt_map,
                                     const Array<double>& src_data, int tgt_nNo);

  void save_results();
};

#endif // PARTITIONED_FSI_H
