// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PARTITIONED_FSI_H
#define PARTITIONED_FSI_H

#include "Simulation.h"
#include "Integrator.h"
#include "Array.h"

/// @brief Configuration for partitioned FSI coupling, read from XML input.
struct PartitionedFSIConfig {
  int max_coupling_iterations = 50;
  double coupling_tolerance = 1e-6;
  double initial_relaxation = 1.0;
  bool use_aitken = true;

  // Equation indices (auto-detected from equation types)
  int fluid_eq_index = -1;
  int solid_eq_index = -1;
  int mesh_eq_index = -1;

  // Face names for the FSI interface
  std::string fluid_interface_face;
  std::string solid_interface_face;
};

/// @brief Partitioned FSI coupling orchestrator.
///
/// Implements Dirichlet-Neumann partitioned coupling with Aitken relaxation:
///   1. Apply interface displacement to fluid/mesh (Dirichlet)
///   2. Solve mesh equation
///   3. Solve fluid equation
///   4. Extract fluid traction at interface
///   5. Solve solid equation with traction (Neumann)
///   6. Extract solid displacement at interface
///   7. Apply Aitken relaxation
///   8. Check coupling convergence
///
/// Related to GitHub issue #431: Implement partitioned FSI in svMultiPhysics
class PartitionedFSI {
public:
  PartitionedFSI(Simulation* simulation, Integrator* integrator,
                 const PartitionedFSIConfig& config);

  /// @brief Execute one time step of partitioned FSI coupling.
  /// Call this after predictor() and set_bc_dir().
  /// @return true if coupling converged within max iterations
  bool step();

private:
  Simulation* simulation_;
  Integrator* integrator_;
  PartitionedFSIConfig config_;

  // Mesh and face references (resolved from config names)
  const mshType* fluid_mesh_ = nullptr;
  const faceType* fluid_face_ = nullptr;
  const faceType* solid_face_ = nullptr;

  // Interface data from previous coupling iteration
  Array<double> disp_prev_;
  Array<double> disp_residual_prev_;

  // Aitken relaxation parameter
  double omega_;

  /// Resolve mesh/face pointers from config names
  void resolve_faces();

  /// Compute Aitken relaxation factor
  double compute_aitken_omega(const Array<double>& residual,
                              const Array<double>& residual_prev,
                              double omega_prev);
};

#endif // PARTITIONED_FSI_H
