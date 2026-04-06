// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FSI_COUPLING_H
#define FSI_COUPLING_H

#include "ComMod.h"
#include "SolutionStates.h"

/// @brief FSI interface data exchange functions for partitioned coupling.
///
/// These functions extract and apply fluid traction and solid displacement
/// at the FSI interface, enabling partitioned (Dirichlet-Neumann) coupling
/// between separately solved fluid and solid equations.
///
/// Related to GitHub issue #431: Implement partitioned FSI in svMultiPhysics

namespace fsi_coupling {

/// @brief Extract consistent nodal traction forces from fluid at FSI interface.
///
/// Computes f(i,a) = integral(sigma_ij * n_j * N_a dGamma) at each face node,
/// where sigma is the Cauchy stress tensor of the fluid. The stress is computed
/// from the velocity gradient and pressure using volume-element shape functions,
/// then integrated over the face using face Gauss quadrature.
///
/// Sign convention: returns the force that the fluid exerts ON the solid,
/// i.e., -(sigma_fluid . n_fluid) where n_fluid points outward from the fluid.
///
/// @param com_mod Common module with global data
/// @param cm_mod Communication module
/// @param fluid_mesh The fluid volume mesh containing the interface face
/// @param fluid_face The fluid-side FSI interface face
/// @param fluid_eq The fluid equation (for domain and viscosity access)
/// @param Yg Solution variables at generalized-alpha level
/// @param Dg Integrated variables at generalized-alpha level
/// @param solutions Solution states at old and current time levels
/// @return Array(nsd, fluid_face.nNo) of consistent nodal forces
Array<double> extract_fluid_traction(
    ComMod& com_mod, const CmMod& cm_mod,
    const mshType& fluid_mesh, const faceType& fluid_face,
    const eqType& fluid_eq,
    const Array<double>& Yg, const Array<double>& Dg,
    const SolutionStates& solutions);

/// @brief Extract solid displacement at interface face nodes.
///
/// @param com_mod Common module
/// @param solid_eq The solid equation (for DOF offset)
/// @param solid_face The solid-side FSI interface face
/// @param solutions Solution states
/// @return Array(nsd, solid_face.nNo) of displacement values
Array<double> extract_solid_displacement(
    const ComMod& com_mod, const eqType& solid_eq,
    const faceType& solid_face, const SolutionStates& solutions);

/// @brief Apply velocity as strong Dirichlet BC on fluid interface nodes.
/// Directly sets Yn at the fluid equation DOF range for the face nodes.
void apply_velocity_on_fluid(
    ComMod& com_mod, const eqType& fluid_eq,
    const faceType& fluid_face,
    const Array<double>& velocity,
    SolutionStates& solutions);

/// @brief Apply pre-computed consistent nodal forces to the solid residual.
///
/// Adds the traction forces directly to com_mod.R at the global node locations
/// corresponding to the solid face. This should be called during the
/// post-assembly callback of step_equation() for the solid equation.
///
/// @param com_mod Common module (R is modified)
/// @param solid_eq The solid equation (for DOF offset)
/// @param solid_face The solid-side FSI interface face
/// @param traction Array(nsd, solid_face.nNo) of consistent nodal forces
void apply_traction_on_solid(
    ComMod& com_mod, const eqType& solid_eq,
    const faceType& solid_face,
    const Array<double>& traction);

/// @brief Apply displacement as strong Dirichlet BC on mesh interface nodes.
///
/// Directly sets the displacement in the solution arrays (An, Yn, Dn) for
/// the mesh equation DOF range at the interface face nodes.
///
/// @param com_mod Common module
/// @param mesh_eq The mesh equation (for DOF offset)
/// @param mesh_face The mesh-side FSI interface face
/// @param displacement Array(nsd, mesh_face.nNo) of displacement values
/// @param solutions Solution states (modified)
void apply_displacement_on_mesh(
    ComMod& com_mod, const eqType& mesh_eq,
    const faceType& mesh_face,
    const Array<double>& displacement,
    SolutionStates& solutions);

/// @brief Transfer nodal data between projected FSI interface faces.
///
/// Uses the shared global node IDs established by set_projector() to map
/// data from source face nodes to target face nodes. The faces must be
/// @brief Enforce Dirichlet BC at face nodes in the assembled linear system.
///
/// Zeros the residual and diagonalizes the system matrix rows for the
/// face nodes, so the linear solve produces zero correction there.
/// Call this in a post_assembly callback for step_equation().
void enforce_dirichlet_on_face(ComMod& com_mod, const faceType& lFa, int nsd);

/// @brief Enforce Dirichlet BC for specific DOFs at face nodes.
///
/// Like enforce_dirichlet_on_face but only modifies DOFs in the range
/// [dof_start, dof_start + num_dofs). Used for enforcing velocity Dirichlet
/// on a fluid interface without freezing the pressure DOF.
void enforce_dirichlet_dofs_on_face(ComMod& com_mod, const faceType& lFa,
                                    int dof_start, int num_dofs);

} // namespace fsi_coupling

#endif // FSI_COUPLING_H
