// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fsi_coupling.h"
#include "all_fun.h"
#include "fluid.h"
#include "fs.h"
#include "nn.h"
#include "utils.h"

#include <unordered_map>

namespace fsi_coupling {

//----------------------------------------------------------------------
// extract_fluid_traction
//----------------------------------------------------------------------
// Adapts the traction computation from post::bpost() (post.cpp:90-339)
// but computes consistent nodal forces via face Gauss integration
// instead of area-weighted nodal traction.
//
Array<double> extract_fluid_traction(
    ComMod& com_mod, const CmMod& cm_mod,
    const mshType& lM, const faceType& lFa,
    const eqType& eq,
    const Array<double>& Yg, const Array<double>& Dg,
    const SolutionStates& solutions)
{
  const int nsd = com_mod.nsd;
  const int eNoN = lM.eNoN;

  // Set up pressure function space (following bpost pattern, post.cpp:149-175)
  // For single function space (P1-P1), pressure uses same shape functions.
  // For Taylor-Hood (P2-P1), pressure uses the lower-order function space.
  fsType fsP;
  if (lM.nFs == 1) {
    fsP.eNoN = lM.fs[0].eNoN;
    fsP.N = lM.fs[0].N;
  } else {
    fsP.eNoN = lM.fs[1].eNoN;
    // For Taylor-Hood, evaluate pressure shape functions at velocity Gauss points
    fsP.nG = lM.fs[0].nG;
    fsP.eType = lM.fs[1].eType;
    fs::alloc_fs(fsP, nsd, nsd);
    fsP.xi = lM.fs[0].xi;
    for (int g = 0; g < fsP.nG; g++) {
      nn::get_gnn(nsd, fsP.eType, fsP.eNoN, g, fsP.xi, fsP.N, fsP.Nx);
    }
  }

  // Result: consistent nodal forces at face nodes
  Array<double> result(nsd, lFa.nNo);

  // Build reverse map: global node ID -> face-local index
  std::unordered_map<int, int> global_to_face;
  global_to_face.reserve(lFa.nNo);
  for (int a = 0; a < lFa.nNo; a++) {
    global_to_face[lFa.gN(a)] = a;
  }

  // Temporary arrays for volume element computation
  Array<double> xl(nsd, eNoN);    // element node coordinates
  Array<double> ul(nsd, eNoN);    // element node velocities
  Vector<double> pl(fsP.eNoN);    // element node pressures
  Array<double> Nx(nsd, eNoN);    // shape function gradients (spatial)
  Array<double> ks(nsd, nsd);     // metric tensor

  // Loop over face elements
  for (int e = 0; e < lFa.nEl; e++) {
    int Ec = lFa.gE(e);  // parent volume element

    // Get fluid domain for this element
    int cEq = com_mod.cEq;
    int cDmn = all_fun::domain(com_mod, lM, cEq, Ec);
    if (cDmn == -1) continue;

    // Load volume element nodal data (following bpost lines 199-218)
    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a, Ec);
      for (int i = 0; i < nsd; i++) {
        xl(i, a) = com_mod.x(i, Ac);
        ul(i, a) = Yg(i, Ac);  // velocity DOFs 0..nsd-1
      }
    }

    // Load pressure at element nodes
    for (int a = 0; a < fsP.eNoN; a++) {
      int Ac = lM.IEN(a, Ec);
      pl(a) = Yg(nsd, Ac);  // pressure is DOF nsd
    }

    // Compute element-averaged stress tensor from volume Gauss points
    // (exact for linear elements, approximate for higher order)
    Array<double> sigma_avg(nsd, nsd);
    double Jac_vol;

    for (int g = 0; g < lM.nG; g++) {
      if (g == 0 || !lM.lShpF) {
        auto lM_Nx = lM.Nx.slice(g);
        nn::gnn(eNoN, nsd, nsd, lM_Nx, xl, Nx, Jac_vol, ks);
      }

      auto N = lM.N.col(g);

      // Velocity gradient: ux(i,j) = du_j/dx_i
      Array<double> ux(nsd, nsd);
      for (int a = 0; a < eNoN; a++) {
        for (int i = 0; i < nsd; i++) {
          for (int j = 0; j < nsd; j++) {
            ux(i, j) += Nx(i, a) * ul(j, a);
          }
        }
      }

      // Pressure at Gauss point
      double p = 0.0;
      for (int a = 0; a < fsP.eNoN; a++) {
        p += fsP.N(a, g) * pl(a);
      }

      // Shear rate: gam = sqrt(2 * e_ij * e_ij)
      double gam = 0.0;
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          gam += (ux(i, j) + ux(j, i)) * (ux(i, j) + ux(j, i));
        }
      }
      gam = sqrt(0.5 * gam);

      // Get viscosity (handles non-Newtonian models)
      double mu, mu_s;
      fluid::get_viscosity(com_mod, eq.dmn[cDmn], gam, mu, mu_s, mu_s);

      // Accumulate stress: sigma = -p*I + mu*(grad_u + grad_u^T)
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          double delta_ij = (i == j) ? 1.0 : 0.0;
          sigma_avg(i, j) += (-p * delta_ij + mu * (ux(i, j) + ux(j, i)))
                             / static_cast<double>(lM.nG);
        }
      }
    }

    // Integrate over FACE Gauss points to get consistent nodal forces
    for (int gf = 0; gf < lFa.nG; gf++) {
      // Compute face normal vector (weighted by face Jacobian)
      Vector<double> nV(nsd);
      auto face_Nx = lFa.Nx.slice(gf);
      nn::gnnb(com_mod, lFa, e, gf, nsd, nsd - 1, lFa.eNoN, face_Nx, nV,
               solutions, consts::MechanicalConfigurationType::reference);
      double Jac_face = sqrt(utils::norm(nV));
      double w = lFa.w(gf) * Jac_face;

      // Unit normal
      for (int i = 0; i < nsd; i++) {
        nV(i) /= Jac_face;
      }

      // Traction: t_i = sigma_ij * n_j
      Vector<double> trac(nsd);
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          trac(i) += sigma_avg(i, j) * nV(j);
        }
      }

      // Accumulate consistent nodal forces at face nodes
      // Sign: solid sees NEGATIVE of fluid outward stress
      auto N = lFa.N.col(gf);
      for (int a = 0; a < lFa.eNoN; a++) {
        int Ac = lFa.IEN(a, e);      // global node ID
        int a_local = global_to_face[Ac];  // face-local node index
        for (int i = 0; i < nsd; i++) {
          result(i, a_local) -= w * N(a) * trac(i);
        }
      }
    }
  }

  // Sum contributions across MPI processes
  // Note: result is indexed by face-local nodes; need to communicate via global arrays
  // For now, use a global temporary array for communication
  Array<double> gResult(nsd, com_mod.tnNo);
  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    for (int i = 0; i < nsd; i++) {
      gResult(i, Ac) = result(i, a);
    }
  }
  all_fun::commu(com_mod, gResult);

  // Copy back to face-local result
  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    for (int i = 0; i < nsd; i++) {
      result(i, a) = gResult(i, Ac);
    }
  }

  return result;
}

//----------------------------------------------------------------------
// extract_solid_displacement
//----------------------------------------------------------------------
Array<double> extract_solid_displacement(
    const ComMod& com_mod, const eqType& solid_eq,
    const faceType& lFa, const SolutionStates& solutions)
{
  const int nsd = com_mod.nsd;
  const int s = solid_eq.s;  // DOF offset for the solid equation
  const auto& Dn = solutions.current.get_displacement();

  Array<double> result(nsd, lFa.nNo);
  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    for (int i = 0; i < nsd; i++) {
      result(i, a) = Dn(i + s, Ac);
    }
  }
  return result;
}

//----------------------------------------------------------------------
// extract_solid_velocity
//----------------------------------------------------------------------
Array<double> extract_solid_velocity(
    const ComMod& com_mod, const eqType& solid_eq,
    const faceType& lFa, const SolutionStates& solutions)
{
  const int nsd = com_mod.nsd;
  const int s = solid_eq.s;
  const auto& Yn = solutions.current.get_velocity();

  Array<double> result(nsd, lFa.nNo);
  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    for (int i = 0; i < nsd; i++) {
      result(i, a) = Yn(i + s, Ac);
    }
  }
  return result;
}

//----------------------------------------------------------------------
// apply_velocity_on_fluid
//----------------------------------------------------------------------
void apply_velocity_on_fluid(
    ComMod& com_mod, const eqType& fluid_eq,
    const faceType& lFa,
    const Array<double>& velocity,
    SolutionStates& solutions)
{
  const int nsd = com_mod.nsd;
  const int s = fluid_eq.s;

  auto& An = solutions.current.get_acceleration();
  auto& Yn = solutions.current.get_velocity();

  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    for (int i = 0; i < nsd; i++) {
      Yn(i + s, Ac) = velocity(i, a);
      An(i + s, Ac) = 0.0;  // zero acceleration for prescribed velocity
    }
  }
}

//----------------------------------------------------------------------
// apply_traction_on_solid
//----------------------------------------------------------------------
void apply_traction_on_solid(
    ComMod& com_mod, const eqType& solid_eq,
    const faceType& lFa,
    const Array<double>& traction)
{
  // The traction array contains consistent nodal forces that are already
  // integrated over the face. Add them directly to the residual R.
  // R stores the RHS of the Newton system, indexed by (dof, global_node).
  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    for (int i = 0; i < traction.nrows(); i++) {
      com_mod.R(i, Ac) += traction(i, a);
    }
  }
}

//----------------------------------------------------------------------
// apply_displacement_on_mesh
//----------------------------------------------------------------------
void apply_displacement_on_mesh(
    ComMod& com_mod, const eqType& mesh_eq,
    const faceType& lFa,
    const Array<double>& displacement,
    SolutionStates& solutions)
{
  const int nsd = com_mod.nsd;
  const int s = mesh_eq.s;
  const double dt = com_mod.dt;
  const double gam = mesh_eq.gam;
  const double beta = mesh_eq.beta;

  auto& An = solutions.current.get_acceleration();
  auto& Yn = solutions.current.get_velocity();
  auto& Dn = solutions.current.get_displacement();
  const auto& Do = solutions.old.get_displacement();
  const auto& Yo = solutions.old.get_velocity();
  const auto& Ao = solutions.old.get_acceleration();

  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    for (int i = 0; i < nsd; i++) {
      double d_new = displacement(i, a);
      double d_old = Do(i + s, Ac);
      double v_old = Yo(i + s, Ac);
      double a_old = Ao(i + s, Ac);

      // Prescribe displacement
      Dn(i + s, Ac) = d_new;

      // Compute acceleration and velocity consistent with Newmark:
      //   Dn = Do + dt*Yo + dt^2*((0.5-beta)*Ao + beta*An)
      //   Yn = Yo + dt*((1-gamma)*Ao + gamma*An)
      double a_new = (d_new - d_old - dt * v_old) / (beta * dt * dt)
                   - (0.5 - beta) / beta * a_old;
      double v_new = v_old + dt * ((1.0 - gam) * a_old + gam * a_new);

      An(i + s, Ac) = a_new;
      Yn(i + s, Ac) = v_new;
    }
  }
}

//----------------------------------------------------------------------
// transfer_face_data
//----------------------------------------------------------------------
Array<double> transfer_face_data(
    const ComMod& com_mod,
    const faceType& source_face, const faceType& target_face,
    const Array<double>& source_data)
{
  const int nrows = source_data.nrows();

  // Build reverse map: global_node_id -> source face local index
  std::unordered_map<int, int> global_to_source;
  global_to_source.reserve(source_face.nNo);
  for (int a = 0; a < source_face.nNo; a++) {
    global_to_source[source_face.gN(a)] = a;
  }

  // Transfer: for each target node, find matching source node via shared global ID
  Array<double> result(nrows, target_face.nNo);
  for (int a = 0; a < target_face.nNo; a++) {
    int global_id = target_face.gN(a);
    auto it = global_to_source.find(global_id);
    if (it != global_to_source.end()) {
      result.set_col(a, source_data.col(it->second));
    }
    // If not found, the node is not on the projected interface -- leave as zero
  }

  return result;
}

//----------------------------------------------------------------------
// regularize_unassembled_nodes
//----------------------------------------------------------------------
void regularize_unassembled_nodes(ComMod& com_mod, const mshType& active_mesh)
{
  const auto& eq = com_mod.eq[com_mod.cEq];
  const int dof = eq.dof;
  const auto& rowPtr = com_mod.rowPtr;
  const auto& colPtr = com_mod.colPtr;
  auto& R = com_mod.R;
  auto& Val = com_mod.Val;

  // Mark nodes belonging to the active mesh
  std::vector<bool> is_active(com_mod.tnNo, false);
  for (int a = 0; a < active_mesh.nNo; a++) {
    is_active[active_mesh.gN(a)] = true;
  }

  // For inactive nodes: zero R, zero all Val entries, set diagonal to 1
  for (int Ac = 0; Ac < com_mod.tnNo; Ac++) {
    if (is_active[Ac]) continue;

    // Zero residual
    for (int i = 0; i < dof; i++) {
      R(i, Ac) = 0.0;
    }

    // Zero entire row in Val and set diagonal to identity
    for (int j = rowPtr(Ac); j <= rowPtr(Ac + 1) - 1; j++) {
      for (int iDof = 0; iDof < dof * dof; iDof++) {
        Val(iDof, j) = 0.0;
      }
      if (colPtr(j) == Ac) {
        for (int i = 0; i < dof; i++) {
          Val(i * dof + i, j) = 1.0;
        }
      }
    }
  }
}

//----------------------------------------------------------------------
// enforce_dirichlet_on_face
//----------------------------------------------------------------------
// Following the pattern of set_bc_undef_neu_l (set_bc.cpp:1882):
// zero residual R and diagonalize Val rows at face nodes.
//
void enforce_dirichlet_on_face(ComMod& com_mod, const faceType& lFa, int nsd)
{
  const auto& eq = com_mod.eq[com_mod.cEq];
  const int dof = eq.dof;
  const auto& rowPtr = com_mod.rowPtr;
  const auto& colPtr = com_mod.colPtr;
  auto& R = com_mod.R;
  auto& Val = com_mod.Val;

  for (int a = 0; a < lFa.nNo; a++) {
    int rowN = lFa.gN(a);

    // Zero residual for all DOFs of this equation at this node
    for (int i = 0; i < dof; i++) {
      R(i, rowN) = 0.0;
    }

    // Diagonalize the row in the system matrix:
    // set diagonal = 1, off-diagonal = 0 for this node's DOFs
    for (int j = rowPtr(rowN); j <= rowPtr(rowN + 1) - 1; j++) {
      int colN = colPtr(j);
      for (int iDof = 0; iDof < dof * dof; iDof++) {
        Val(iDof, j) = 0.0;
      }
      if (colN == rowN) {
        // Set diagonal entries to 1
        for (int i = 0; i < dof; i++) {
          Val(i * dof + i, j) = 1.0;
        }
      }
    }
  }
}

//----------------------------------------------------------------------
// enforce_dirichlet_dofs_on_face
//----------------------------------------------------------------------
// Selective DOF enforcement: only modifies DOFs in [dof_start, dof_start+num_dofs).
// Used for fluid interface where velocity DOFs (0..nsd-1) should be fixed
// but pressure DOF (nsd) should remain free.
//
void enforce_dirichlet_dofs_on_face(ComMod& com_mod, const faceType& lFa,
                                    int dof_start, int num_dofs)
{
  const auto& eq = com_mod.eq[com_mod.cEq];
  const int dof = eq.dof;
  const auto& rowPtr = com_mod.rowPtr;
  const auto& colPtr = com_mod.colPtr;
  auto& R = com_mod.R;
  auto& Val = com_mod.Val;

  for (int a = 0; a < lFa.nNo; a++) {
    int rowN = lFa.gN(a);

    // Zero residual for specified DOFs only
    for (int i = dof_start; i < dof_start + num_dofs; i++) {
      R(i, rowN) = 0.0;
    }

    // Modify Val: zero specified DOF rows and set their diagonal to 1
    for (int j = rowPtr(rowN); j <= rowPtr(rowN + 1) - 1; j++) {
      int colN = colPtr(j);
      // Zero the rows for specified DOFs
      for (int i = dof_start; i < dof_start + num_dofs; i++) {
        for (int k = 0; k < dof; k++) {
          Val(i * dof + k, j) = 0.0;
        }
      }
      // Set diagonal entries for specified DOFs
      if (colN == rowN) {
        for (int i = dof_start; i < dof_start + num_dofs; i++) {
          Val(i * dof + i, j) = 1.0;
        }
      }
    }
  }
}

} // namespace fsi_coupling
