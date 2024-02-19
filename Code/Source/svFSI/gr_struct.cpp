/* Copyright (c) Stanford University, The Regents of the University of
 * California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// Functions for solving nonlinear structural mechanics
// problems (pure displacement-based formulation).
//
// Replicates the Fortran functions in 'STRUCT.f'.

#include "gr_struct.h"
#include "sv_struct.h"

#include "all_fun.h"
#include "consts.h"
#include "gr_equilibrated.h"
#include "lhsa.h"
#include "mat_fun.h"
#include "mat_fun_carray.h"
#include "mat_models.h"
#include "nn.h"
#include "utils.h"
#include "vtk_xml.h"

#ifdef WITH_TRILINOS
#include "trilinos_linear_solver.h"
#endif

namespace gr {

/// @brief Loop solid elements and assemble into global matrices
void construct_gr(ComMod &com_mod, const mshType &lM, const Array<double> &Dg,
                  bool eval_fd) {
  using namespace consts;

  // Get dimensions
  const int eNoN = lM.eNoN;
  const int dof = com_mod.dof;

  // Initialize residual and tangent
  Vector<int> ptr(eNoN);
  Array<double> lR(dof, eNoN);
  Array3<double> lK(dof * dof, eNoN, eNoN);

  // Make editable copy
  Array<double> e_Dg(Dg);

  // Loop over all elements of mesh
  for (int e = 0; e < lM.nEl; e++) {
    // Reset
    ptr = 0;
    lR = 0.0;
    lK = 0.0;

    // Evaluate solid equations
    if (eval_fd) {
      eval_gr_fd_ele(e, com_mod, lM, e_Dg, ptr, lR, lK);
    } else {
      // Update G&R internal variables
      eval_gr(e, com_mod, lM, Dg, ptr, lR, lK, false, false);

      // Compute stress and tangent
      eval_gr(e, com_mod, lM, Dg, ptr, lR, lK, true, true);
    }

    // Assemble into global residual and tangent
    lhsa_ns::do_assem(com_mod, eNoN, ptr, lK, lR);
  }
}

/// @brief Finite difference on each element
void eval_gr_fd_ele(const int &e, ComMod &com_mod, const mshType &lM,
                    Array<double> &Dg, Vector<int> &ptr, Array<double> &lR,
                    Array3<double> &lK) {
  // Get dimensions
  const int eNoN = lM.eNoN;
  const int dof = com_mod.dof;
  const int tDof = com_mod.tDof;
  const int tnNo = com_mod.tnNo;

  // Set step size for finite difference
  const double eps = 1.0e-10;

  // Time integration parameters
  int Ac;
  int cEq = com_mod.cEq;
  auto &eq = com_mod.eq[cEq];

  // Initialize residual and tangent
  Array<double> lRp(dof, eNoN);
  Array<double> dlR(dof, eNoN);
  Array3<double> lK_dummy(dof * dof, eNoN, eNoN);
  lRp = 0.0;
  lK_dummy = 0.0;

  // Central evaluation
  eval_gr(e, com_mod, lM, Dg, ptr, lR, lK_dummy, true, false);

  // Finite differences
  for (int i = 0; i < dof; ++i) {
    for (int a = 0; a < eNoN; ++a) {
      dlR = 0.0;

      // Global node ID
      Ac = lM.IEN(a, e);

      // Perturb
      Dg(i, Ac) += eps;

      // Displacement
      lRp = 0.0;
      eval_gr(e, com_mod, lM, Dg, ptr, lRp, lK_dummy, true, false);
      dlR += (lRp - lR) / eps;

      // Restore
      Dg(i, Ac) -= eps;

      // Assign to tangent matrix
      for (int j = 0; j < dof; ++j) {
        for (int b = 0; b < eNoN; ++b) {
          lK(j * dof + i, b, a) = dlR(j, b);
        }
      }
    }
  }
}

void construct_gr_fd_global(ComMod &com_mod, const mshType &lM,
                            const Array<double> &Dg) {
  // Get dimensions
  const int eNoN = lM.eNoN;
  const int dof = com_mod.dof;
  const int tDof = com_mod.tDof;
  const int tnNo = com_mod.tnNo;

  // Set step size for finite difference
  const double eps = 1.0e-8;

  // Make editable copy
  Array<double> e_Dg(tDof, tnNo);
  e_Dg = Dg;

  // Initialize residual and tangent
  Vector<int> ptr(eNoN);
  Array<double> lR(dof, eNoN);
  Array3<double> lK(dof * dof, eNoN, eNoN);

  // Residual evaluation
  eval_gr_fd_global(com_mod, lM, Dg, 1.0 / eps);

  // Loop global nodes
  for (int Ac = 0; Ac < tnNo; ++Ac) {
    // Central evaluation
    eval_gr_fd_global(com_mod, lM, Dg, 1.0 / eps, Ac);

    // Loop DOFs
    for (int i = 0; i < dof; ++i) {
      // Perturb solution vectors
      e_Dg(i, Ac) += eps;

      // Perturbed evaluations
      eval_gr_fd_global(com_mod, lM, e_Dg, 1.0 / eps, Ac, i);

      // Restore solution vectors
      e_Dg(i, Ac) = Dg(i, Ac);
    }
  }
}

void eval_gr_fd_global(ComMod &com_mod, const mshType &lM,
                       const Array<double> &Dg, const double eps, const int dAc,
                       const int dj) {
  using namespace consts;

  // Check if the residual should be assembled
  const bool residual = (dAc == -1) && (dj == -1);

  // Check if this is the central evaluation
  const bool central = (dAc != -1) && (dj == -1);

  // Get dimensions
  const int eNoN = lM.eNoN;
  const int dof = com_mod.dof;

  // Initialize residual and tangent
  Vector<int> ptr_dummy(eNoN);
  Array<double> lR_dummy(dof, eNoN);
  Array3<double> lK_dummy(dof * dof, eNoN, eNoN);
  ptr_dummy = 0;
  lR_dummy = 0.0;
  lK_dummy = 0.0;

  // Select element set to evaluate
  std::set<int> elements;
  if (residual) {
    // Residual evaluation: evaluate all elements
    for (int i = 0; i < lM.nEl; ++i) {
      elements.insert(i);
    }
  } else {
    // Pick elements according to smoothing algorithm
    elements = lM.map_node_ele[1].at(dAc);
  }

  // Update internal G&R variables without evaluating stress or tangent
  for (int e : elements) {
    eval_gr(e, com_mod, lM, Dg, ptr_dummy, lR_dummy, lK_dummy, false, false);
  }

  // Index of collagen mass fraction multiplier
  const int igr = 37;

  // Initialize arrays
  Vector<double> grInt_a(lM.gnNo);
  Vector<double> grInt_n(lM.gnNo);
  grInt_a = 0.0;
  grInt_n = 0.0;

  int Ac;
  double w;
  double val;
  Vector<double> N(lM.eNoN);

  // Project: integration point -> nodes
  for (int g = 0; g < lM.nG; g++) {
    w = lM.w(g);
    N = lM.N.col(g);
    for (int e : elements) {
      val = com_mod.grInt(e, g, igr);
      for (int a = 0; a < lM.eNoN; a++) {
        Ac = lM.IEN(a, e);
        // todo: add jacobian
        grInt_n(Ac) += w * N(a) * val;
        grInt_a(Ac) += w * N(a);
      }
    }
  }

  // Project: nodes -> integration points
  for (int e : elements) {
    for (int g = 0; g < lM.nG; g++) {
      N = lM.N.col(g);
      val = 0.0;
      for (int a = 0; a < lM.eNoN; a++) {
        Ac = lM.IEN(a, e);
        val += N(a) * grInt_n(Ac) / grInt_a(Ac);
      }
      com_mod.grInt(e, g, igr) = val;
    }
  }

  // Store internal G&R variables
  if (residual) {
    com_mod.grInt_orig = com_mod.grInt;
  }

  // Initialzie arrays for Finite Difference (FD)
  Vector<int> ptr_row(eNoN);
  Vector<int> ptr_col(1);
  Array<double> lR(dof, eNoN);
  Array3<double> lK(dof * dof, eNoN, 1);

  // Assemble only the FD node
  ptr_col = dAc;

  // Loop over all elements of mesh
  for (int e : elements) {
    // Reset
    ptr_row = 0;
    lR = 0.0;
    lK = 0.0;

    // Evaluate solid equations (with smoothed internal G&R variables)
    eval_gr(e, com_mod, lM, Dg, ptr_row, lR, lK_dummy, true, false);

    // Assemble into global residual
    if (residual) {
      lhsa_ns::do_assem_residual(com_mod, lM.eNoN, ptr_row, lR);
      continue;
    }

    // Components of FD: central and difference
    for (int a = 0; a < eNoN; ++a) {
      for (int i = 0; i < dof; ++i) {
        if (central) {
          for (int j = 0; j < dof; ++j) {
            lK(i * dof + j, a, 0) = -lR(i, a) * eps;
          }
        } else {
          lK(i * dof + dj, a, 0) = lR(i, a) * eps;
        }
      }
    }

    // Assemble into global tangent
    lhsa_ns::do_assem_tangent(com_mod, lM.eNoN, 1, ptr_row, ptr_col, lK);
  }

  // Restore internal G&R variables
  com_mod.grInt = com_mod.grInt_orig;
}

/// @brief
void eval_gr(const int &e, ComMod &com_mod, const mshType &lM,
             const Array<double> &Dg, Vector<int> &ptr, Array<double> &lR,
             Array3<double> &lK, const bool eval_s, const bool eval_cc) {
  using namespace consts;

  const int nsd = com_mod.nsd;
  const int tDof = com_mod.tDof;
  int eNoN = lM.eNoN;
  assert(nsd == 3);

  // STRUCT: dof = nsd
  Vector<double> N(eNoN), gr_int_g(com_mod.nGrInt), gr_props_g(lM.n_gr_props);
  Array<double> xl(nsd, eNoN), dl(tDof, eNoN), Nx(nsd, eNoN),
      gr_props_l(lM.n_gr_props, eNoN);

  // Create local copies
  gr_int_g = 0.0;
  gr_props_l = 0.0;
  gr_props_g = 0.0;

  for (int a = 0; a < eNoN; a++) {
    int Ac = lM.IEN(a, e);
    ptr(a) = Ac;

    for (int i = 0; i < nsd; i++) {
      xl(i, a) = com_mod.x(i, Ac);
    }

    for (int i = 0; i < tDof; i++) {
      dl(i, a) = Dg(i, Ac);
    }

    if (lM.gr_props.size() != 0) {
      for (int igr = 0; igr < lM.n_gr_props; igr++) {
        gr_props_l(igr, a) = lM.gr_props(igr, Ac);
      }
    }
  }

  // Gauss integration
  double Jac{0.0};
  Array<double> ksix(nsd, nsd);

  for (int g = 0; g < lM.nG; g++) {
    if (g == 0 || !lM.lShpF) {
      auto Nx_g = lM.Nx.slice(g);
      nn::gnn(eNoN, nsd, nsd, Nx_g, xl, Nx, Jac, ksix);
      if (utils::is_zero(Jac)) {
        throw std::runtime_error("[construct_dsolid] Jacobian for element " +
                                 std::to_string(e) + " is < 0.");
      }
    }
    double w = lM.w(g) * Jac;
    N = lM.N.col(g);

    // Get internal growth and remodeling variables
    if (com_mod.grEq) {
      for (int i = 0; i < com_mod.nGrInt; i++) {
        gr_int_g(i) = com_mod.grInt(e, g, i);
      }
    }

    struct_3d_carray(com_mod, eNoN, w, N, Nx, dl, gr_int_g, gr_props_l, lR, lK,
                     eval_s, eval_cc);

    // Set internal growth and remodeling variables
    if (com_mod.grEq) {
      // todo mrp089: add a function like rslice for vectors to Array3
      for (int i = 0; i < com_mod.nGrInt; i++) {
        com_mod.grInt(e, g, i) = gr_int_g(i);
      }
    }
  }
}

void struct_3d_carray(ComMod &com_mod, const int eNoN, const double w,
                      const Vector<double> &N, const Array<double> &Nx,
                      const Array<double> &dl, Vector<double> &gr_int_g,
                      Array<double> &gr_props_l, Array<double> &lR,
                      Array3<double> &lK, const bool eval_s,
                      const bool eval_cc) {
  using namespace consts;
  using namespace mat_fun;

  const int dof = com_mod.dof;
  int cEq = com_mod.cEq;
  auto &eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto &dmn = eq.dmn[cDmn];

  // Set parameters
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;
  int indices[] = {i, j, k};

  // Inertia, body force and deformation tensor (F)
  double F[3][3] = {};
  Vector<double> gr_props_g(gr_props_l.nrows());

  F[0][0] = 1.0;
  F[1][1] = 1.0;
  F[2][2] = 1.0;

  for (int a = 0; a < eNoN; a++) {
    // Deformation gradient
    for (int row = 0; row < 3; row++) {
      for (int col = 0; col < 3; col++) {
        F[row][col] += Nx(col, a) * dl(indices[row], a);
      }
    }

    // G&R material properties
    for (int igr = 0; igr < gr_props_l.nrows(); igr++) {
      gr_props_g(igr) += gr_props_l(igr, a) * N(a);
    }
  }

  double Jac = mat_fun_carray::mat_det<3>(F);
  double Fi[3][3];
  mat_fun_carray::mat_inv<3>(F, Fi);

  // Initialize tensor indexing.
  mat_fun_carray::ten_init(3);

  // 2nd Piola-Kirchhoff tensor (S) and material stiffness tensor in
  // Voigt notationa (Dm)
  double S[3][3];
  double Dm[6][6];
  get_pk2cc<3>(com_mod, dmn, F, gr_int_g, gr_props_g, S, Dm, eval_s, eval_cc);

  if (!eval_s && !eval_cc) {
    return;
  }

  // 1st Piola-Kirchhoff tensor (P)
  double P[3][3];
  mat_fun_carray::mat_mul<3>(F, S, P);

  // Local residual
  for (int a = 0; a < eNoN; a++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        lR(i, a) += w * Nx(j, a) * P[i][j];
      }
    }
  }

  // Skip evaluation if tangent not required
  if (!eval_cc) {
    return;
  }

  // Auxilary quantities for computing stiffness tensor
  Array3<double> Bm(6, 3, eNoN);
  for (int a = 0; a < eNoN; a++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        Bm(i, j, a) = Nx(i, a) * F[j][i];
      }
    }

    Bm(3, 0, a) = (Nx(0, a) * F[0][1] + F[0][0] * Nx(1, a));
    Bm(3, 1, a) = (Nx(0, a) * F[1][1] + F[1][0] * Nx(1, a));
    Bm(3, 2, a) = (Nx(0, a) * F[2][1] + F[2][0] * Nx(1, a));

    Bm(4, 0, a) = (Nx(1, a) * F[0][2] + F[0][1] * Nx(2, a));
    Bm(4, 1, a) = (Nx(1, a) * F[1][2] + F[1][1] * Nx(2, a));
    Bm(4, 2, a) = (Nx(1, a) * F[2][2] + F[2][1] * Nx(2, a));

    Bm(5, 0, a) = (Nx(2, a) * F[0][0] + F[0][2] * Nx(0, a));
    Bm(5, 1, a) = (Nx(2, a) * F[1][0] + F[1][2] * Nx(0, a));
    Bm(5, 2, a) = (Nx(2, a) * F[2][0] + F[2][2] * Nx(0, a));
  }

  // Local stiffness tensor
  double NxSNx, BmDBm;
  Array<double> DBm(6, 3);

  for (int b = 0; b < eNoN; b++) {
    for (int a = 0; a < eNoN; a++) {
      // Geometric stiffness
      NxSNx = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          NxSNx += Nx(i, a) * S[i][j] * Nx(j, b);
        }
      }

      // Material Stiffness (Bt*D*B)
      mat_fun_carray::mat_mul6x3<3>(Dm, Bm.rslice(b), DBm);
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          BmDBm = 0.0;
          for (int k = 0; k < 6; k++) {
            BmDBm += Bm(k, i, a) * DBm(k, j);
          }
          int index = i * dof + j;
          lK(index, a, b) += w * BmDBm;
          if (index == 0 || index == 4 || index == 8) {
            lK(index, a, b) += w * NxSNx;
          }
        }
      }
    }
  }
}

}; // namespace gr