/* Copyright (c) Stanford University, The Regents of the University of California, and others.
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
#include "lhsa.h"
#include "mat_fun.h"
#include "mat_fun_carray.h"
#include "mat_models.h"
#include "nn.h"
#include "utils.h"
#include "vtk_xml.h"
#include "gr_equilibrated.h"

#ifdef WITH_TRILINOS
#include "trilinos_linear_solver.h"
#endif

namespace gr {

/// @brief Loop solid elements and assemble into global matrices
void construct_gr_fd_ele(ComMod& com_mod, CepMod& cep_mod, CmMod& cm_mod, 
                      const mshType& lM, const Array<double>& Ag, 
                      const Array<double>& Yg, const Array<double>& Dg)
{
  using namespace consts;

  // Get dimensions
  const int eNoN = lM.eNoN;
  const int dof = com_mod.dof;

  // Initialize residual and tangent
  Vector<int> ptr(eNoN);
  Array<double> lR(dof, eNoN);
  Array3<double> lK(dof * dof, eNoN, eNoN);

  // Loop over all elements of mesh
  for (int e = 0; e < lM.nEl; e++) {
    // Reset
    ptr = 0;
    lR = 0.0;
    lK = 0.0;

    // Evaluate solid equations
    // eval_gr_fd_ele(e, com_mod, cep_mod, lM, Ag, Yg, Dg, ptr, lR, lK);
    eval_dsolid(e, com_mod, cep_mod, lM, Ag, Yg, Dg, ptr, lR, lK);

    // Assemble into global residual and tangent
    lhsa_ns::do_assem(com_mod, eNoN, ptr, lK, lR);
  }
}

/// @brief Finite difference on each element
void eval_gr_fd_ele(const int &e, ComMod &com_mod, CepMod &cep_mod,
                       const mshType &lM, const Array<double> &Ag,
                       const Array<double> &Yg, const Array<double> &Dg,
                       Vector<int> &ptr, Array<double> &lR, Array3<double> &lK,
                       const bool eval)
{
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
  auto& eq = com_mod.eq[cEq];
  const double dt = com_mod.dt;
  const double fd_eps = eq.af * eq.beta * dt * dt / eps;
  const double fy_eps = eq.af * eq.gam * dt / eps;
  const double fa_eps = eq.am / eps;

  // Make editable copy
  Array<double> e_Ag(tDof,tnNo); 
  Array<double> e_Yg(tDof,tnNo); 
  Array<double> e_Dg(tDof,tnNo);
  e_Ag = Ag;
  e_Yg = Yg;
  e_Dg = Dg;

  // Initialize residual and tangent
  Array<double> lRp(dof, eNoN);
  Array<double> dlR(dof, eNoN);
  Array3<double> lK_dummy(dof * dof, eNoN, eNoN);
  lRp = 0.0;
  lK_dummy = 0.0;

  // Central evaluation
  eval_dsolid(e, com_mod, cep_mod, lM, Ag, Yg, Dg, ptr, lR, lK_dummy);

  // Finite differences
  for (int i = 0; i < dof; ++i) {
    for (int a = 0; a < eNoN; ++a) {
      dlR = 0.0;
      
      // Global node ID
      Ac = lM.IEN(a, e);

      // Perturb
      e_Ag(i, Ac) += eps;
      e_Yg(i, Ac) += eps;
      e_Dg(i, Ac) += eps;

      // Aceleration
      lRp = 0.0;
      eval_dsolid(e, com_mod, cep_mod, lM, e_Ag, Yg, Dg, ptr, lRp, lK_dummy);
      dlR += (lRp - lR) * fa_eps;

      // Velocity
      lRp = 0.0;
      eval_dsolid(e, com_mod, cep_mod, lM, Ag, e_Yg, Dg, ptr, lRp, lK_dummy);
      dlR += (lRp - lR) * fy_eps;

      // Displacement
      lRp = 0.0;
      eval_dsolid(e, com_mod, cep_mod, lM, Ag, Yg, e_Dg, ptr, lRp, lK_dummy);
      dlR += (lRp - lR) * fd_eps;

      // Restore
      e_Ag(i, Ac) = Ag(i, Ac);
      e_Yg(i, Ac) = Yg(i, Ac);
      e_Dg(i, Ac) = Dg(i, Ac);

      // Assign to tangent matrix
      for (int j = 0; j < dof; ++j) {
        for (int b = 0; b < eNoN; ++b) {
          lK(j * dof + i, b, a) = dlR(j, b);
        }
      }
    }
  }
}

void construct_gr_fd_global(ComMod& com_mod, CepMod& cep_mod, CmMod& cm_mod, 
             const mshType& lM, const Array<double>& Ag, 
             const Array<double>& Yg, const Array<double>& Dg)
{
  // Get dimensions
  const int eNoN = lM.eNoN;
  const int dof = com_mod.dof;
  const int tDof = com_mod.tDof;
  const int tnNo = com_mod.tnNo;

  // Set step size for finite difference
  const double eps = 1.0e-8;

  // Time integration parameters
  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  const double dt = com_mod.dt;
  const double fd_eps = eq.af * eq.beta * dt * dt / eps;
  const double fy_eps = eq.af * eq.gam * dt / eps;
  const double fa_eps = eq.am / eps;

  // Make editable copy
  Array<double> e_Ag(tDof,tnNo); 
  Array<double> e_Yg(tDof,tnNo); 
  Array<double> e_Dg(tDof,tnNo);
  e_Ag = Ag;
  e_Yg = Yg;
  e_Dg = Dg;

  // Initialize residual and tangent
  Vector<int> ptr(eNoN);
  Array<double> lR(dof, eNoN);
  Array3<double> lK(dof * dof, eNoN, eNoN);

  // Central evaluation
  eval_gr_fd_global(com_mod, cep_mod, cm_mod, lM, Ag, Yg, Dg, fa_eps + fy_eps + fd_eps);

  // Loop global nodes
  for (int Ac = 0; Ac < tnNo; ++Ac) {
    // Central evaluation
    eval_gr_fd_global(com_mod, cep_mod, cm_mod, lM, Ag, Yg, Dg, fa_eps + fy_eps + fd_eps, Ac);

    // Loop DOFs
    for (int i = 0; i < dof; ++i) {
      // Perturb solution vectors
      e_Ag(i, Ac) += eps;
      e_Yg(i, Ac) += eps;
      e_Dg(i, Ac) += eps;

      // Perturbed evaluations
      eval_gr_fd_global(com_mod, cep_mod, cm_mod, lM, e_Ag, Yg, Dg, fa_eps, Ac, i);
      eval_gr_fd_global(com_mod, cep_mod, cm_mod, lM, Ag, e_Yg, Dg, fy_eps, Ac, i);
      eval_gr_fd_global(com_mod, cep_mod, cm_mod, lM, Ag, Yg, e_Dg, fd_eps, Ac, i);

      // Restore solution vectors
      e_Ag(i, Ac) = Ag(i, Ac);
      e_Yg(i, Ac) = Yg(i, Ac);
      e_Dg(i, Ac) = Dg(i, Ac);
    }
  }
// com_mod.Val.print("Val");
}

void eval_gr_fd_global(ComMod& com_mod, CepMod& cep_mod, CmMod& cm_mod, 
             const mshType& lM, const Array<double>& Ag, 
             const Array<double>& Yg, const Array<double>& Dg,
             const double eps, const int dAc, const int dj)
{
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

  // Smooth internal G&R variables
  enum Smoothing { none, element, elementnode };
  Smoothing smooth = elementnode;

  // Select element set to evaluate
  std::set<int> elements;

  // Residual evaluation: evaluate all elements
  if (residual) {
    for (int i = 0; i < lM.nEl; ++i) {
      elements.insert(i);
    }
  }

  // Pick elements according to smoothing algorithm
  else {
    switch(smooth) {
      case none: 
      case element: {
        elements = lM.map_node_ele[0].at(dAc);
        break;
      }
      case elementnode: {
        elements = lM.map_node_ele[1].at(dAc);
        break;
      }
    }
  }

  // Update internal G&R variables without assembly
  for (int e : elements) {
    eval_dsolid(e, com_mod, cep_mod, lM, Ag, Yg, Dg, ptr_dummy, lR_dummy, lK_dummy, false);
  }

  // Index of Lagrange multiplier
  std::vector<int> gr_variables = {37};

  switch(smooth) {
    // No smoothing
    case none: {
      break;
    }

    // Average over Gauss points in element
    case element: {
      for (int igr : gr_variables) {
        for (int e : elements) {
          double avg = 0.0;
          for (int g = 0; g < lM.nG; g++) {
            avg += com_mod.grInt(e, g, igr);
          }
          avg /= lM.nG;
          for (int g = 0; g < lM.nG; g++) {
            com_mod.grInt(e, g, igr) = avg;
          }
        }
      }
      break;
    }

    case elementnode: {
      // Initialize arrays
      // todo: this could be of size elements.size() * lM.eNoN
      Array<double> grInt_a(lM.gnNo, com_mod.nGrInt);
      Array<double> grInt_n(lM.gnNo, com_mod.nGrInt);
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
          for (int igr : gr_variables) {
            val = com_mod.grInt(e, g, igr);
            for (int a = 0; a < lM.eNoN; a++) {
              Ac = lM.IEN(a, e);
              // todo: add jacobian
              grInt_n(Ac, igr) += w * N(a) * val;
              grInt_a(Ac, igr) += w * N(a);
            }
          }
      }
      }

      // Project: nodes -> integration points
      for (int igr : gr_variables) {
        for (int e : elements) {
          for (int g = 0; g < lM.nG; g++) {
            N = lM.N.col(g);
            val = 0.0;
            for (int a = 0; a < lM.eNoN; a++) {
              Ac = lM.IEN(a, e);
              val += N(a) * grInt_n(Ac, igr) / grInt_a(Ac, igr);
            }
            com_mod.grInt(e, g, igr) = val;
          }
        }
      }
      break;
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
    eval_dsolid(e, com_mod, cep_mod, lM, Ag, Yg, Dg, ptr_row, lR, lK_dummy);

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
            lK(i * dof + j, a, 0) = - lR(i, a) * eps;
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
void eval_dsolid(const int &e, ComMod &com_mod, CepMod &cep_mod,
                 const mshType &lM, const Array<double> &Ag,
                 const Array<double> &Yg, const Array<double> &Dg,
                 Vector<int> &ptr, Array<double> &lR, Array3<double> &lK,
                 const bool eval)
{
  using namespace consts;

  #define n_debug_construct_dsolid
  #ifdef debug_construct_dsolid
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& cem = cep_mod.cem;
  const int nsd  = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int nsymd = com_mod.nsymd;
  auto& pS0 = com_mod.pS0;
  auto& pSn = com_mod.pSn;
  auto& pSa = com_mod.pSa;
  bool pstEq = com_mod.pstEq;

  int eNoN = lM.eNoN;
  int nFn = lM.nFn;
  if (nFn == 0) {
    nFn = 1;
  }

  #ifdef debug_construct_dsolid
  dmsg << "lM.nEl: " << lM.nEl;
  dmsg << "eNoN: " << eNoN;
  dmsg << "nsymd: " << nsymd;
  dmsg << "nFn: " << nFn;
  dmsg << "lM.nG: " << lM.nG;
  #endif

  // STRUCT: dof = nsd
  Vector<double> pSl(nsymd), ya_l(eNoN), N(eNoN), gr_int_g(com_mod.nGrInt), gr_props_g(lM.n_gr_props);
  Array<double> xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), 
                bfl(nsd,eNoN), fN(nsd,nFn), pS0l(nsymd,eNoN), Nx(nsd,eNoN),
                gr_props_l(lM.n_gr_props,eNoN);
  
  // Create local copies
  fN  = 0.0;
  pS0l = 0.0;
  ya_l = 0.0;
  gr_int_g = 0.0;
  gr_props_l = 0.0;
  gr_props_g = 0.0;

  for (int a = 0; a < eNoN; a++) {
    int Ac = lM.IEN(a,e);
    ptr(a) = Ac;

    for (int i = 0; i < nsd; i++) {
      xl(i,a) = com_mod.x(i,Ac);
      bfl(i,a) = com_mod.Bf(i,Ac);
    }

    for (int i = 0; i < tDof; i++) {
      al(i,a) = Ag(i,Ac);
      dl(i,a) = Dg(i,Ac);
      yl(i,a) = Yg(i,Ac);
    }

    if (lM.fN.size() != 0) {
      for (int iFn = 0; iFn < nFn; iFn++) {
        for (int i = 0; i < nsd; i++) {
          fN(i,iFn) = lM.fN(i+nsd*iFn,e);
        }
      }
    }

    if (pS0.size() != 0) { 
      pS0l.set_col(a, pS0.col(Ac));
    }

    if (cem.cpld) {
      ya_l(a) = cem.Ya(Ac);
    }

    if (lM.gr_props.size() != 0) {
      for (int igr = 0; igr < lM.n_gr_props; igr++) {
        gr_props_l(igr,a) = lM.gr_props(igr,Ac);
      }
    }
  }

  // Gauss integration
  double Jac{0.0};
  Array<double> ksix(nsd,nsd);

  for (int g = 0; g < lM.nG; g++) {
    if (g == 0 || !lM.lShpF) {
      auto Nx_g = lM.Nx.slice(g);
      nn::gnn(eNoN, nsd, nsd, Nx_g, xl, Nx, Jac, ksix);
      if (utils::is_zero(Jac)) {
        throw std::runtime_error("[construct_dsolid] Jacobian for element " + std::to_string(e) + " is < 0.");
      }
    }
    double w = lM.w(g) * Jac;
    N = lM.N.col(g);
    pSl = 0.0;

    // Get internal growth and remodeling variables
    if (com_mod.grEq) {
      // todo mrp089: add a function like rslice for vectors to Array3
      for (int i = 0; i < com_mod.nGrInt; i++) {
          gr_int_g(i) = com_mod.grInt(e,g,i);
      }
    }

    if (nsd == 3) {
      struct_3d_carray(com_mod, cep_mod, eNoN, nFn, w, N, Nx, al, yl, dl, bfl, fN, pS0l, pSl, ya_l, gr_int_g, gr_props_l, lR, lK, eval);
      // struct_3d_carray(com_mod, cep_mod, eNoN, nFn, w, N, Nx, al, yl, dl, bfl, fN, pS0l, pSl, ya_l, gr_int_g, gr_props_l, lR, lK);

#if 0
        if (e == 0 && g == 0) {
          Array3<double>::write_enabled = true;
          Array<double>::write_enabled = true;
          lR.write("lR");
          lK.write("lK");
          exit(0);
        }
#endif
    } else {
      std::terminate();
    }

    // Set internal growth and remodeling variables
    if (com_mod.grEq) {
      // todo mrp089: add a function like rslice for vectors to Array3
      for (int i = 0; i < com_mod.nGrInt; i++) {
          com_mod.grInt(e,g,i) = gr_int_g(i);
      }
    }

    // Prestress
    if (pstEq) {
      for (int a = 0; a < eNoN; a++) {
        int Ac = ptr(a);
        pSa(Ac) = pSa(Ac) + w*N(a);
        for (int i = 0; i < pSn.nrows(); i++) {
          pSn(i,Ac) = pSn(i,Ac) + w*N(a)*pSl(i);
        }
      }
    }
  } 
}

void struct_3d_carray(ComMod& com_mod, CepMod& cep_mod, const int eNoN, const int nFn, const double w, 
    const Vector<double>& N, const Array<double>& Nx, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& dl, const Array<double>& bfl, const Array<double>& fN, const Array<double>& pS0l, 
    Vector<double>& pSl, const Vector<double>& ya_l, Vector<double>& gr_int_g, Array<double>& gr_props_l, 
    Array<double>& lR, Array3<double>& lK, const bool eval) 
{
  using namespace consts;
  using namespace mat_fun;

  const int dof = com_mod.dof;
  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  // Set parameters
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  // Inertia, body force and deformation tensor (F)
  double F[3][3]={}; 
  Vector<double> gr_props_g(gr_props_l.nrows());

  F[0][0] = 1.0;
  F[1][1] = 1.0;
  F[2][2] = 1.0;
  double ya_g = 0.0;

  for (int a = 0; a < eNoN; a++) {
    F[0][0] += Nx(0,a)*dl(i,a);
    F[0][1] += Nx(1,a)*dl(i,a);
    F[0][2] += Nx(2,a)*dl(i,a);
    F[1][0] += Nx(0,a)*dl(j,a);
    F[1][1] += Nx(1,a)*dl(j,a);
    F[1][2] += Nx(2,a)*dl(j,a);
    F[2][0] += Nx(0,a)*dl(k,a);
    F[2][1] += Nx(1,a)*dl(k,a);
    F[2][2] += Nx(2,a)*dl(k,a);

    ya_g = ya_g + N(a)*ya_l(a);

    for (int igr = 0; igr < gr_props_l.nrows(); igr++) {
      gr_props_g(igr) += gr_props_l(igr,a) * N(a);
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
  double phic;
  get_pk2cc<3>(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, gr_int_g, gr_props_g, S, Dm, phic);

  // 1st Piola-Kirchhoff tensor (P)
  double P[3][3]; 
  mat_fun_carray::mat_mul<3>(F, S, P);

  // Local residual
  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) + w*(Nx(0,a)*P[0][0] + Nx(1,a)*P[0][1] + Nx(2,a)*P[0][2]);
    lR(1,a) = lR(1,a) + w*(Nx(0,a)*P[1][0] + Nx(1,a)*P[1][1] + Nx(2,a)*P[1][2]);
    lR(2,a) = lR(2,a) + w*(Nx(0,a)*P[2][0] + Nx(1,a)*P[2][1] + Nx(2,a)*P[2][2]);
  }

  // Auxilary quantities for computing stiffness tensor
  //
  Array3<double> Bm(6,3,eNoN);

  for (int a = 0; a < eNoN; a++) {
    Bm(0,0,a) = Nx(0,a)*F[0][0];
    Bm(0,1,a) = Nx(0,a)*F[1][0];
    Bm(0,2,a) = Nx(0,a)*F[2][0];

    Bm(1,0,a) = Nx(1,a)*F[0][1];
    Bm(1,1,a) = Nx(1,a)*F[1][1];
    Bm(1,2,a) = Nx(1,a)*F[2][1];

    Bm(2,0,a) = Nx(2,a)*F[0][2];
    Bm(2,1,a) = Nx(2,a)*F[1][2];
    Bm(2,2,a) = Nx(2,a)*F[2][2];

    Bm(3,0,a) = (Nx(0,a)*F[0][1] + F[0][0]*Nx(1,a));
    Bm(3,1,a) = (Nx(0,a)*F[1][1] + F[1][0]*Nx(1,a));
    Bm(3,2,a) = (Nx(0,a)*F[2][1] + F[2][0]*Nx(1,a));

    Bm(4,0,a) = (Nx(1,a)*F[0][2] + F[0][1]*Nx(2,a));
    Bm(4,1,a) = (Nx(1,a)*F[1][2] + F[1][1]*Nx(2,a));
    Bm(4,2,a) = (Nx(1,a)*F[2][2] + F[2][1]*Nx(2,a));

    Bm(5,0,a) = (Nx(2,a)*F[0][0] + F[0][2]*Nx(0,a));
    Bm(5,1,a) = (Nx(2,a)*F[1][0] + F[1][2]*Nx(0,a));
    Bm(5,2,a) = (Nx(2,a)*F[2][0] + F[2][2]*Nx(0,a));
  }

  // Local stiffness tensor
  double NxSNx, BmDBm;
  Array<double> DBm(6,3);

  for (int b = 0; b < eNoN; b++) {
    for (int a = 0; a < eNoN; a++) {
      // Geometric stiffness
      NxSNx = Nx(0,a)*S[0][0]*Nx(0,b) + Nx(1,a)*S[1][0]*Nx(0,b) +
              Nx(2,a)*S[2][0]*Nx(0,b) + Nx(0,a)*S[0][1]*Nx(1,b) +
              Nx(1,a)*S[1][1]*Nx(1,b) + Nx(2,a)*S[2][1]*Nx(1,b) +
              Nx(0,a)*S[0][2]*Nx(2,b) + Nx(1,a)*S[1][2]*Nx(2,b) +
              Nx(2,a)*S[2][2]*Nx(2,b);

      NxSNx = NxSNx;

      // Material Stiffness (Bt*D*B)
      mat_fun_carray::mat_mul6x3<3>(Dm, Bm.rslice(b), DBm);

      // dM1/du1
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,0,a)*DBm(0,0) + Bm(1,0,a)*DBm(1,0) +
              Bm(2,0,a)*DBm(2,0) + Bm(3,0,a)*DBm(3,0) +
              Bm(4,0,a)*DBm(4,0) + Bm(5,0,a)*DBm(5,0);

      lK(0,a,b) = lK(0,a,b) + w*(NxSNx + BmDBm);

      // dM1/du2
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,0,a)*DBm(0,1) + Bm(1,0,a)*DBm(1,1) +
              Bm(2,0,a)*DBm(2,1) + Bm(3,0,a)*DBm(3,1) +
              Bm(4,0,a)*DBm(4,1) + Bm(5,0,a)*DBm(5,1);

      lK(1,a,b) = lK(1,a,b) + w*BmDBm;

      // dM1/du3
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,0,a)*DBm(0,2) + Bm(1,0,a)*DBm(1,2) +
              Bm(2,0,a)*DBm(2,2) + Bm(3,0,a)*DBm(3,2) +
              Bm(4,0,a)*DBm(4,2) + Bm(5,0,a)*DBm(5,2);

      lK(2,a,b) = lK(2,a,b) + w*BmDBm;

      // dM2/du1
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,1,a)*DBm(0,0) + Bm(1,1,a)*DBm(1,0) +
              Bm(2,1,a)*DBm(2,0) + Bm(3,1,a)*DBm(3,0) +
              Bm(4,1,a)*DBm(4,0) + Bm(5,1,a)*DBm(5,0);

      lK(dof+0,a,b) = lK(dof+0,a,b) + w*BmDBm;

      // dM2/du2
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,1,a)*DBm(0,1) + Bm(1,1,a)*DBm(1,1) +
              Bm(2,1,a)*DBm(2,1) + Bm(3,1,a)*DBm(3,1) +
              Bm(4,1,a)*DBm(4,1) + Bm(5,1,a)*DBm(5,1);

      lK(dof+1,a,b) = lK(dof+1,a,b) + w*(NxSNx + BmDBm);

      // dM2/du3
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,1,a)*DBm(0,2) + Bm(1,1,a)*DBm(1,2) +
              Bm(2,1,a)*DBm(2,2) + Bm(3,1,a)*DBm(3,2) +
              Bm(4,1,a)*DBm(4,2) + Bm(5,1,a)*DBm(5,2);

      lK(dof+2,a,b) = lK(dof+2,a,b) + w*BmDBm;

      // dM3/du1
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,2,a)*DBm(0,0) + Bm(1,2,a)*DBm(1,0) +
              Bm(2,2,a)*DBm(2,0) + Bm(3,2,a)*DBm(3,0) +
              Bm(4,2,a)*DBm(4,0) + Bm(5,2,a)*DBm(5,0);

      lK(2*dof+0,a,b) = lK(2*dof+0,a,b) + w*BmDBm;
 
      // dM3/du2
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,2,a)*DBm(0,1) + Bm(1,2,a)*DBm(1,1) +
              Bm(2,2,a)*DBm(2,1) + Bm(3,2,a)*DBm(3,1) +
              Bm(4,2,a)*DBm(4,1) + Bm(5,2,a)*DBm(5,1);

     lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + w*BmDBm;

      // dM3/du3
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,2,a)*DBm(0,2) + Bm(1,2,a)*DBm(1,2) +
              Bm(2,2,a)*DBm(2,2) + Bm(3,2,a)*DBm(3,2) +
              Bm(4,2,a)*DBm(4,2) + Bm(5,2,a)*DBm(5,2);

      lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + w*(NxSNx + BmDBm);
    }
  }

}

};