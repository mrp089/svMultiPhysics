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
#include "mat_models_carray.h"
#include "nn.h"
#include "utils.h"
#include "vtk_xml.h"

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
    eval_gr_fd_ele(e, com_mod, cep_mod, lM, Ag, Yg, Dg, ptr, lR, lK);

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
  struct_ns::eval_dsolid(e, com_mod, cep_mod, lM, Ag, Yg, Dg, ptr, lR, lK_dummy);

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
      struct_ns::eval_dsolid(e, com_mod, cep_mod, lM, e_Ag, Yg, Dg, ptr, lRp, lK_dummy);
      dlR += (lRp - lR) * fa_eps;

      // Velocity
      lRp = 0.0;
      struct_ns::eval_dsolid(e, com_mod, cep_mod, lM, Ag, e_Yg, Dg, ptr, lRp, lK_dummy);
      dlR += (lRp - lR) * fy_eps;

      // Displacement
      lRp = 0.0;
      struct_ns::eval_dsolid(e, com_mod, cep_mod, lM, Ag, Yg, e_Dg, ptr, lRp, lK_dummy);
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
    struct_ns::eval_dsolid(e, com_mod, cep_mod, lM, Ag, Yg, Dg, ptr_dummy, lR_dummy, lK_dummy, false);
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
    struct_ns::eval_dsolid(e, com_mod, cep_mod, lM, Ag, Yg, Dg, ptr_row, lR, lK_dummy);

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

};