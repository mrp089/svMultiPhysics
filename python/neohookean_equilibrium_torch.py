"""
neohookean_equilibrium_torch.py
===============================

A self-contained, **differentiable**, **vectorized** PyTorch implementation of
the core svMultiPhysics machinery for the static solid-mechanics problem:

  * tetrahedral mesh with linear (TET4) shape functions
  * follower pressure BC on one surface  (endocardium)
  * traction-free (homogeneous Neumann) surface (epicardium)
  * fixed (homogeneous Dirichlet, u = 0) surface (base)
  * static problem
  * (nearly) incompressible neo-Hookean material, ST91 volumetric penalty

It is meant to be dropped straight into a neural-network training loop: the
global equilibrium residual

    R(U, theta) = F_internal(U, theta) - F_external_pressure(U)

is differentiable w.r.t. both the nodal displacement field U (e.g. a network
output) and the material parameters theta = (C10, Kpen). At equilibrium R = 0, so

    L_eq(U, theta) = || R_free(U, theta) ||^2

is a mechanical-equilibrium loss.

The formulas mirror the svMultiPhysics C++ source:
  * Quadrature / shape data .......... Code/Source/solver/nn_elem_gip.h, nn_elem_gnn.h
  * Parametric -> physical mapping ... Code/Source/solver/nn.cpp  (nn::gnn)
  * Element residual + tangent ....... Code/Source/solver/sv_struct.cpp (struct_3d)
  * Neo-Hookean stress + elasticity .. Code/Source/solver/mat_models.cpp
                                       (compute_pk2cc, bar_to_iso, compute_svol_p)
  * Follower pressure load ........... Code/Source/solver/sv_struct.cpp (b_struct_3d),
                                       Code/Source/solver/eq_assem.cpp  (b_neu_folw_p)

What's where
------------
  pk2_stress()                 batched neo-Hookean 2nd PK stress S(F)
  equilibrium_residual()       vectorized global residual R(U, theta)  (the loss core)
  equilibrium_loss()           ||R_free||^2
  assemble_sparse_tangent()    the SPARSE global Jacobian K = dR/dU (+ residual)
  solve_forward()              Newton solver using the sparse Jacobian
  TorchMesh / box_mesh()       mesh container + demo mesh generator

Sparse Jacobian
---------------
The tangent is assembled element-by-element: each element/face tangent block is
computed with autograd (torch.func.jacrev + vmap over the element batch), then
the small dense blocks are scattered into a single torch.sparse_coo_tensor. This
is O(nElements) storage/work, unlike a dense d R/d U which is O(nDof^2). The
follower-pressure tangent is unsymmetric (as svMultiPhysics notes), so the linear
solve uses a sparse *direct* LU. PyTorch's torch.sparse.spsolve is CUDA-only
(cuDSS), so on CPU we use SciPy's sparse LU as the linear-algebra backend; the
Jacobian itself is a torch sparse tensor either way.

Faithfulness note (Gauss loop)
------------------------------
For a *linear* tetrahedron F is constant over the element, so S, P = F.S and the
internal-force integrand are constant; the 4-point Gauss sum collapses to
(element volume) x integrand with element volume = Jac/6 = sum_g w_g * Jac.
Likewise the follower-pressure area-normal g1 x g2 is constant over a linear
triangle and sum_g w_g N_a = 1/6 per node. We use these closed forms (the exact
value of the Gauss sum, not an approximation) so assembly is a few einsums with
no Python-level Gauss loop. For higher-order elements the loop must be restored.
"""

from __future__ import annotations

import warnings

import numpy as np
import torch
from torch.func import vmap, jacrev

# torch sparse CSR/COO support is "beta" and prints a UserWarning on every tensor
# construction; the usage here is on well-formed matrices, so silence the noise.
warnings.filterwarnings("ignore", message=".*Sparse CSR tensor support is in beta.*")
warnings.filterwarnings("ignore", message=".*Sparse invariant checks.*")

NSD = 3

# Constant parametric shape derivatives dN_a/dxi_j (nn_elem_gnn.h).
_TET4_DNDXI = torch.tensor(
    [[1.0, 0.0, 0.0, -1.0],
     [0.0, 1.0, 0.0, -1.0],
     [0.0, 0.0, 1.0, -1.0]], dtype=torch.float64)          # (3, 4)

_TRI3_DNDXI = torch.tensor(
    [[1.0, 0.0, -1.0],
     [0.0, 1.0, -1.0]], dtype=torch.float64)                # (2, 3)


# =============================================================================
# Mesh container
# =============================================================================


class TorchMesh:
    """Connectivity + reference geometry as torch tensors, with per-element
    physical shape derivatives Nx and volumes precomputed (these depend only on
    the reference configuration).

    Parameters
    ----------
    nodes      : (nNo, 3) reference coordinates.
    tets       : (nEl, 4) TET4 connectivity (positively oriented).
    endo_faces : (nF, 3)  TRI3 pressure-surface connectivity, wound so
                 cross(g1, g2) points OUT of the solid.
    base_nodes : (nB,)    fixed nodes (homogeneous Dirichlet, all 3 comps).
    """

    def __init__(self, nodes, tets, endo_faces, base_nodes,
                 dtype=torch.float64, device="cpu"):
        self.dtype, self.device = dtype, device
        self.nodes = torch.as_tensor(nodes, dtype=dtype, device=device)
        self.tets = torch.as_tensor(tets, dtype=torch.long, device=device)
        self.endo_faces = torch.as_tensor(np.asarray(endo_faces), dtype=torch.long,
                                          device=device)
        self.base_nodes = torch.as_tensor(np.asarray(base_nodes), dtype=torch.long,
                                          device=device)
        self.nNo = self.nodes.shape[0]
        self.ndof = NSD * self.nNo

        self._precompute_volume_geometry()
        self._build_dof_maps()

    def _precompute_volume_geometry(self):
        """Per-tet physical derivatives Nx (nEl,3,4) and volume vol (nEl,).
        Mirrors nn::gnn: xXi(i,j)=sum_a X(a,i) dN_a/dxi_j ; Nx=(xXi^-1)^T dN/dxi ;
        vol = det(xXi)/6 (= sum of the four Gauss weights 1/24 times Jac)."""
        dN = _TET4_DNDXI.to(self.dtype).to(self.device)        # (3,4)
        Xe = self.nodes[self.tets]                             # (nEl,4,3)
        Jac = torch.linalg.det(torch.einsum("eai,ja->eij", Xe, dN))
        # Auto-reorient inverted tets (Jac < 0) by swapping two local nodes;
        # this only reorders connectivity, not the node-indexed fields, so it is
        # safe and necessary because the signed weight w = w_g * Jac must be > 0.
        inv = Jac < 0
        if inv.any():
            self.tets[inv] = self.tets[inv][:, [0, 1, 3, 2]]
            Xe = self.nodes[self.tets]
        xXi = torch.einsum("eai,ja->eij", Xe, dN)              # (nEl,3,3)
        Jac = torch.linalg.det(xXi)
        xiX = torch.linalg.inv(xXi)
        self.Nx = torch.einsum("eji,ja->eia", xiX, dN)         # (nEl,3,4)
        self.vol = Jac / 6.0
        if (self.vol <= 0).any():
            raise RuntimeError("TorchMesh: degenerate element (zero volume).")
        self.Xf = self.nodes[self.endo_faces]                  # (nF,3,3) ref coords

    def _build_dof_maps(self):
        fixed = torch.zeros(self.ndof, dtype=torch.bool, device=self.device)
        for c in range(NSD):
            fixed[NSD * self.base_nodes + c] = True
        self.fixed_mask = fixed
        self.free_mask = ~fixed
        self.free_idx = torch.nonzero(self.free_mask, as_tuple=False).flatten()
        # global dof -> free dof index, or -1 for fixed dofs
        g2f = torch.full((self.ndof,), -1, dtype=torch.long, device=self.device)
        g2f[self.free_idx] = torch.arange(self.free_idx.numel(), device=self.device)
        self.g2f = g2f
        # global dof index per (element, local node, component): (nEl, 12)
        self.elem_dofs = (NSD * self.tets[:, :, None]
                          + torch.arange(NSD, device=self.device)).reshape(-1, 12)
        if self.endo_faces.numel() > 0:
            self.face_dofs = (NSD * self.endo_faces[:, :, None]
                              + torch.arange(NSD, device=self.device)).reshape(-1, 9)
        else:
            self.face_dofs = torch.zeros((0, 9), dtype=torch.long, device=self.device)


# =============================================================================
# Constitutive law: batched neo-Hookean 2nd PK stress (compute_pk2cc)
# =============================================================================


def pk2_stress(F, C10, Kpen, vol_model="ST91"):
    """Batched neo-Hookean 2nd Piola-Kirchhoff stress S for F of shape (...,3,3).

    Mirrors mat_models.cpp compute_pk2cc (stIso_nHook) + bar_to_iso +
    compute_svol_p. Only S is returned; the tangent comes from autograd. C10/Kpen
    may be floats or tensors broadcastable to the batch (so a network can predict
    spatially varying material fields).
    """
    nsd = NSD
    I = torch.eye(nsd, dtype=F.dtype, device=F.device)
    J = torch.linalg.det(F)                                  # (...,)
    C = F.transpose(-1, -2) @ F                              # (...,3,3)
    Ci = torch.linalg.inv(C)
    J2d = J ** (-2.0 / nsd)

    C10 = torch.as_tensor(C10, dtype=F.dtype, device=F.device)
    Kpen = torch.as_tensor(Kpen, dtype=F.dtype, device=F.device)

    # ---- isochoric part: S_bar = 2 C10 I, CC_bar = 0 -> bar_to_iso ----------
    trC = C.diagonal(dim1=-2, dim2=-1).sum(-1)
    S_bar_scale = 2.0 * C10
    r1 = J2d * (S_bar_scale * trC) / nsd
    # S_iso = J2d S_bar - r1 Ci = (J2d 2 C10) I - r1 Ci
    S = (J2d * S_bar_scale)[..., None, None] * I - r1[..., None, None] * Ci

    # ---- volumetric penalty (compute_svol_p) -------------------------------
    if vol_model == "ST91":
        p = 0.5 * Kpen * (J - 1.0 / J)
    elif vol_model == "Quad":
        p = Kpen * (J - 1.0)
    elif vol_model == "M94":
        p = Kpen * (1.0 - 1.0 / J)
    else:
        raise ValueError(f"unknown vol_model {vol_model!r}")
    S = S + (p * J)[..., None, None] * Ci
    return S


# =============================================================================
# Per-element / per-face force kernels (single element; batched via vmap)
# =============================================================================


def _tet_internal_force(Ue, Nx, vol, C10, Kpen, vol_model):
    """Internal nodal force of ONE tet. Ue:(4,3) disp, Nx:(3,4), vol scalar.
    Returns f:(4,3). (Gauss sum collapsed: see module docstring.)"""
    I = torch.eye(NSD, dtype=Ue.dtype, device=Ue.device)
    grad_u = torch.einsum("ja,ai->ij", Nx, Ue)              # (3,3)
    F = I + grad_u
    S = pk2_stress(F, C10, Kpen, vol_model)                 # (3,3)
    P = F @ S                                                # 1st PK
    return vol * torch.einsum("ik,ka->ai", P, Nx)           # (4,3)


def _face_pressure_force(Uf, Xf, pressure):
    """Follower-pressure residual contribution of ONE face. Uf:(3,3) disp,
    Xf:(3,3) ref coords. Returns f:(3,3) (same +p/6 * (g1 x g2) on each node).

    This is the residual contribution (assemble adds it to R). It equals
    +p * (sum_g w_g N_a) * (g1 x g2) with sum_g w_g N_a = 1/6 per node, and
    g1, g2 the in-plane deformed tangents -- the surface form of Nanson's
    follower load used in b_struct_3d."""
    tri0 = _TRI3_DNDXI[0].to(Uf.dtype).to(Uf.device)
    tri1 = _TRI3_DNDXI[1].to(Uf.dtype).to(Uf.device)
    x = Xf + Uf                                             # current coords (3,3)
    g1 = torch.einsum("a,ac->c", tri0, x)
    g2 = torch.einsum("a,ac->c", tri1, x)
    nA = torch.linalg.cross(g1, g2)                         # area-normal
    f_node = (pressure / 6.0) * nA                          # (3,)
    return f_node.unsqueeze(0).expand(3, NSD)              # (3,3)


# =============================================================================
# Global equilibrium residual  R(U, theta)   (vectorized, differentiable)
# =============================================================================


def equilibrium_residual(mesh: TorchMesh, U, C10, Kpen, pressure, vol_model="ST91"):
    """Global residual R = F_internal - F_external_pressure, shape (nNo, 3).
    Differentiable w.r.t. U, C10, Kpen, pressure."""
    U = U.reshape(mesh.nNo, NSD)
    R = torch.zeros(mesh.nNo, NSD, dtype=U.dtype, device=U.device)

    # ---- internal force (vectorized over tets) -----------------------------
    Ue = U[mesh.tets]                                       # (nEl,4,3)
    Nx = mesh.Nx
    grad_u = torch.einsum("eja,eai->eij", Nx, Ue)
    F = torch.eye(NSD, dtype=U.dtype, device=U.device) + grad_u
    S = pk2_stress(F, C10, Kpen, vol_model)
    P = F @ S
    f_int = mesh.vol[:, None, None] * torch.einsum("eik,eka->eai", P, Nx)
    R = R.index_add(0, mesh.tets.reshape(-1), f_int.reshape(-1, NSD))

    # ---- follower pressure (vectorized over faces) -------------------------
    if mesh.endo_faces.numel() > 0:
        tri0 = _TRI3_DNDXI[0].to(U.dtype).to(U.device)
        tri1 = _TRI3_DNDXI[1].to(U.dtype).to(U.device)
        x = (mesh.nodes + U)[mesh.endo_faces]              # (nF,3,3)
        g1 = torch.einsum("a,fac->fc", tri0, x)
        g2 = torch.einsum("a,fac->fc", tri1, x)
        nA = torch.linalg.cross(g1, g2)
        f_pres = ((pressure / 6.0) * nA)[:, None, :].expand(-1, 3, -1)
        R = R.index_add(0, mesh.endo_faces.reshape(-1), f_pres.reshape(-1, NSD))
    return R


def equilibrium_loss(mesh: TorchMesh, U, C10, Kpen, pressure, vol_model="ST91",
                     reduction="mean"):
    """Mechanical-equilibrium loss L = ||R_free||^2 (free = non-Dirichlet dofs).
    Differentiable w.r.t. the network output U and/or material params C10, Kpen."""
    R = equilibrium_residual(mesh, U, C10, Kpen, pressure, vol_model).reshape(-1)
    Rf = R[mesh.free_mask]
    return Rf.pow(2).mean() if reduction == "mean" else Rf.pow(2).sum()


# =============================================================================
# Sparse global Jacobian  K = dR/dU
# =============================================================================


def assemble_sparse_tangent(mesh: TorchMesh, U, C10, Kpen, pressure,
                            vol_model="ST91", restrict_free=True):
    """Assemble the global residual and the SPARSE consistent tangent K = dR/dU.

    Each element/face tangent block is computed by autograd (jacrev) and the
    blocks are scattered into a torch.sparse_coo_tensor -- O(nElements) work and
    storage, never forming the dense nDof x nDof Jacobian.

    Parameters
    ----------
    restrict_free : if True (default) return K over the free dofs only
                    (nFree x nFree) and R over the free dofs (nFree,), ready for
                    the Newton solve. If False, return the full (ndof x ndof) /
                    (ndof,) system.

    Returns
    -------
    R : residual vector (free or full).
    K : torch sparse COO tangent (coalesced).
    """
    U = U.reshape(mesh.nNo, NSD)
    R_full = equilibrium_residual(mesh, U, C10, Kpen, pressure, vol_model).reshape(-1)

    rows_list, cols_list, vals_list = [], [], []

    # ---- element tangents: Ke[e] = d f_int(Ue) / d Ue  -> (nEl,4,3,4,3) ------
    Ue = U[mesh.tets]                                       # (nEl,4,3)
    tet_jac = jacrev(lambda u, nx, v: _tet_internal_force(
        u, nx, v, C10, Kpen, vol_model), argnums=0)
    Ke = vmap(tet_jac, in_dims=(0, 0, 0))(Ue, mesh.Nx, mesh.vol)
    Ke = Ke.reshape(-1, 12, 12)                            # (nEl,12,12)
    gd = mesh.elem_dofs                                    # (nEl,12)
    rows_list.append(gd[:, :, None].expand(-1, 12, 12).reshape(-1))
    cols_list.append(gd[:, None, :].expand(-1, 12, 12).reshape(-1))
    vals_list.append(Ke.reshape(-1))

    # ---- face tangents (follower pressure; unsymmetric) ---------------------
    if mesh.endo_faces.numel() > 0:
        Uf = U[mesh.endo_faces]                            # (nF,3,3)
        face_jac = jacrev(lambda u, xf: _face_pressure_force(u, xf, pressure),
                          argnums=0)
        Kf = vmap(face_jac, in_dims=(0, 0))(Uf, mesh.Xf)
        Kf = Kf.reshape(-1, 9, 9)                          # (nF,9,9)
        fd = mesh.face_dofs                                # (nF,9)
        rows_list.append(fd[:, :, None].expand(-1, 9, 9).reshape(-1))
        cols_list.append(fd[:, None, :].expand(-1, 9, 9).reshape(-1))
        vals_list.append(Kf.reshape(-1))

    rows = torch.cat(rows_list)
    cols = torch.cat(cols_list)
    vals = torch.cat(vals_list)

    if not restrict_free:
        K = torch.sparse_coo_tensor(torch.stack([rows, cols]), vals,
                                    (mesh.ndof, mesh.ndof)).coalesce()
        return R_full, K

    # Restrict to free dofs: drop any triplet touching a fixed dof, remap indices.
    fr = mesh.g2f[rows]
    fc = mesh.g2f[cols]
    keep = (fr >= 0) & (fc >= 0)
    n_free = mesh.free_idx.numel()
    K = torch.sparse_coo_tensor(torch.stack([fr[keep], fc[keep]]), vals[keep],
                                (n_free, n_free)).coalesce()
    return R_full[mesh.free_mask], K


def _sparse_lu_solve(K_coo, b):
    """Solve K x = b for a torch sparse-COO K using SciPy's sparse LU (CPU
    direct solver; handles the unsymmetric follower-pressure tangent). Returns a
    torch tensor. Used only inside the (no_grad) Newton iterations."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    K = K_coo.coalesce()
    idx = K.indices().cpu().numpy()
    val = K.values().detach().cpu().numpy()
    n = K.shape[0]
    A = sp.coo_matrix((val, (idx[0], idx[1])), shape=(n, n)).tocsc()
    x = spla.spsolve(A, b.detach().cpu().numpy())
    return torch.as_tensor(x, dtype=b.dtype, device=b.device)


# =============================================================================
# Newton forward solver (sparse Jacobian)
# =============================================================================


def solve_forward(mesh: TorchMesh, C10, Kpen, pressure, vol_model="ST91",
                  n_load_steps=5, max_newton=50, tol=1e-9, verbose=False):
    """Solve R(U)=0 for the displacement field by Newton's method using the
    sparse assembled Jacobian. Homogeneous Dirichlet dofs are eliminated;
    pressure is ramped over `n_load_steps` increments. Returns U:(nNo,3)."""
    free = mesh.free_idx
    U = torch.zeros(mesh.nNo, NSD, dtype=mesh.dtype, device=mesh.device)

    for ls in range(1, n_load_steps + 1):
        p_ls = pressure * ls / n_load_steps
        for it in range(max_newton):
            # NB: do NOT wrap in torch.no_grad() -- torch.func.jacrev (used inside
            # assemble_sparse_tangent) needs grad tracking and silently returns
            # zero derivatives under no_grad. U has requires_grad=False, so no
            # autograd graph is retained across iterations anyway.
            Rf, Kff = assemble_sparse_tangent(mesh, U, C10, Kpen, p_ls,
                                              vol_model, restrict_free=True)
            rnorm = torch.linalg.norm(Rf).item()
            if verbose:
                print(f"  load {ls}/{n_load_steps}  newton {it}  "
                      f"|R|={rnorm:.3e}  nnz={Kff._nnz()}")
            if rnorm < tol:
                break
            du = _sparse_lu_solve(Kff, -Rf)
            with torch.no_grad():
                U.reshape(-1)[free] += du
            if torch.linalg.norm(du).item() < tol:
                break
        else:
            raise RuntimeError(f"Newton did not converge at load step {ls}.")
    return U


# =============================================================================
# Demo mesh generator
# =============================================================================


def box_mesh(n=2, L=1.0):
    """Structured n x n x n unit-length-L box, split into 6 positively-oriented
    TET4 per voxel. Returns (nodes, tets, endo_faces, base_nodes) as numpy arrays:
        endo_faces : the x = 0 face (pressure surface), normals out (-x).
        base_nodes : the z = 0 face (fixed). Remaining faces are traction-free.
    """
    coords = np.linspace(0.0, L, n + 1)
    nodes = np.array([(x, y, z) for z in coords for y in coords for x in coords],
                     dtype=float)

    def nid(i, j, k):
        return i + (n + 1) * (j + (n + 1) * k)

    HEX_TETS = [(0, 1, 3, 7), (0, 1, 7, 5), (0, 5, 7, 4),
                (0, 3, 2, 7), (0, 2, 6, 7), (0, 6, 4, 7)]
    tets = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                c = [nid(i + dx, j + dy, k + dz)
                     for dz in (0, 1) for dy in (0, 1) for dx in (0, 1)]
                for a, b, cc, d in HEX_TETS:
                    tets.append((c[a], c[b], c[cc], c[d]))
    tets = np.array(tets, dtype=int)

    # Enforce positive Jacobian (det(dX/dxi) > 0): swap last two local nodes of
    # any inverted element. svMultiPhysics expects positively oriented elements.
    dN = _TET4_DNDXI.numpy()
    for e in range(tets.shape[0]):
        xl = nodes[tets[e]].T
        if np.linalg.det(xl @ dN.T) < 0.0:
            tets[e, [2, 3]] = tets[e, [3, 2]]

    # x = 0 face -> two triangles per cell, wound so cross(g1,g2) points -x (out).
    endo = []
    for k in range(n):
        for j in range(n):
            n00 = nid(0, j, k); n10 = nid(0, j + 1, k)
            n01 = nid(0, j, k + 1); n11 = nid(0, j + 1, k + 1)
            endo.append((n00, n01, n11))
            endo.append((n00, n11, n10))
    endo = np.array(endo, dtype=int)

    # Orient faces so cross(g1,g2) points away from the mesh centroid.
    tri = _TRI3_DNDXI.numpy()
    centroid = nodes.mean(axis=0)
    for f in range(endo.shape[0]):
        x = nodes[endo[f]]
        nrm = np.cross(tri[0] @ x, tri[1] @ x)
        if np.dot(nrm, x.mean(axis=0) - centroid) < 0.0:
            endo[f] = endo[f][[0, 2, 1]]

    base = np.array([nid(i, j, 0) for j in range(n + 1) for i in range(n + 1)],
                    dtype=int)
    return nodes, tets, endo, base


# =============================================================================
# Demo / self-test
# =============================================================================

if __name__ == "__main__":
    torch.set_printoptions(precision=4, sci_mode=True)
    dtype = torch.float64

    nodes, tets, endo, base = box_mesh(n=3, L=1.0)
    mesh = TorchMesh(nodes, tets, endo, base, dtype=dtype)
    C10_true, Kpen, pressure = 1.0e4, 1.0e6, 2.0e3
    print(f"Mesh: {mesh.nNo} nodes, {mesh.tets.shape[0]} tets, "
          f"{mesh.endo_faces.shape[0]} pressure faces, "
          f"{mesh.free_idx.numel()} free dofs.")

    # --- (a) sparse Jacobian vs dense autograd jacobian (validation) ----------
    from torch.autograd.functional import jacobian
    rng = np.random.default_rng(0)
    U0 = torch.zeros(mesh.nNo, NSD, dtype=dtype)
    U0.reshape(-1)[mesh.free_idx] = torch.tensor(
        1e-2 * rng.standard_normal(mesh.free_idx.numel()), dtype=dtype)

    Rf, Kff = assemble_sparse_tangent(mesh, U0, C10_true, Kpen, pressure)
    Kdense = Kff.to_dense()

    def res_free(uf):
        U = torch.zeros(mesh.ndof, dtype=dtype).index_copy(0, mesh.free_idx, uf)
        return equilibrium_residual(mesh, U, C10_true, Kpen, pressure).reshape(-1)[
            mesh.free_idx]
    Kref = jacobian(res_free, U0.reshape(-1)[mesh.free_idx], vectorize=True)
    rel = (Kdense - Kref).norm() / Kref.norm()
    print(f"(a) sparse tangent vs autograd jacobian, rel. error: {rel:.2e}  "
          f"(nnz={Kff._nnz()}, dense would be {mesh.free_idx.numel()**2} entries)")

    # --- (b) forward solve (Newton with sparse Jacobian) ----------------------
    print("\n(b) Forward solve (sparse Newton):")
    U = solve_forward(mesh, C10_true, Kpen, pressure, n_load_steps=5, verbose=True)
    print(f"    max |displacement| = {U.abs().max().item():.4e}")
    print(f"    equilibrium_loss at solution: "
          f"{equilibrium_loss(mesh, U, C10_true, Kpen, pressure).item():.3e}")

    # --- (c) physics loss differentiable w.r.t. U (network output) ------------
    print("\n(c) Physics loss is differentiable w.r.t. U:")
    U_net = U.clone().requires_grad_(True)
    equilibrium_loss(mesh, U_net, C10_true, Kpen, pressure).backward()
    print(f"    dL/dU computed, ||grad|| = {U_net.grad.norm().item():.3e}")

    # --- (d) inverse: identify C10 from displacements by gradient descent ------
    print("\n(d) Inverse: identify C10 from displacements (gradient descent):")
    log_C10 = torch.tensor(np.log(5.0e3), dtype=dtype, requires_grad=True)
    opt = torch.optim.Adam([log_C10], lr=0.05)
    U_meas = U.detach()
    for step in range(400):
        opt.zero_grad()
        loss = equilibrium_loss(mesh, U_meas, torch.exp(log_C10), Kpen, pressure)
        loss.backward()
        opt.step()
        if step % 80 == 0 or step == 399:
            print(f"    step {step:3d}  C10={torch.exp(log_C10).item():.4e}  "
                  f"loss={loss.item():.3e}")
    print(f"    true  C10 = {C10_true:.4e}")
    print(f"    found C10 = {torch.exp(log_C10).item():.4e}")
