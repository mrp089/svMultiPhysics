"""
tet10_identifiability.py
========================

Full synthetic loop testing whether QUADRATIC (TET10) elements make the
per-element neo-Hookean stiffness identifiable from DISPLACEMENTS ALONE
(no stress input -- a mockup of image-based material identification).

Pipeline
--------
  1. Build a TET10 box mesh (linear box -> add edge-midpoint nodes).
  2. Plant a KNOWN heterogeneous C10 field: healthy bulk + a stiff spherical
     "infarct" inclusion.
  3. Forward-solve static neo-Hookean equilibrium under a known (dead) pressure
     load with the planted field -> synthetic displacement field u_data.
  4. INVERSE: recover per-element C10 from u_data alone, with NO regularizer,
     by minimizing the equilibrium-residual ||R_free(u_data, c)||^2.
  5. DEGRADE the "imaging": replace the mid-edge displacements by linear
     interpolation of the corner displacements (i.e. only vertex-resolution data)
     and re-run the inverse -> identifiability should collapse.

This isolates the identifiability question from discretization error (forward and
inverse use the same TET10 discretization -- a consistent "inverse crime", which
is exactly what you want to test identifiability).

Material: neo-Hookean, nu = 0.48 given, Kpen = C10 * 4(1+nu)/(3(1-2nu)) per element
(so one unknown C10 per element). Quadrature & node ordering mirror svMultiPhysics
TET10 (15-pt Keast rule). Uses pk2_stress / _sparse_lu_solve from the TET4 module.
"""

import numpy as np
import torch
from torch.func import vmap, jacrev

from neohookean_equilibrium_torch import pk2_stress, _sparse_lu_solve, box_mesh, NSD

torch.set_default_dtype(torch.float64)
NU = 0.48
KFAC = 4.0 * (1.0 + NU) / (3.0 * (1.0 - 2.0 * NU))      # Kpen = KFAC * C10

# svMultiPhysics TET10 edges (local node 4..9), and shape functions ----------
EDGES = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]


def tet10_N(xi):
    """TET10 shape functions (svMultiPhysics ordering). xi=(r,s,t) -> (10,)."""
    r, s, t = xi[0], xi[1], xi[2]
    u = 1.0 - r - s - t
    return torch.stack([
        r * (2 * r - 1), s * (2 * s - 1), t * (2 * t - 1), u * (2 * u - 1),
        4 * r * s, 4 * s * t, 4 * r * t, 4 * r * u, 4 * s * u, 4 * t * u])


def tet10_quadrature():
    """15-point degree-5 Keast rule (svMultiPhysics nn_elem_gip.h TET10)."""
    w = torch.tensor([0.0302836780970890] + [0.0060267857142860] * 4
                     + [0.0116452490860290] * 4 + [0.0109491415613860] * 6)
    s0 = 0.25
    a, b = 1.0 / 3.0, 0.0
    c, d = 0.0909090909090910, 0.7272727272727270
    e, f = 0.0665501535736640, 0.4334498464263360
    xi = torch.tensor([
        [s0, s0, s0],
        [b, a, a], [a, b, a], [a, a, b], [a, a, a],
        [d, c, c], [c, d, c], [c, c, d], [c, c, c],
        [e, e, f], [e, f, e], [e, f, f], [f, f, e], [f, e, f], [f, e, e]])
    return w, xi


def tri6_surface_coeff():
    """3-pt rule integral coefficients c_a = sum_g w_g M_a(zeta_g) for the TRI6
    face, node order [c0,c1,c2, m01,m12,m20]. Consistent dead-pressure load on a
    flat face is f_a = p * |cross(e1,e2)| * c_a * n_hat."""
    zeta = torch.tensor([[2 / 3, 1 / 6], [1 / 6, 2 / 3], [1 / 6, 1 / 6]])
    wg = torch.full((3,), 1.0 / 6.0)
    coeff = torch.zeros(6)
    for g in range(3):
        r, s = zeta[g]; u = 1 - r - s
        M = torch.tensor([r * (2 * r - 1), s * (2 * s - 1), u * (2 * u - 1),
                          4 * r * s, 4 * s * u, 4 * u * r])
        coeff += wg[g] * M
    return coeff                                          # (6,)


# ---------------------------------------------------------------------------
# Mesh: TET4 box -> TET10
# ---------------------------------------------------------------------------


def build_tet10(nodes4, tets4):
    ev = np.stack([np.sort(tets4[:, list(e)], axis=1) for e in EDGES], axis=1)  # (nEl,6,2)
    uniq, inv = np.unique(ev.reshape(-1, 2), axis=0, return_inverse=True)
    mids = 0.5 * (nodes4[uniq[:, 0]] + nodes4[uniq[:, 1]])
    nV = nodes4.shape[0]
    edge_node = (nV + inv).reshape(-1, 6)
    nodes10 = np.vstack([nodes4, mids])
    tets10 = np.hstack([tets4, edge_node]).astype(np.int64)
    key = uniq[:, 0].astype(np.int64) * (nV + 1) + uniq[:, 1]
    lookup = {int(k): i + nV for i, k in enumerate(key)}     # edge -> mid node id
    return nodes10, tets10, nV, lookup


def tri6_faces(tris3, nV, lookup):
    out = []
    for a, b, c in tris3:
        def mid(x, y):
            x, y = (x, y) if x < y else (y, x)
            return lookup[int(x) * (nV + 1) + int(y)]
        out.append([a, b, c, mid(a, b), mid(b, c), mid(c, a)])
    return np.array(out, dtype=np.int64)


# ---------------------------------------------------------------------------
# Precompute geometry (reference config): Nx = dN/dX, w*detJ per Gauss point
# ---------------------------------------------------------------------------


class Tet10Mesh:
    def __init__(self, nodes10, tets10, base_nodes, faces6, p_unit_dir=(1.0, 0, 0)):
        self.nodes = torch.tensor(nodes10)
        self.tets = torch.tensor(tets10)
        self.nNo = nodes10.shape[0]
        self.nEl = tets10.shape[0]
        self.ndof = 3 * self.nNo

        wq, xiq = tet10_quadrature()
        self.nG = wq.shape[0]
        dNdxi = torch.stack([torch.autograd.functional.jacobian(tet10_N, xiq[g])
                             for g in range(self.nG)])        # (nG,10,3)

        Xe = self.nodes[self.tets]                            # (nEl,10,3)
        xXi = torch.einsum("eai,gaj->egij", Xe, dNdxi)        # (nEl,nG,3,3)
        detJ = torch.linalg.det(xXi)
        if (detJ <= 0).any():
            # flip inverted linear corners (rebuild would be needed); guard only
            raise RuntimeError("TET10: non-positive Jacobian at some Gauss point.")
        invJ = torch.linalg.inv(xXi)
        # Nx[e,g,a,i] = sum_j dNdxi[g,a,j] invJ[e,g,j,i]
        self.Nx = torch.einsum("gaj,egji->egai", dNdxi, invJ)  # (nEl,nG,10,3)
        self.wdet = wq[None, :] * detJ                         # (nEl,nG)

        # free / fixed dofs
        fixed = torch.zeros(self.ndof, dtype=torch.bool)
        for cdof in range(3):
            fixed[3 * torch.tensor(base_nodes) + cdof] = True
        self.free_mask = ~fixed
        self.free_idx = torch.nonzero(self.free_mask).flatten()
        g2f = torch.full((self.ndof,), -1, dtype=torch.long)
        g2f[self.free_idx] = torch.arange(self.free_idx.numel())
        self.g2f = g2f
        self.elem_dofs = (3 * self.tets[:, :, None] + torch.arange(3)).reshape(-1, 30)

        # dead pressure load (reference config): unit-pressure nodal force vector
        coeff = tri6_surface_coeff()
        ndir = torch.tensor(p_unit_dir)
        Fext = torch.zeros(self.nNo, 3)
        X = self.nodes
        for fa in faces6:
            x = X[torch.tensor(fa)]
            area2 = torch.linalg.norm(torch.cross(x[1] - x[0], x[2] - x[0], dim=0))
            Fext.index_add_(0, torch.tensor(fa), area2 * coeff[:, None] * ndir[None, :])
        self.Fext_unit = Fext                                 # multiply by pressure


# ---------------------------------------------------------------------------
# Internal force (vectorized, Gauss loop) and per-element tangent
# ---------------------------------------------------------------------------


def internal_force(mesh, u, c):
    """Global internal nodal force (nNo,3). c = per-element C10 (nEl,)."""
    Ue = u[mesh.tets]                                         # (nEl,10,3)
    R = torch.zeros(mesh.nNo, 3, dtype=u.dtype)
    I = torch.eye(3, dtype=u.dtype)
    for g in range(mesh.nG):
        Nxg = mesh.Nx[:, g]                                   # (nEl,10,3)
        F = I + torch.einsum("eai,eaj->eij", Ue, Nxg)
        S = pk2_stress(F, c, KFAC * c)
        P = F @ S
        fe = mesh.wdet[:, g, None, None] * torch.einsum("eaj,eij->eai", Nxg, P)
        R = R.index_add(0, mesh.tets.reshape(-1), fe.reshape(-1, 3))
    return R


def residual_free(mesh, u, c, pressure):
    R = internal_force(mesh, u, c) - pressure * mesh.Fext_unit
    return R.reshape(-1)[mesh.free_mask]


def _elem_force(ue, Nxe, wde, c):
    """Internal force of ONE TET10 element. ue(10,3), Nxe(nG,10,3), wde(nG)."""
    I = torch.eye(3, dtype=ue.dtype)
    f = torch.zeros(10, 3, dtype=ue.dtype)
    for g in range(Nxe.shape[0]):
        Nxg = Nxe[g]
        F = I + torch.einsum("ai,aj->ij", ue, Nxg)
        S = pk2_stress(F, c, KFAC * c)
        P = F @ S
        f = f + wde[g] * torch.einsum("aj,ij->ai", Nxg, P)
    return f


def assemble_sparse_tangent(mesh, u, c):
    """Sparse free-dof tangent K = d R_free / d u_free (internal part; the dead
    pressure load has zero tangent)."""
    Ue = u[mesh.tets]
    tet_jac = jacrev(lambda ue, nx, wd, ce: _elem_force(ue, nx, wd, ce), argnums=0)
    Ke = vmap(tet_jac, in_dims=(0, 0, 0, 0))(Ue, mesh.Nx, mesh.wdet, c)
    Ke = Ke.reshape(-1, 30, 30)
    gd = mesh.elem_dofs
    rows = gd[:, :, None].expand(-1, 30, 30).reshape(-1)
    cols = gd[:, None, :].expand(-1, 30, 30).reshape(-1)
    fr, fc = mesh.g2f[rows], mesh.g2f[cols]
    keep = (fr >= 0) & (fc >= 0)
    n = mesh.free_idx.numel()
    return torch.sparse_coo_tensor(torch.stack([fr[keep], fc[keep]]),
                                   Ke.reshape(-1)[keep], (n, n)).coalesce()


def solve_forward(mesh, c, pressure, n_steps=6, tol=1e-8, verbose=True):
    free = mesh.free_idx
    u = torch.zeros(mesh.nNo, 3)
    for ls in range(1, n_steps + 1):
        p = pressure * ls / n_steps
        for it in range(40):
            Rf = residual_free(mesh, u, c, p)
            rn = torch.linalg.norm(Rf).item()
            if rn < tol:
                break
            K = assemble_sparse_tangent(mesh, u, c)
            du = _sparse_lu_solve(K, -Rf)
            with torch.no_grad():
                u.reshape(-1)[free] += du
        if verbose:
            print(f"    load {ls}/{n_steps} p={p:.3f}  |R|={rn:.2e}  max|u|={u.abs().max():.4f}")
    return u.detach()


# ---------------------------------------------------------------------------
# Inverse: recover per-element C10 from displacement only, NO regularizer
# ---------------------------------------------------------------------------


def invert(mesh, u_data, pressure, c0, iters=400, lr=0.05, verbose=True, tag=""):
    with torch.no_grad():
        scale = residual_free(mesh, u_data, torch.full((mesh.nEl,), 1e-3 * c0),
                              pressure).norm().item()
    theta = torch.full((mesh.nEl,), np.log(c0), requires_grad=True)
    opt = torch.optim.Adam([theta], lr=lr)
    for it in range(iters):
        opt.zero_grad()
        c = torch.exp(theta)
        loss = residual_free(mesh, u_data, c, pressure).pow(2).sum() / scale ** 2
        loss.backward()
        opt.step()
        if verbose and (it % 100 == 0 or it == iters - 1):
            print(f"    {tag} it {it:4d}  data={loss.item():.3e}")
    return torch.exp(theta).detach()


def report(name, c_rec, c_true):
    r = c_rec / c_true
    corr = np.corrcoef(c_rec.numpy(), c_true.numpy())[0, 1]
    relerr = (c_rec - c_true).norm() / c_true.norm()
    print(f"  [{name}]  corr(c_rec,c_true)={corr:.4f}  rel-err={relerr:.3f}  "
          f"recovered ratio min/median/max = {r.min():.2f}/{r.median():.2f}/{r.max():.2f}")
    return corr, float(relerr)


if __name__ == "__main__":
    np.set_printoptions(precision=3, suppress=True)

    # ---- 1. mesh ----------------------------------------------------------
    n, L = 8, 8.0
    nodes4, tets4, endo3, base4 = box_mesh(n=n, L=L)
    nodes10, tets10, nV, lookup = build_tet10(nodes4, tets4)
    base10 = np.nonzero(nodes10[:, 2] < 1e-9)[0]              # z=0 face fixed
    faces6 = tri6_faces(endo3, nV, lookup)                    # x=0 pressure faces
    mesh = Tet10Mesh(nodes10, tets10, base10, faces6, p_unit_dir=(1.0, 0, 0))
    print(f"TET10 box: {mesh.nNo} nodes ({nV} vertices + {mesh.nNo-nV} edge mids), "
          f"{mesh.nEl} elements, free dofs {mesh.free_idx.numel()}")
    print(f"per-element unknowns {mesh.nEl}  vs  equations {mesh.free_idx.numel()}  "
          f"-> {mesh.nEl/mesh.free_idx.numel():.2f}x  "
          f"({'OVER-determined' if mesh.nEl < mesh.free_idx.numel() else 'under'})")

    # ---- 2. planted heterogeneous C10 (stiff spherical inclusion) ---------
    c0, ratio, rad = 1.0, 5.0, 2.0
    ec = nodes10[tets10].mean(1)                              # element centroids
    d = np.linalg.norm(ec - np.array([L / 2, L / 2, L / 2]), axis=1)
    c_true = torch.tensor(np.where(d < rad, c0 * ratio, c0))
    print(f"planted: bulk C10={c0}, inclusion C10={c0*ratio} (r={rad}), "
          f"{int((c_true>c0).sum())} stiff elements")

    # ---- 3. forward solve -> synthetic displacements ----------------------
    print("forward solve (planted field):")
    u_data = solve_forward(mesh, c_true, pressure=0.6)

    # ---- 4. inverse from FULL (quadratic-resolution) displacement ---------
    print("inverse from full TET10 displacement (no regularizer):")
    c_full = invert(mesh, u_data, 0.6, c0, tag="full")
    report("TET10 full-res", c_full, c_true)

    # ---- 5. degrade imaging: mid-edge dofs = linear interp of corners -----
    # emulates displacement measured only at vertex spacing (no sub-element info)
    u_lin = u_data.clone()
    # rebuild edge endpoints to average corner displacements into midedge nodes
    ev = np.stack([np.sort(tets4[:, list(e)], axis=1) for e in EDGES], axis=1).reshape(-1, 2)
    uniq = np.unique(ev, axis=0)
    mid_ids = np.array([lookup[int(a) * (nV + 1) + int(b)] for a, b in uniq])
    u_lin[mid_ids] = 0.5 * (u_data[uniq[:, 0]] + u_data[uniq[:, 1]])
    print(f"degraded imaging: mid-edge u set to corner average "
          f"(changed {np.abs((u_lin-u_data).numpy()).max():.2e} max)")
    print("inverse from vertex-resolution displacement (no regularizer):")
    c_lin = invert(mesh, u_lin, 0.6, c0, tag="lin")
    report("TET10 vertex-res", c_lin, c_true)

    np.savez("tet10_identifiability.npz", c_true=c_true.numpy(),
             c_full=c_full.numpy(), c_lin=c_lin.numpy())
    print("saved tet10_identifiability.npz")
