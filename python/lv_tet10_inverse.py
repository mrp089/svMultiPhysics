"""
lv_tet10_inverse.py
===================

Recover the per-element neo-Hookean stiffness C10 of the REAL PEXA12 infarct LV
from the MEASURED displacement field alone, using a coarse TET10 (quadratic)
discretization -- no stress input, no regularizer.

Inputs (from lv_tet10_prep.py): a coarse TET10 mesh, the fine svMultiPhysics
infarct displacement interpolated onto it, the endocardial TRI6 faces, the fixed
base nodes, and the applied endocardial pressure.

Method: the measured displacement u is fixed; the only unknown is the per-element
field C10. The endocardial FOLLOWER pressure load is evaluated once at the
measured deformed configuration (so it is a known constant force; no follower
tangent needed). We then minimize the equilibrium residual

    J(c) = || R_internal_free(u, c) - F_pressure_free ||^2

over c (one value per element). Because the known pressure fixes the force scale,
this recovers ABSOLUTE C10. This is a genuine inverse on real data -- the coarse
TET10 model differs from the fine TET4 ground truth, so there is a discretization
floor (no inverse crime).

Run with torch + numpy + scipy.
"""

import numpy as np
import torch

from tet10_identifiability import Tet10Mesh, internal_force, KFAC
from neohookean_equilibrium_torch import pk2_stress

torch.set_default_dtype(torch.float64)


def internal_force_masked(mesh, u, c, good):
    """Internal nodal force with degenerate (det F<=0) elements zeroed out: their
    F is replaced by I (stress-free) so pk2_stress never sees a negative Jacobian,
    and their force contribution is masked to 0."""
    Ue = u[mesh.tets]
    I = torch.eye(3, dtype=u.dtype)
    R = torch.zeros(mesh.nNo, 3, dtype=u.dtype)
    gm = good[:, None, None]                            # bool (nEl,1,1)
    gf = good.double()[:, None, None]
    for g in range(mesh.nG):
        Nxg = mesh.Nx[:, g]
        F = I + torch.einsum("eai,eaj->eij", Ue, Nxg)
        F = torch.where(gm, F, I)                       # neutralize bad elements
        S = pk2_stress(F, c, KFAC * c)
        P = F @ S
        fe = (mesh.wdet[:, g, None, None] * gf
              * torch.einsum("eaj,eij->eai", Nxg, P))
        R = R.index_add(0, mesh.tets.reshape(-1), fe.reshape(-1, 3))
    return R


def tri6_shape_and_grad():
    """TRI6 shape M(zeta) (6,) and dM/dzeta (6,2) at a 6-point degree-4 rule.
    Node order [c0,c1,c2,m01,m12,m20]; weights sum to the param area 1/2."""
    a, b = 0.445948490915965, 0.091576213509771
    zeta = torch.tensor([[a, a], [1 - 2 * a, a], [a, 1 - 2 * a],
                         [b, b], [1 - 2 * b, b], [b, 1 - 2 * b]])
    w = torch.tensor([0.111690794839005] * 3 + [0.054975871827661] * 3)

    def M(z):
        r, s = z[0], z[1]; u = 1 - r - s
        return torch.stack([r * (2 * r - 1), s * (2 * s - 1), u * (2 * u - 1),
                            4 * r * s, 4 * s * u, 4 * u * r])
    Mg = torch.stack([M(zeta[g]) for g in range(6)])                  # (6nodes? -> 6gp,6)
    dMg = torch.stack([torch.autograd.functional.jacobian(M, zeta[g])
                       for g in range(6)])                            # (6gp,6,2)
    return w, Mg, dMg


def follower_pressure_load(mesh, faces, u, pressure):
    """Endocardial follower-pressure nodal force, evaluated at measured u.
    R_a = -p * sum_g w_g M_a (g1 x g2),  g1,g2 = deformed in-plane tangents."""
    w, Mg, dMg = tri6_shape_and_grad()
    X = mesh.nodes[faces]
    x = X + u[faces]                                                  # (nF,6,3) current
    F = torch.zeros(mesh.nNo, 3)
    for g in range(w.shape[0]):
        g1 = torch.einsum("a,fac->fc", dMg[g, :, 0], x)
        g2 = torch.einsum("a,fac->fc", dMg[g, :, 1], x)
        nA = torch.linalg.cross(g1, g2)                              # (nF,3)
        contrib = -pressure * w[g] * Mg[g][None, :, None] * nA[:, None, :]
        F = F.index_add(0, faces.reshape(-1), contrib.reshape(-1, 3))
    return F


def main():
    d = np.load("lv_tet10_infarct.npz")
    faces = torch.tensor(d["endo_faces"])
    mesh = Tet10Mesh(d["nodes"], d["tets"], d["base_nodes"], d["endo_faces"])
    u = torch.tensor(d["u_meas"])
    p = float(d["pressure"]); c0 = float(d["C10_healthy"])
    print(f"coarse TET10: {mesh.nNo} nodes, {mesh.nEl} elements, free dofs "
          f"{mesh.free_idx.numel()}  (unknowns/eqns = {mesh.nEl/mesh.free_idx.numel():.2f}x)")
    print(f"pressure={p:.0f} dyne/cm^2  healthy C10={c0:.4e}")

    Fext = follower_pressure_load(mesh, faces, u, p)                  # constant load
    Ff = Fext.reshape(-1)[mesh.free_mask]

    # Guard: coarse interpolation of the measured field can invert a few badly
    # shaped elements (det F <= 0); the neo-Hookean model is undefined there, so
    # zero out their stress contribution (their C10 stays at the prior).
    I = torch.eye(3)
    Ue = u[mesh.tets]
    detmin = torch.full((mesh.nEl,), float("inf"))
    for g in range(mesh.nG):
        F = I + torch.einsum("eai,eaj->eij", Ue, mesh.Nx[:, g])
        detmin = torch.minimum(detmin, torch.linalg.det(F))
    strain = []
    for g in range(mesh.nG):
        F = I + torch.einsum("eai,eaj->eij", Ue, mesh.Nx[:, g])
        strain.append((F - I).reshape(mesh.nEl, -1).norm(dim=1))
    strain = torch.stack(strain).max(0).values
    # mask degenerate (det F<=0) and unphysically-strained (aliasing/sliver) elems
    good = (detmin > 0.05) & (strain < 0.7)
    print(f"valid elements (det F>0.05 & strain<0.7): {int(good.sum())}/{mesh.nEl} "
          f"({100*good.double().mean():.1f}%)")
    Rint0 = internal_force_masked(mesh, u, torch.full((mesh.nEl,), c0), good
                                  ).reshape(-1)[mesh.free_mask]
    print(f"DIAG @ uniform healthy C10:  ||R_int||/||F_press|| = "
          f"{(Rint0.norm()/Ff.norm()).item():.2f}   (==1 would be equilibrium)")
    print(f"  element strain ||F-I||  p50/p95/max = "
          f"{strain.median():.2f}/{torch.quantile(strain,0.95):.2f}/{strain.max():.2f}")
    print(f"  follower load ||F_press|| = {Ff.norm():.3e}")

    def resid_free(c):
        return internal_force_masked(mesh, u, c, good).reshape(-1)[mesh.free_mask] - Ff

    # ---- inverse: per-element C10 from displacement, NO regularizer ----------
    scale = Ff.norm().item()
    theta = torch.full((mesh.nEl,), np.log(c0), requires_grad=True)
    opt = torch.optim.Adam([theta], lr=0.05)
    for it in range(500):
        opt.zero_grad()
        loss = resid_free(torch.exp(theta)).pow(2).sum() / scale ** 2
        loss.backward(); opt.step()
        if it % 100 == 0 or it == 499:
            with torch.no_grad():
                c = torch.exp(theta)
                print(f"  it {it:4d}  rel-resid^2={loss.item():.3e}  "
                      f"C10 p50/p95/max = {c.median():.3e}/{torch.quantile(c,0.95):.3e}/{c.max():.3e}")
    c = torch.exp(theta).detach()
    np.save("c10_tet10_infarct.npy", c.numpy())

    # ---- report: bulk + infarct localization --------------------------------
    bulk = c.median().item()
    ec = d["nodes"][d["tets"]].mean(1)
    ratio = (c / bulk).numpy()
    stiff = ratio > 1.3
    print(f"\nbulk (median) C10 = {bulk:.4e}  ({bulk/c0:.2f}x healthy)")
    print(f"stiff (>1.3x bulk): {100*stiff.mean():.1f}% of elements")
    if stiff.any():
        wgt = ratio[stiff] - 1
        ctr = np.average(ec[stiff], axis=0, weights=wgt)
        vol = mesh.wdet.sum(1).numpy()[stiff].sum()
        print(f"  infarct centroid = [{ctr[0]:.2f}, {ctr[1]:.2f}, {ctr[2]:.2f}] cm  "
              f"(diagnostic said ~[-0.85, 3.76, -11.70])")
        print(f"  stiff volume = {vol:.2f} cm^3 (sphere-equiv r={ (3*vol/4/np.pi)**(1/3):.2f} cm)")
        print(f"  peak ratio = {ratio.max():.1f}x healthy")
    print("saved c10_tet10_infarct.npy")


if __name__ == "__main__":
    main()
