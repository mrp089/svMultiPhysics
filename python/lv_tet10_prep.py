"""
lv_tet10_prep.py
================

Prepare a COARSE TET10 problem from the REAL PEXA12 infarct data so that
quadratic-element per-element identification can be run on measured displacements.

Why coarse: the svMultiPhysics ground truth is a LINEAR TET4 solution sampled at
the fine mesh vertices. TET10 on the *same* mesh gains nothing (mid-edge u = exact
corner average -> the degraded case). The fix is TET10 elements COARSER than the
data spacing (~0.14 cm here), so each mid-edge node samples an independent
displacement from the fine field. We build a ~0.4 cm coarse TET10 mesh, interpolate
the fine displacement onto it, and classify boundary faces for the BCs.

Output: lv_tet10_infarct.npz  (nodes, TET10 connectivity, measured displacement,
endocardial TRI6 faces, base nodes, pressure, healthy C10).

Run with pyvista + tetgen + scipy.
"""

import numpy as np
import pyvista as pv
import tetgen
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import cKDTree

ROOT = "/Users/pfaller/repos/PEXA12_LV_Simulations"
PEAK = 107                                     # peak-pressure timestep
DECIMATE_TRIS = 1200                           # coarse surface target
MAXVOL = 0.08                                  # tetgen max element volume (cm^3)

# svMultiPhysics TET10 local edges (nodes 4..9)
EDGES = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
DN_TET4 = np.array([[1, 0, 0, -1], [0, 1, 0, -1], [0, 0, 1, -1.0]])


def build_tet10(nodes4, tets4):
    ev = np.stack([np.sort(tets4[:, list(e)], axis=1) for e in EDGES], axis=1)
    uniq, inv = np.unique(ev.reshape(-1, 2), axis=0, return_inverse=True)
    mids = 0.5 * (nodes4[uniq[:, 0]] + nodes4[uniq[:, 1]])
    nV = nodes4.shape[0]
    tets10 = np.hstack([tets4, (nV + inv).reshape(-1, 6)]).astype(np.int64)
    nodes10 = np.vstack([nodes4, mids])
    key = uniq[:, 0].astype(np.int64) * (nV + 1) + uniq[:, 1]
    lookup = {int(k): i + nV for i, k in enumerate(key)}
    return nodes10, tets10, nV, lookup


def tet10_N(xi):
    """TET10 shape functions (svMultiPhysics order). xi (...,3) -> (...,10)."""
    r, s, t = xi[..., 0], xi[..., 1], xi[..., 2]
    u = 1.0 - r - s - t
    return np.stack([r * (2 * r - 1), s * (2 * s - 1), t * (2 * t - 1),
                     u * (2 * u - 1), 4 * r * s, 4 * s * t, 4 * r * t,
                     4 * r * u, 4 * s * u, 4 * t * u], axis=-1)


def tet10_quad():
    """15-point degree-5 Keast rule (svMultiPhysics)."""
    w = np.array([0.0302836780970890] + [0.0060267857142860] * 4
                 + [0.0116452490860290] * 4 + [0.0109491415613860] * 6)
    s0 = 0.25; a, b = 1 / 3, 0.0; c, d = 0.0909090909090910, 0.7272727272727270
    e, f = 0.0665501535736640, 0.4334498464263360
    xi = np.array([[s0, s0, s0], [b, a, a], [a, b, a], [a, a, b], [a, a, a],
                   [d, c, c], [c, d, c], [c, c, d], [c, c, c],
                   [e, e, f], [e, f, e], [e, f, f], [f, f, e], [f, e, f], [f, e, e]])
    return w, xi


def l2_project(nodes10, tets10, fine):
    """Galerkin L2 projection of the fine displacement onto the coarse TET10
    space: solve (int N^T N) u_hat = int N^T u_fine, i.e. M u_hat = b, where the
    integrals use the 15-pt rule and u_fine is sampled at the coarse Gauss points
    (NOT just the 10 nodes -> no aliasing). Returns u_hat (nNo,3)."""
    w, xiq = tet10_quad()
    Ng = tet10_N(xiq)                                    # (15,10)
    # dN/dxi by central differences (only needed for detJ)
    h = 1e-6
    dNg = np.stack([(tet10_N(xiq + np.eye(3)[j] * h)
                     - tet10_N(xiq - np.eye(3)[j] * h)) / (2 * h)
                    for j in range(3)], axis=-1)         # (15,10,3)
    Xe = nodes10[tets10]                                 # (nEl,10,3)
    detJ = np.linalg.det(np.einsum("eai,gaj->egij", Xe, dNg))   # (nEl,15)
    We = w[None, :] * detJ                               # (nEl,15) quad weights
    Xg = np.einsum("ga,eac->egc", Ng, Xe)               # (nEl,15,3) phys gauss pts

    # sample u_fine at all coarse Gauss points; kernel-fill any outside fine mesh
    pts = Xg.reshape(-1, 3)
    smp = pv.PolyData(pts).sample(fine, tolerance=1e-3)
    ug = np.asarray(smp["Displacement"], float)
    bad = ~np.asarray(smp["vtkValidPointMask"]).astype(bool)
    if bad.any():
        ker = pv.PolyData(pts[bad]).interpolate(fine, radius=0.35, sharpness=4,
                                                strategy="closest_point")
        ug[bad] = np.asarray(ker["Displacement"], float)
    ug = ug.reshape(Xe.shape[0], 15, 3)

    # element mass matrices and rhs
    Me = np.einsum("eg,ga,gb->eab", We, Ng, Ng)         # (nEl,10,10)
    be = np.einsum("eg,ga,egc->eac", We, Ng, ug)        # (nEl,10,3)
    nNo = nodes10.shape[0]
    rows = np.repeat(tets10, 10, axis=1).reshape(-1)
    cols = np.tile(tets10, 10).reshape(-1)
    M = sp.coo_matrix((Me.reshape(-1), (rows, cols)), shape=(nNo, nNo)).tocsc()
    b = np.zeros((nNo, 3))
    np.add.at(b, tets10.reshape(-1), be.reshape(-1, 3))
    u_hat = np.column_stack([spla.spsolve(M, b[:, c]) for c in range(3)])
    return u_hat, bad.mean()


def tri6(tri3, nV, lookup):
    a, b, c = tri3
    def mid(x, y):
        x, y = (x, y) if x < y else (y, x)
        return lookup[int(x) * (nV + 1) + int(y)]
    return [a, b, c, mid(a, b), mid(b, c), mid(c, a)]


def main():
    # --- fine reference mesh + measured displacement (infarct, peak) ---------
    fine = pv.read(f"{ROOT}/passive_infl_NH_model_infarct_1.5radius/"
                   f"result_infarct_radius1.5_{PEAK:03d}.vtu")
    press = np.loadtxt(f"{ROOT}/pressure_interpolated.dat", skiprows=1)[PEAK, 1]

    # --- coarse tet mesh from the (decimated) LV surface ---------------------
    vm = pv.read(f"{ROOT}/volume_mesh_5000.mesh_cm.vtu")
    surf = vm.extract_surface(algorithm="dataset_surface").triangulate().clean()
    s = surf.decimate(1 - DECIMATE_TRIS / surf.n_cells).clean()
    tg = tetgen.TetGen(s.points, s.faces.reshape(-1, 4)[:, 1:])
    tg.tetrahedralize(order=1, mindihedral=20, minratio=1.5, maxvolume=MAXVOL)
    nodes4 = np.asarray(tg.node, float)
    tets4 = np.asarray(tg.elem, np.int64)
    # ensure positive Jacobian
    det = np.linalg.det((nodes4[tets4].transpose(0, 2, 1) @ DN_TET4.T))
    tets4[det < 0] = tets4[det < 0][:, [0, 1, 3, 2]]
    print(f"coarse mesh: {nodes4.shape[0]} vertices, {tets4.shape[0]} tets")

    nodes10, tets10, nV, lookup = build_tet10(nodes4, tets4)
    print(f"TET10: {nodes10.shape[0]} nodes ({nV} corners + {nodes10.shape[0]-nV} mids)")

    # --- boundary faces of the coarse mesh -> classify endo/epi/base ---------
    F = np.sort(np.concatenate([tets4[:, [0, 1, 2]], tets4[:, [0, 1, 3]],
                                tets4[:, [0, 2, 3]], tets4[:, [1, 2, 3]]]), axis=1)
    order = np.lexsort(F.T); Fs = F[order]
    same = np.all(Fs[1:] == Fs[:-1], axis=1)
    interior = np.zeros(len(Fs), bool); interior[:-1] |= same; interior[1:] |= same
    bfaces = Fs[~interior]                        # (nB,3) boundary triangles
    bc = nodes4[bfaces].mean(1)
    trees = {n: cKDTree(pv.read(f"{ROOT}/mesh_surfaces/{n}.vtp").points)
             for n in ["endocardium", "epicardium", "endoepiconnection"]}
    dist = np.stack([trees[n].query(bc)[0] for n in
                     ["endocardium", "epicardium", "endoepiconnection"]], axis=1)
    label = dist.argmin(1)                         # 0=endo,1=epi,2=base
    print(f"boundary faces: {(label==0).sum()} endo, {(label==1).sum()} epi, "
          f"{(label==2).sum()} base")

    endo_faces6 = np.array([tri6(t, nV, lookup) for t in bfaces[label == 0]])
    base_tris6 = np.array([tri6(t, nV, lookup) for t in bfaces[label == 2]])
    base_nodes = np.unique(base_tris6)

    # orient endo TRI6 so cross(g1,g2) points toward cavity centroid
    cav = nodes10[np.unique(endo_faces6)].mean(0)
    x = nodes10[endo_faces6[:, :3]]
    nrm = np.cross(x[:, 0] - x[:, 2], x[:, 1] - x[:, 2])
    flip = np.einsum("ij,ij->i", nrm, cav - x.mean(1)) < 0
    endo_faces6[flip] = endo_faces6[flip][:, [0, 2, 1, 5, 4, 3]]

    # --- fit (L2-project) the fine displacement onto the coarse TET10 space ---
    # Best coarse representation of the WHOLE fine field (vs point-sampling at
    # nodes, which aliases). Solves the consistent-mass projection M u_hat = b.
    u_meas, frac_out = l2_project(nodes10, tets10, fine)
    print(f"L2 projection: {100*frac_out:.1f}% of Gauss pts outside fine mesh "
          f"(kernel-filled),  max|u|={np.abs(u_meas).max():.4f} cm")

    np.savez_compressed(
        "lv_tet10_infarct.npz", nodes=nodes10, tets=tets10,
        endo_faces=endo_faces6, base_nodes=base_nodes, u_meas=u_meas,
        pressure=press, C10_healthy=0.25e6 / 1.48)
    print(f"pressure={press:.1f} dyne/cm^2  -> wrote lv_tet10_infarct.npz")


if __name__ == "__main__":
    main()
