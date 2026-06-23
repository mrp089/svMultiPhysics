"""
inverse_passive_inflation.py
============================

Passive-inflation INVERSE problem on the real PEXA12 left-ventricle mesh,
identifying a PER-ELEMENT neo-Hookean stiffness field C10(e).

Given (fixed) inputs
--------------------
  * Poisson ratio nu = 0.48 (given). The volumetric penalty is DERIVED from it:
    with C10 = 0.25 E/(1+nu) and Kpen = E/(3(1-2nu)) (svMultiPhysics defaults),
    eliminating E gives, per element,
        Kpen(e) = C10(e) * 4(1+nu) / (3(1-2nu))      (= 49.33 * C10 at nu=0.48).
    So a single per-element unknown C10(e) sets both the isochoric and the
    volumetric stiffness.

Identification
--------------
The NH data is quasi-static (no viscosity, negligible inertia), so each saved
state is a static equilibrium R(U, C10) = F_internal + F_pressure = 0. The whole
internal force is linear in the per-element field C10(e), so we fit it by
minimizing, over several load levels t, the regularized equilibrium residual

    J(c) = sum_t || R_free(U_t, c; p_t) ||^2 / s_t^2
           + lam_s * smoothness(log c)        (graph-Laplacian on element faces)
           + lam_p * || log c - log c0 ||^2   (Tikhonov toward the global fit)

Per-element identification from displacements alone is ill-posed; the smoothness
+ prior terms pick the smoothest minimal-deviation field that explains the data.
With a homogeneous material the field comes back ~uniform; with the infarct it
shows a localized stiff region.

Usage (Python with torch + numpy):
    python inverse_passive_inflation.py [npz] [out_npy]
    # default: lv_passive_inflation_infarct.npz  ->  c10_estimate_infarct.npy
"""

import sys
import numpy as np
import torch

import neohookean_equilibrium_torch as M


def kpen_factor(nu):
    return 4.0 * (1.0 + nu) / (3.0 * (1.0 - 2.0 * nu))


def load(npz, dtype=torch.float64):
    d = np.load(npz)
    mesh = M.TorchMesh(d["nodes"], d["tets"], d["endo_faces"], d["base_nodes"],
                       dtype=dtype)
    U = torch.as_tensor(d["displacements"], dtype=dtype)         # (nT,nNo,3)
    P = torch.as_tensor(d["pressures"], dtype=dtype)             # (nT,)
    edges = torch.as_tensor(d["edges"], dtype=torch.long)        # (nE,2)
    return mesh, U, P, edges, float(d["nu"]), float(d["C10_true"])


def free_resid(mesh, U, c10, kfac, p):
    """Equilibrium residual on the free dofs, with Kpen derived per element."""
    kpen = kfac * c10
    R = M.equilibrium_residual(mesh, U, c10, kpen, p).reshape(-1)
    return R[mesh.free_mask]


def fit_global(mesh, U, P, kfac, c0, iters=250, lr=0.05):
    """Robust scalar baseline C10 (well-posed), as init/prior for the field."""
    logc = torch.tensor(np.log(c0), dtype=torch.float64, requires_grad=True)
    opt = torch.optim.Adam([logc], lr=lr)
    for it in range(iters):
        opt.zero_grad()
        c10 = torch.exp(logc)
        loss = sum(free_resid(mesh, U[k], c10, kfac, P[k]).pow(2).mean()
                   for k in range(U.shape[0]))
        loss.backward()
        opt.step()
    return float(torch.exp(logc))


def fit_per_element(mesh, U, P, edges, kfac, c_prior, iters=500, lr=0.02,
                    lam_s=5e-2, lam_p=1e-3, verbose=True):
    """Per-element C10 field by regularized residual minimization.

    c_prior anchors the bulk (use the known healthy stiffness, NOT the global fit
    which is biased low by the infarct). lam_s (smoothness on log C10, graph
    Laplacian over element faces) regularizes the null space; lam_p (Tikhonov to
    c_prior) pins data-insensitive (low-strain) elements to the bulk value.
    Both regularizers are per-item means so the weights are scale-intuitive.
    """
    import os
    lam_s = float(os.environ.get("LAM_S", lam_s))
    lam_p = float(os.environ.get("LAM_P", lam_p))
    iters = int(os.environ.get("ITERS", iters))
    nT = U.shape[0]
    # per-load normalization s_t = || pressure load ||_free  (residual at c->0)
    with torch.no_grad():
        tiny = torch.full((mesh.tets.shape[0],), 1e-3 * c_prior, dtype=torch.float64)
        scales = [free_resid(mesh, U[k], tiny, kfac, P[k]).norm().item()
                  for k in range(nT)]

    logc0 = np.log(c_prior)
    theta = torch.full((mesh.tets.shape[0],), logc0, dtype=torch.float64,
                       requires_grad=True)
    e0, e1 = edges[:, 0], edges[:, 1]
    opt = torch.optim.Adam([theta], lr=lr)
    print(f"  (lam_s={lam_s}, lam_p={lam_p}, iters={iters}, prior C10={c_prior:.3e})")
    for it in range(iters):
        opt.zero_grad()
        c10 = torch.exp(theta)
        data = sum(free_resid(mesh, U[k], c10, kfac, P[k]).pow(2).sum()
                   / scales[k] ** 2 for k in range(nT))      # sum of rel-resid^2
        smooth = (theta[e0] - theta[e1]).pow(2).mean()
        prior = (theta - logc0).pow(2).mean()
        loss = data + lam_s * smooth + lam_p * prior
        loss.backward()
        opt.step()
        if verbose and (it % 100 == 0 or it == iters - 1):
            with torch.no_grad():
                c = torch.exp(theta)
                print(f"  it {it:4d}  data={float(data):.4e}  smooth={float(smooth):.3e}"
                      f"  C10[p5/med/p95/max]=[{torch.quantile(c,0.05):.3e} "
                      f"{c.median():.3e} {torch.quantile(c,0.95):.3e} {c.max():.3e}]")
    return torch.exp(theta).detach().numpy()


def main():
    npz = sys.argv[1] if len(sys.argv) > 1 else "lv_passive_inflation_infarct.npz"
    out = sys.argv[2] if len(sys.argv) > 2 else "c10_estimate_infarct.npy"
    torch.set_printoptions(precision=4, sci_mode=True)

    mesh, U, P, edges, nu, C10_true = load(npz)
    kfac = kpen_factor(nu)
    print(f"dataset: {npz}")
    print(f"  {mesh.nNo} nodes, {mesh.tets.shape[0]} tets, loads={list(P.numpy().astype(int))} dyne/cm^2")
    print(f"  nu={nu} (given) -> Kpen = {kfac:.3f} * C10;  healthy C10_true={C10_true:.4e}\n")

    print("Global scalar baseline:")
    c0 = fit_global(mesh, U, P, kfac, c0=C10_true)
    print(f"  global C10 = {c0:.4e}  ({(c0/C10_true-1)*100:+.2f}% vs healthy)")
    print("  (biased low when an infarct is present; bulk prior uses healthy value)\n")

    print("Per-element field (regularized residual minimization):")
    # anchor the prior at the known healthy bulk stiffness, not the biased global fit
    c10 = fit_per_element(mesh, U, P, edges, kfac, c_prior=C10_true)
    np.save(out, c10)

    # ---- report: locate the stiffened (infarct) region ----------------------
    elem_c = mesh.nodes.numpy()[mesh.tets.numpy()].mean(axis=1)   # element centroids
    ratio = c10 / C10_true                       # vs healthy bulk
    stiff = ratio > 1.3                          # >30% stiffer than healthy bulk
    print(f"\nsaved per-element C10 -> {out}")
    print(f"  field C10  min/median/max = {c10.min():.3e} / {np.median(c10):.3e} / {c10.max():.3e}")
    print(f"  elements >1.3x baseline: {stiff.sum()} ({100*stiff.mean():.2f}% of mesh)")
    if stiff.any():
        c = elem_c[stiff]
        w = (ratio[stiff] - 1.0)
        center = np.average(c, axis=0, weights=w)
        # effective radius from stiff-region volume (sphere-equivalent)
        vol = mesh.vol.numpy()[stiff].sum()
        reff = (3.0 * vol / (4.0 * np.pi)) ** (1.0 / 3.0)
        print(f"  stiff-region centroid (cm) = [{center[0]:.2f} {center[1]:.2f} {center[2]:.2f}]")
        print(f"  stiff-region volume = {vol:.3f} cm^3  (sphere-equiv radius {reff:.2f} cm)")


if __name__ == "__main__":
    main()
