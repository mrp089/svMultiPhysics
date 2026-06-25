"""
inverse_all_timesteps.py
========================

Fit the per-element TET4 neo-Hookean stiffness C10 to the WHOLE 200-step
passive-inflation history of the real infarct LV, with NO regularizer -- testing
whether many (correlated) load levels add enough independent information to make
the per-element field identifiable on their own.

Loss = sum over load steps t of || R_free(U_t, c; p_t) ||^2 / scale_t^2, optimized
by minibatch Adam over the 200 steps (so all 200 are used without paying a 200x
per-step cost). nu=0.48 given -> Kpen = 49.33*C10 per element.

Run with torch + numpy.  Compares against the healthy C10 and the known infarct
location from the earlier diagnostic.
"""

import os
import numpy as np
import torch

import neohookean_equilibrium_torch as M

LAM_S = float(os.environ.get("LAM_S", 0.0))     # default: NO regularization
BATCH = int(os.environ.get("BATCH", 8))
STEPS = int(os.environ.get("STEPS", 2000))


def main():
    torch.set_printoptions(precision=4, sci_mode=True)
    g = np.load("lv_passive_inflation_infarct.npz")
    mesh = M.TorchMesh(g["nodes"], g["tets"], g["endo_faces"], g["base_nodes"])
    edges = torch.as_tensor(g["edges"], dtype=torch.long)
    C10h = float(g["C10_true"]); kfac = 4 * 1.48 / (3 * 0.04)

    a = np.load("lv_all200_infarct.npz")
    U = torch.as_tensor(a["displacements"])             # (200,nNo,3) float32
    P = torch.as_tensor(a["pressures"], dtype=torch.float64)
    pmax = P.max().item()
    use = torch.nonzero(P > 0.02 * pmax).flatten()      # skip near-zero-load steps
    print(f"mesh {mesh.nNo} nodes, {mesh.tets.shape[0]} elems; using {use.numel()}/200 "
          f"load steps (p>{0.02*pmax:.0f}); LAM_S={LAM_S}")

    # per-step normalization scale_t = ||pressure load||_free
    with torch.no_grad():
        tiny = torch.full((mesh.tets.shape[0],), 1e-3 * C10h)
        scale = {int(t): M.equilibrium_residual(
            mesh, U[t].double(), tiny, kfac * tiny, P[t]).reshape(-1)[mesh.free_mask
            ].norm().item() for t in use.tolist()}

    theta = torch.full((mesh.tets.shape[0],), np.log(C10h), requires_grad=True)
    e0, e1 = edges[:, 0], edges[:, 1]
    opt = torch.optim.Adam([theta], lr=0.05)
    gen = torch.Generator().manual_seed(0)
    for it in range(STEPS):
        opt.zero_grad()
        batch = use[torch.randint(use.numel(), (BATCH,), generator=gen)]
        c = torch.exp(theta)
        data = 0.0
        for t in batch.tolist():
            R = M.equilibrium_residual(mesh, U[t].double(), c, kfac * c, P[t])
            data = data + R.reshape(-1)[mesh.free_mask].pow(2).sum() / scale[t] ** 2
        loss = data / BATCH
        if LAM_S:
            loss = loss + LAM_S * (theta[e0] - theta[e1]).pow(2).mean()
        loss.backward(); opt.step()
        if it % 200 == 0 or it == STEPS - 1:
            with torch.no_grad():
                cc = torch.exp(theta)
                print(f"  it {it:4d}  data={float(data)/BATCH:.3e}  "
                      f"C10 p5/p50/p95/max = [{torch.quantile(cc,0.05):.3e} "
                      f"{cc.median():.3e} {torch.quantile(cc,0.95):.3e} {cc.max():.3e}]")
    c = torch.exp(theta).detach()
    np.save("c10_all200.npy", c.numpy())

    # ---- report: bulk + infarct localization --------------------------------
    bulk = c.median().item()
    ec = g["nodes"][g["tets"]].mean(1)
    ratio = (c / bulk).numpy()
    stiff = ratio > 1.3
    print(f"\nbulk (median) C10 = {bulk:.4e} ({bulk/C10h:.2f}x healthy)")
    print(f"stiff (>1.3x bulk): {100*stiff.mean():.1f}% of elements")
    if stiff.any():
        ctr = np.average(ec[stiff], axis=0, weights=ratio[stiff] - 1)
        vol = mesh.vol.numpy()[stiff].sum()
        print(f"  infarct centroid = [{ctr[0]:.2f}, {ctr[1]:.2f}, {ctr[2]:.2f}] cm "
              f"(diagnostic ~[-0.85, 3.76, -11.70])")
        print(f"  stiff volume = {vol:.2f} cm^3 (sphere-equiv r={(3*vol/4/np.pi)**(1/3):.2f} cm)")
        print(f"  peak ratio = {ratio.max():.1f}x")
    print("saved c10_all200.npy")


if __name__ == "__main__":
    main()
