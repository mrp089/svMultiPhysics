"""
region_sweep.py
===============

Sweep the region size K (number of voxel partitions on which C10 is held
constant) and minimize the *balance-equation residual* with NO smoothness
regularizer at each K. The region count acts as the regularization knob:

  * small voxels  -> large K  -> almost NO regularization (high resolution,
                                  approaches the ill-posed per-element limit,
                                  more variance / noise)
  * large voxels  -> small K  -> HIGH regularization (few parameters, very
                                  smooth, biased, low contrast)

Each fit is a unique minimizer as long as K << #nodal-balance-equations, so the
sweep isolates the effect of the regularization strength alone. Saves all
per-element fields to region_sweep.npz for combined_vtu.py to pack into one VTU.

Run with a torch+numpy Python. ~minutes (run in background).
"""

import numpy as np
import torch

import neohookean_equilibrium_torch as M

# fine (~no regularization) -> coarse (high regularization)
H_LIST = [0.2, 0.35, 0.6, 1.0, 1.6, 2.6]      # voxel edge in cm
ITERS = 500

d = np.load("lv_passive_inflation_infarct.npz")
mesh = M.TorchMesh(d["nodes"], d["tets"], d["endo_faces"], d["base_nodes"])
U = torch.as_tensor(d["displacements"])
P = torch.as_tensor(d["pressures"])
C10h = float(d["C10_true"])
nu = 0.48
kfac = 4 * (1 + nu) / (3 * (1 - 2 * nu))
nFree = int(mesh.free_mask.sum())
ctr = d["nodes"][d["tets"]].mean(1)

# per-load normalization (pressure-load norm)
with torch.no_grad():
    sc = [M.equilibrium_residual(mesh, U[k], 1e-3 * C10h, kfac * 1e-3 * C10h,
          P[k]).reshape(-1)[mesh.free_mask].norm().item() for k in range(len(P))]


def fit(region, K):
    theta = torch.full((K,), np.log(C10h), dtype=torch.float64, requires_grad=True)
    opt = torch.optim.Adam([theta], lr=0.05)
    for it in range(ITERS):
        opt.zero_grad()
        c = torch.exp(theta)[region]
        loss = sum(M.equilibrium_residual(mesh, U[k], c, kfac * c, P[k])
                   .reshape(-1)[mesh.free_mask].pow(2).sum() / sc[k] ** 2
                   for k in range(len(P)))
        loss.backward()
        opt.step()
    return torch.exp(theta).detach()[region].numpy(), float(loss)


fields, Ks, hs = [], [], []
for H in H_LIST:
    vox = np.floor((ctr - ctr.min(0)) / H).astype(np.int64)
    key = vox[:, 0] + 1000 * (vox[:, 1] + 1000 * vox[:, 2])
    _, region = np.unique(key, return_inverse=True)
    K = int(region.max() + 1)
    c, data = fit(torch.as_tensor(region), K)
    r = c / C10h
    fields.append(c); Ks.append(K); hs.append(H)
    print(f"h={H:4.2f}cm  K={K:5d}  data={data:.3e}  "
          f"median={np.median(c):.3e}  peak_ratio={r.max():.1f}x  "
          f">1.3x={100*(r>1.3).mean():.1f}%")

np.savez_compressed("region_sweep.npz", h=np.array(hs), K=np.array(Ks),
                    fields=np.stack(fields), C10_true=C10h)
print("saved region_sweep.npz")
