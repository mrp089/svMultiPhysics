"""Well-posed, regularizer-free per-region C10: minimize the (nodal) balance
residual with C10 piecewise-constant on a voxel partition of the elements.
Unknowns K = #occupied voxels << #equations -> unique minimizer, no smoothness."""
import sys, numpy as np, torch
import neohookean_equilibrium_torch as M

H = float(sys.argv[1]) if len(sys.argv)>1 else 1.0   # voxel edge (cm)
d=np.load("lv_passive_inflation_infarct.npz")
mesh=M.TorchMesh(d["nodes"],d["tets"],d["endo_faces"],d["base_nodes"])
U=torch.as_tensor(d["displacements"]); P=torch.as_tensor(d["pressures"])
C10h=float(d["C10_true"]); kfac=4*(1+0.48)/(3*(1-2*0.48))

# voxel-partition elements by centroid
ctr=d["nodes"][d["tets"]].mean(1)
vox=np.floor((ctr-ctr.min(0))/H).astype(np.int64)
key=vox[:,0]+1000*(vox[:,1]+1000*vox[:,2])
_,region=np.unique(key,return_inverse=True)
K=int(region.max()+1); region=torch.as_tensor(region)
print(f"h={H}cm  regions K={K}  vs nFree={int(mesh.free_mask.sum())}  "
      f"({'OVER-determined, unique' if K<mesh.free_mask.sum() else 'still underdetermined'})")

theta=torch.full((K,),np.log(C10h),dtype=torch.float64,requires_grad=True)
with torch.no_grad():
    sc=[M.equilibrium_residual(mesh,U[k],1e-3*C10h,kfac*1e-3*C10h,P[k]).reshape(-1)[mesh.free_mask].norm().item() for k in range(len(P))]
opt=torch.optim.Adam([theta],lr=0.05)
for it in range(400):
    opt.zero_grad(); c=torch.exp(theta)[region]
    loss=sum(M.equilibrium_residual(mesh,U[k],c,kfac*c,P[k]).reshape(-1)[mesh.free_mask].pow(2).sum()/sc[k]**2 for k in range(len(P)))
    loss.backward(); opt.step()
c=torch.exp(theta).detach()
celem=c[region].numpy()
np.save(f"c10_region_h{H}.npy", celem)
r=celem/C10h
hi=r>1.3
cc=np.average(ctr[hi],axis=0,weights=(r[hi]-1)) if hi.any() else [0,0,0]
vol=mesh.vol.numpy()[hi].sum()
print(f"  data(rel-resid^2 sum)={float(loss):.3e}  region C10 med/max={c.median():.3e}/{c.max():.3e}")
print(f"  >1.3x healthy: {100*hi.mean():.1f}%  centroid=[{cc[0]:.2f},{cc[1]:.2f},{cc[2]:.2f}] "
      f"r_eff={(3*vol/4/np.pi)**(1/3):.2f}cm  peak ratio={r.max():.1f}x")
