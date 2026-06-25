"""
load_all200.py
==============

Extract ALL 200 timesteps of the PEXA12 infarct passive-inflation simulation
(displacement field + applied endocardial pressure per step) into one compact
.npz, so the per-element TET4 inverse can be fit against the whole load history.

Reuses the mesh / surfaces / element-adjacency already stored in
lv_passive_inflation_infarct.npz. Run with pyvista.
"""

import numpy as np
import pyvista as pv

ROOT = "/Users/pfaller/repos/PEXA12_LV_Simulations"
SUB = "passive_infl_NH_model_infarct_1.5radius"
PREFIX = "result_infarct_radius1.5_"


def main():
    p_all = np.loadtxt(f"{ROOT}/pressure_interpolated.dat", skiprows=1)[:, 1]
    disps, press = [], []
    for ts in range(1, 201):
        r = pv.read(f"{ROOT}/{SUB}/{PREFIX}{ts:03d}.vtu")
        disps.append(np.asarray(r.point_data["Displacement"], np.float32))
        press.append(float(p_all[min(ts, len(p_all) - 1)]))
        if ts % 40 == 0:
            print(f"  read {ts}/200  (p={press[-1]:.0f}, max|u|={np.abs(disps[-1]).max():.3f})")
    np.savez_compressed("lv_all200_infarct.npz",
                        displacements=np.stack(disps),       # (200,nNo,3) float32
                        pressures=np.array(press))
    print(f"wrote lv_all200_infarct.npz  ({np.stack(disps).shape}, "
          f"pmax={max(press):.0f} dyne/cm^2)")


if __name__ == "__main__":
    main()
