"""
write_estimate_vtu.py
=====================

Attach the estimated per-element neo-Hookean material to the reference volume
mesh and write a .vtu (cell data), for visualization in ParaView.

Run with a Python that has pyvista:
    python3 write_estimate_vtu.py  c10_estimate_infarct.npy  estimated_material_infarct.vtu

Cell arrays written:
  C10_estimated      per-element C10 (dyne/cm^2)
  E_estimated        Young's modulus = 4 C10 (1+nu)   (nu=0.48 given)
  C10_ratio_healthy  C10 / healthy bulk value (1.6892e5) -> infarct contrast
"""

import sys
import numpy as np
import pyvista as pv

ROOT = "/Users/pfaller/repos/PEXA12_LV_Simulations"
NU = 0.48
C10_HEALTHY = 0.25 * 1.0e6 / (1.0 + NU)          # 1.6892e5 dyne/cm^2


def main():
    npy = sys.argv[1] if len(sys.argv) > 1 else "c10_estimate_infarct.npy"
    out = sys.argv[2] if len(sys.argv) > 2 else "estimated_material_infarct.vtu"

    c10 = np.load(npy)
    vm = pv.read(f"{ROOT}/volume_mesh_5000.mesh_cm.vtu")
    assert vm.n_cells == c10.size, f"cell count {vm.n_cells} != field {c10.size}"

    vm.cell_data["C10_estimated"] = c10
    vm.cell_data["E_estimated"] = 4.0 * c10 * (1.0 + NU)
    vm.cell_data["C10_ratio_healthy"] = c10 / C10_HEALTHY
    vm.save(out)

    r = c10 / C10_HEALTHY
    print(f"wrote {out}  ({vm.n_cells} cells)")
    print(f"  C10  median={np.median(c10):.3e}  p95={np.quantile(c10,0.95):.3e}  "
          f"max={c10.max():.3e} dyne/cm^2")
    print(f"  ratio vs healthy: median={np.median(r):.2f}  p95={np.quantile(r,0.95):.2f}  "
          f"frac>1.3x={100*(r>1.3).mean():.1f}%")


if __name__ == "__main__":
    main()
