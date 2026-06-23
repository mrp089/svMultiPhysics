"""
combined_vtu.py
===============

Pack every field from region_sweep.npz into ONE .vtu as separate cell arrays, so
the impact of the regularization strength (region size K) can be visualized by
scrubbing through the arrays in ParaView.

Run with a Python that has pyvista:
    python3 combined_vtu.py
"""

import numpy as np
import pyvista as pv

ROOT = "/Users/pfaller/repos/PEXA12_LV_Simulations"
NU = 0.48
OUT = "regularization_sweep.vtu"


def main():
    s = np.load("region_sweep.npz")
    hs, Ks, fields, C10h = s["h"], s["K"], s["fields"], float(s["C10_true"])
    vm = pv.read(f"{ROOT}/volume_mesh_5000.mesh_cm.vtu")

    # order coarse->fine so array list reads high-regularization -> low
    order = np.argsort(Ks)
    for i in order:
        h, K, c = hs[i], int(Ks[i]), fields[i]
        tag = f"K{K:05d}_h{h:.2f}"
        vm.cell_data[f"C10_{tag}"] = c
        vm.cell_data[f"ratio_{tag}"] = c / C10h        # contrast vs healthy bulk
    vm.save(OUT)

    print(f"wrote {OUT}  ({vm.n_cells} cells, {len(hs)} regularization levels)")
    print(f"{'K':>7} {'h(cm)':>6} {'median':>10} {'peak/healthy':>13} {'>1.3x %':>8}")
    for i in order:
        r = fields[i] / C10h
        print(f"{int(Ks[i]):>7} {hs[i]:>6.2f} {np.median(fields[i]):>10.3e} "
              f"{r.max():>13.1f} {100*(r>1.3).mean():>8.1f}")
    print("\nIn ParaView: color by ratio_* arrays; small K = high regularization "
          "(smooth, low contrast), large K = low regularization (sharp, noisier).")


if __name__ == "__main__":
    main()
