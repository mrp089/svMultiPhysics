"""
load_lv_data.py
===============

Extract a PEXA12 LV passive-inflation dataset into a compact .npz that the
PyTorch inverse solver can consume (so the heavy VTK reading is done once and the
torch step needs no VTK dependency).

Run with a Python that has pyvista + scipy (e.g. the system interpreter):

    python3 load_lv_data.py                 # both datasets (healthy + infarct)

Data source (cgs units): /Users/pfaller/repos/PEXA12_LV_Simulations
  * volume_mesh_5000.mesh_cm.vtu          reference mesh (31436 nodes, linear tets)
  * <dataset>/<prefix>NNN.vtu             Displacement field per timestep
  * mesh_surfaces/{endocardium,epicardium,endoepiconnection}.vtp  boundary faces
  * pressure_interpolated.dat             endocardial follower pressure vs time

BCs (solver_Neohookean.xml): endocardium -> follower pressure; endoepiconnection
("top") -> fixed Dirichlet; epicardium -> traction free.

The .npz also stores element face-adjacency edges, used by the per-element
inverse for a smoothness (graph-Laplacian) regularizer.
"""

import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

ROOT = "/Users/pfaller/repos/PEXA12_LV_Simulations"

# result_NNN <-> pressure row NNN. Several load levels along the inflation
# (peak pressure is row 107) improve conditioning of the per-element inverse.
TIMESTEPS = [40, 70, 90, 107]

DATASETS = {
    "lv_passive_inflation.npz": ("passive_infl_NH_model", "result_"),
    "lv_passive_inflation_infarct.npz":
        ("passive_infl_NH_model_infarct_1.5radius", "result_infarct_radius1.5_"),
}


def element_face_adjacency(tets):
    """Element adjacency via shared triangular faces: returns edges (nE,2) of
    element-index pairs that share a face. Used for smoothness regularization."""
    nEl = tets.shape[0]
    faces = np.concatenate([tets[:, [0, 1, 2]], tets[:, [0, 1, 3]],
                            tets[:, [0, 2, 3]], tets[:, [1, 2, 3]]], axis=0)
    faces = np.sort(faces, axis=1)
    elem = np.tile(np.arange(nEl), 4)
    order = np.lexsort(faces.T)
    f, e = faces[order], elem[order]
    same = np.all(f[1:] == f[:-1], axis=1)
    return np.stack([e[:-1][same], e[1:][same]], axis=1).astype(np.int64)


def build_common():
    """Mesh, surfaces, adjacency shared by all datasets."""
    vm = pv.read(f"{ROOT}/volume_mesh_5000.mesh_cm.vtu")
    nodes = np.asarray(vm.points, dtype=np.float64)
    tets = vm.cells_dict[10].astype(np.int64)
    tree = cKDTree(nodes)

    def surf(vtp):
        s = pv.read(f"{ROOT}/mesh_surfaces/{vtp}.vtp")
        d, idx = tree.query(s.points)
        assert d.max() < 1e-10, f"{vtp}: surface points do not match volume nodes"
        tris = idx[s.faces.reshape(-1, 4)[:, 1:]].astype(np.int64)
        return tris, idx

    endo_faces, _ = surf("endocardium")
    _, base_idx = surf("endoepiconnection")
    base_nodes = np.unique(base_idx).astype(np.int64)

    # Orient endo triangles so cross(g1,g2) points toward the cavity centroid
    # (solid outward normal) -> positive pressure inflates, matching svMultiPhysics.
    cav = nodes[np.unique(endo_faces)].mean(axis=0)
    x = nodes[endo_faces]
    nrm = np.cross(x[:, 0] - x[:, 2], x[:, 1] - x[:, 2])
    flip = np.einsum("ij,ij->i", nrm, cav - x.mean(axis=1)) < 0.0
    endo_faces[flip] = endo_faces[flip][:, [0, 2, 1]]

    edges = element_face_adjacency(tets)
    print(f"common: {nodes.shape[0]} nodes, {tets.shape[0]} tets, "
          f"{endo_faces.shape[0]} endo tris, {base_nodes.size} fixed nodes, "
          f"{edges.shape[0]} element-adjacency edges")
    return nodes, tets, endo_faces, base_nodes, edges


def main():
    nodes, tets, endo_faces, base_nodes, edges = build_common()
    praw = np.loadtxt(f"{ROOT}/pressure_interpolated.dat", skiprows=1)
    p_val = praw[:, 1]
    nu = 0.48                                    # given Poisson ratio
    # ground-truth healthy neo-Hookean (E=1e6): C10 = 0.25 E/(1+nu)
    C10_true = 0.25 * 1.0e6 / (1.0 + nu)

    for out, (sub, prefix) in DATASETS.items():
        disps, pressures = [], []
        for ts in TIMESTEPS:
            r = pv.read(f"{ROOT}/{sub}/{prefix}{ts:03d}.vtu")
            u = np.asarray(r.point_data["Displacement"], dtype=np.float64)
            disps.append(u)
            pressures.append(float(p_val[ts]))
        np.savez_compressed(
            out, nodes=nodes, tets=tets, endo_faces=endo_faces,
            base_nodes=base_nodes, edges=edges,
            timesteps=np.array(TIMESTEPS), pressures=np.array(pressures),
            displacements=np.stack(disps), nu=nu, C10_true=C10_true)
        mx = np.abs(disps[-1]).max()
        print(f"wrote {out}  ({sub})  peak max|u|={mx:.4e} cm")


if __name__ == "__main__":
    main()
