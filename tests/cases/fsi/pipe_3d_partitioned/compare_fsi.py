#!/usr/bin/env python3
"""Compare monolithic vs partitioned FSI results.

Usage:
    python compare_fsi.py [--step STEP] [--mono DIR] [--part DIR]

Plots:
    1. Flow rate at inlet over time
    2. Radial displacement at solid inlet face over time
    3. Radial displacement of solid along z-axis at given time step
    4. Centerline pressure along z-axis at given time step
    5. Centerline velocity along z-axis at given time step
"""

import argparse
import glob
import os
import numpy as np
import matplotlib.pyplot as plt

try:
    import meshio
except ImportError:
    import subprocess
    subprocess.check_call(["pip3", "install", "meshio", "-q"])
    import meshio


def find_result_files(result_dir, prefix):
    """Find all result VTU files sorted by step number."""
    pattern = os.path.join(result_dir, f"{prefix}_*.vtu")
    files = sorted(glob.glob(pattern), key=lambda f: int(f.split("_")[-1].split(".")[0]))
    return files


def get_step_number(filename):
    return int(filename.split("_")[-1].split(".")[0])


def read_mesh(filename):
    return meshio.read(filename)


def compute_flow_rate(mesh, face_name=None):
    """Integrate velocity dot normal over inlet face.
    For simplicity, sum axial velocity * nodal area at z=0."""
    pts = mesh.points
    vel = mesh.point_data.get("Velocity", None)
    if vel is None:
        return 0.0

    # Find inlet nodes (z ≈ 0)
    z = pts[:, 2]
    z_min = z.min()
    inlet_mask = np.abs(z - z_min) < 1e-6 * (z.max() - z.min() + 1e-30)

    if inlet_mask.sum() == 0:
        return 0.0

    # Axial velocity (z-component) at inlet, averaged
    v_axial = vel[inlet_mask, 2]
    return np.mean(v_axial)


def compute_radial_disp_at_inlet(mesh, disp_key="Displacement"):
    """Mean radial displacement of the outermost ring at the inlet face."""
    pts = mesh.points
    disp = mesh.point_data.get(disp_key, None)
    if disp is None:
        return 0.0

    z = pts[:, 2]
    z_min = z.min()
    z_tol = 1e-6 * (z.max() - z.min() + 1e-30)
    inlet_mask = np.abs(z - z_min) < z_tol

    if inlet_mask.sum() == 0:
        return 0.0

    # Find the outermost ring: nodes at the maximum radius on the inlet face
    r = np.sqrt(pts[inlet_mask, 0]**2 + pts[inlet_mask, 1]**2)
    r_max = r.max()
    outer_ring = r > 0.9 * r_max  # outermost 10% of radial range

    x = pts[inlet_mask, 0][outer_ring]
    y = pts[inlet_mask, 1][outer_ring]
    dx = disp[inlet_mask, 0][outer_ring]
    dy = disp[inlet_mask, 1][outer_ring]
    rr = np.sqrt(x**2 + y**2)
    rr[rr < 1e-30] = 1e-30
    radial = (x * dx + y * dy) / rr
    return np.mean(radial)


def extract_centerline(mesh, field_name, nsd=3):
    """Extract field values along the centerline (r ≈ 0)."""
    pts = mesh.points
    data = mesh.point_data.get(field_name, None)
    if data is None:
        return np.array([]), np.array([])

    r = np.sqrt(pts[:, 0]**2 + pts[:, 1]**2)
    r_max = r.max()
    cl_mask = r < 0.1 * r_max  # within 10% of center

    if cl_mask.sum() == 0:
        return np.array([]), np.array([])

    z = pts[cl_mask, 2]
    vals = data[cl_mask]
    order = np.argsort(z)
    return z[order], vals[order]


def extract_solid_radial_disp_along_z(mesh, disp_key="Displacement"):
    """Extract radial displacement along z at the outer wall ring."""
    pts = mesh.points
    disp = mesh.point_data.get(disp_key, None)
    if disp is None:
        return np.array([]), np.array([])

    # Find outermost nodes (largest radius in the mesh)
    r = np.sqrt(pts[:, 0]**2 + pts[:, 1]**2)
    r_max = r.max()
    inner_mask = r > 0.9 * r_max

    if inner_mask.sum() == 0:
        return np.array([]), np.array([])

    z = pts[inner_mask, 2]
    x, y = pts[inner_mask, 0], pts[inner_mask, 1]
    dx, dy = disp[inner_mask, 0], disp[inner_mask, 1]
    rr = np.sqrt(x**2 + y**2)
    rr[rr < 1e-30] = 1e-30
    radial = (x * dx + y * dy) / rr

    # Average over circumference at each z
    z_unique = np.unique(np.round(z, 8))
    z_avg = []
    d_avg = []
    for zz in sorted(z_unique):
        mask = np.abs(z - zz) < 1e-6
        z_avg.append(zz)
        d_avg.append(np.mean(radial[mask]))

    return np.array(z_avg), np.array(d_avg)


def main():
    parser = argparse.ArgumentParser(description="Compare monolithic vs partitioned FSI")
    parser.add_argument("--step", type=int, default=5, help="Time step for z-axis plots")
    parser.add_argument("--mono", default="../pipe_3d/1-procs", help="Monolithic results directory")
    parser.add_argument("--part", default="1-procs", help="Partitioned results directory")
    parser.add_argument("--dt", type=float, default=1e-4, help="Time step size")
    args = parser.parse_args()

    plt.rcParams.update({"font.size": 10, "figure.figsize": (8, 5)})

    # Find result files
    mono_files = find_result_files(args.mono, "result")
    part_fluid_files = find_result_files(args.part, "result_fluid")
    part_solid_files = find_result_files(args.part, "result_solid")

    if not mono_files:
        print(f"No monolithic results found in {args.mono}")
        return
    if not part_fluid_files:
        print(f"No partitioned fluid results found in {args.part}")
        return

    print(f"Monolithic: {len(mono_files)} steps")
    print(f"Partitioned: {len(part_fluid_files)} fluid, {len(part_solid_files)} solid steps")

    # ---- Time history plots ----
    n_steps = min(len(mono_files), len(part_fluid_files))
    steps = []
    mono_flow, part_flow = [], []
    mono_disp_inlet, part_disp_inlet = [], []

    for i in range(n_steps):
        step = get_step_number(mono_files[i])
        steps.append(step)

        m = read_mesh(mono_files[i])
        mono_flow.append(compute_flow_rate(m))
        mono_disp_inlet.append(compute_radial_disp_at_inlet(m, "FS_Displacement"))

        pf = read_mesh(part_fluid_files[i])
        part_flow.append(compute_flow_rate(pf))

        if i < len(part_solid_files):
            ps = read_mesh(part_solid_files[i])
            part_disp_inlet.append(compute_radial_disp_at_inlet(ps, "Displacement"))
        else:
            part_disp_inlet.append(0.0)

    time = np.array(steps) * args.dt

    # Plot 1: Flow rate at inlet
    fig, ax = plt.subplots()
    ax.plot(time * 1e3, mono_flow, "b-o", ms=3, label="Monolithic")
    ax.plot(time * 1e3, part_flow, "r-s", ms=3, label="Partitioned")
    ax.set_xlabel("Time [ms]")
    ax.set_ylabel("Mean axial velocity at inlet")
    ax.set_title("Inlet flow rate")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("compare_flow_rate.pdf")
    print("Saved compare_flow_rate.pdf")

    # Plot 2: Radial displacement at solid inlet
    fig, ax = plt.subplots()
    ax.plot(time * 1e3, mono_disp_inlet, "b-o", ms=3, label="Monolithic")
    ax.plot(time * 1e3, part_disp_inlet, "r-s", ms=3, label="Partitioned")
    ax.set_xlabel("Time [ms]")
    ax.set_ylabel("Mean radial displacement at inlet")
    ax.set_title("Solid displacement at inlet face")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("compare_disp_inlet.pdf")
    print("Saved compare_disp_inlet.pdf")

    # ---- Spatial plots at given time step ----
    step = args.step
    mono_file = os.path.join(args.mono, f"result_{step:03d}.vtu")
    part_fluid_file = os.path.join(args.part, f"result_fluid_{step:03d}.vtu")
    part_solid_file = os.path.join(args.part, f"result_solid_{step:03d}.vtu")

    if not os.path.exists(mono_file):
        print(f"Monolithic result not found: {mono_file}")
        return
    if not os.path.exists(part_fluid_file):
        print(f"Partitioned fluid result not found: {part_fluid_file}")
        return

    m = read_mesh(mono_file)

    # Plot 3: Radial displacement along z
    fig, ax = plt.subplots()
    mz, md = extract_solid_radial_disp_along_z(m, "FS_Displacement")
    ax.plot(mz, md, "b-o", ms=3, label="Monolithic")
    if os.path.exists(part_solid_file):
        ps = read_mesh(part_solid_file)
        pz, pd = extract_solid_radial_disp_along_z(ps, "Displacement")
        ax.plot(pz, pd, "r-s", ms=3, label="Partitioned")
    ax.set_xlabel("z")
    ax.set_ylabel("Radial displacement")
    ax.set_title(f"Solid radial displacement along z (step {step})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("compare_disp_z.pdf")
    print("Saved compare_disp_z.pdf")

    # Plot 4: Centerline pressure
    fig, ax = plt.subplots()
    mz, mp = extract_centerline(m, "Pressure")
    pf = read_mesh(part_fluid_file)
    pz, pp = extract_centerline(pf, "Pressure")
    if len(mz) > 0:
        ax.plot(mz, mp, "b-o", ms=3, label="Monolithic")
    if len(pz) > 0:
        ax.plot(pz, pp, "r-s", ms=3, label="Partitioned")
    ax.set_xlabel("z")
    ax.set_ylabel("Pressure")
    ax.set_title(f"Centerline pressure (step {step})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("compare_pressure_z.pdf")
    print("Saved compare_pressure_z.pdf")

    # Plot 5: Centerline axial velocity
    fig, ax = plt.subplots()
    mz, mv = extract_centerline(m, "Velocity")
    pz, pv = extract_centerline(pf, "Velocity")
    if len(mz) > 0:
        ax.plot(mz, mv[:, 2] if mv.ndim > 1 else mv, "b-o", ms=3, label="Monolithic")
    if len(pz) > 0:
        ax.plot(pz, pv[:, 2] if pv.ndim > 1 else pv, "r-s", ms=3, label="Partitioned")
    ax.set_xlabel("z")
    ax.set_ylabel("Axial velocity")
    ax.set_title(f"Centerline axial velocity (step {step})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("compare_velocity_z.pdf")
    print("Saved compare_velocity_z.pdf")

    plt.close("all")
    print("Done.")


if __name__ == "__main__":
    main()
