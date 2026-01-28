#!/usr/bin/env python3

# Copyright (c) 2025 ACTS-Project
# This file is part of ACTS.
# See LICENSE for details.

import argparse
import time
from pathlib import Path

import acts
import matplotlib.pyplot as plt
import numpy as np
from acts import MagneticFieldContext
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


def create_toroidal_field():
    """Create a toroidal field with default configuration"""
    config = acts.ToroidalField.Config()
    return acts.ToroidalField(config)


def benchmark_single_evaluations(field, num_evaluations=10000):
    """
    Benchmark single magnetic field evaluations at random points.

    This function tests the performance of individual field evaluations across
    a representative sample of detector geometry positions. It measures the
    time required for single getField() calls and provides statistics on
    evaluation speed.
    """
    print("\n=== Single Field Evaluation Benchmark ===")
    print(f"Number of evaluations: {num_evaluations}")

    np.random.seed(42)  # For reproducible results

    # Generate points in cylindrical coordinates, then convert
    r = np.random.uniform(0.1, 5.0, num_evaluations)  # 0.1m to 5m radius
    phi = np.random.uniform(0, 2 * np.pi, num_evaluations)
    z = np.random.uniform(-10.0, 10.0, num_evaluations)  # -10m to +10m in z

    x = r * np.cos(phi)
    y = r * np.sin(phi)

    positions = np.column_stack((x, y, z))

    # Warm up - evaluate a few fields first
    ctx = MagneticFieldContext()
    cache = field.makeCache(ctx)
    for i in range(10):
        pos = acts.Vector3(positions[i])
        field.getField(pos, cache)

    # Benchmark single evaluations
    start_time = time.perf_counter()

    for i in range(num_evaluations):
        pos = acts.Vector3(positions[i])
        field.getField(pos, cache)

    end_time = time.perf_counter()
    total_time = end_time - start_time

    avg_time_per_eval = total_time / num_evaluations
    evaluations_per_second = num_evaluations / total_time

    print(f"Total time: {total_time:.4f} seconds")
    print(f"Average time per evaluation: " f"{avg_time_per_eval*1e6:.2f} microseconds")
    print(f"Evaluations per second: {evaluations_per_second:.0f}")

    return avg_time_per_eval, evaluations_per_second


def benchmark_vectorized_evaluations(field, batch_sizes=None):
    """Benchmark magnetic field evaluations with different batch sizes"""
    if batch_sizes is None:
        batch_sizes = [1, 10, 100, 1000, 10000]
    print("\n=== Vectorized Field Evaluation Benchmark ===")

    results = []

    for batch_size in batch_sizes:
        print(f"\nBatch size: {batch_size}")

        # Generate random test points for this batch size
        np.random.seed(42)
        r = np.random.uniform(0.1, 5.0, batch_size)
        phi = np.random.uniform(0, 2 * np.pi, batch_size)
        z = np.random.uniform(-10.0, 10.0, batch_size)

        x = r * np.cos(phi)
        y = r * np.sin(phi)
        positions = np.column_stack((x, y, z))

        # Create cache for this batch
        ctx = MagneticFieldContext()
        cache = field.makeCache(ctx)

        # Warm up
        for i in range(min(10, batch_size)):
            pos = acts.Vector3(positions[i])
            field.getField(pos, cache)

        # Benchmark this batch size
        # Adjust iterations based on batch size
        num_iterations = max(1, 1000 // batch_size)

        start_time = time.perf_counter()

        for _iteration in range(num_iterations):
            for i in range(batch_size):
                pos = acts.Vector3(positions[i])
                field.getField(pos, cache)

        end_time = time.perf_counter()

        total_evaluations = num_iterations * batch_size
        total_time = end_time - start_time
        avg_time_per_eval = total_time / total_evaluations
        evaluations_per_second = total_evaluations / total_time

        print(f"  Total evaluations: {total_evaluations}")
        print(f"  Total time: " f"{total_time:.4f} seconds")
        print(
            f"  Average time per evaluation: "
            f"{avg_time_per_eval*1e6:.2f} microseconds"
        )
        print(f"  Evaluations per second: {evaluations_per_second:.0f}")

        results.append(
            {
                "batch_size": batch_size,
                "avg_time_us": avg_time_per_eval * 1e6,
                "eval_per_sec": evaluations_per_second,
                "total_evaluations": total_evaluations,
            }
        )

    return results


def benchmark_spatial_distribution(field, grid_resolution=50):
    """Benchmark field evaluation across different spatial regions"""
    print("\n=== Spatial Distribution Benchmark ===")
    print(f"Grid resolution: {grid_resolution}x{grid_resolution} points")

    # Test different spatial regions
    regions = [
        {"name": "Barrel Toroid", "r_range": (1.0, 3.0), "z_range": (-2.0, 2.0)},
        {"name": "Forward Endcap", "r_range": (0.5, 4.0), "z_range": (2.0, 8.0)},
        {"name": "Backward Endcap", "r_range": (0.5, 4.0), "z_range": (-8.0, -2.0)},
        {"name": "Central", "r_range": (0.1, 1.0), "z_range": (-1.0, 1.0)},
    ]

    region_results = []

    for region in regions:
        print(f"\n--- {region['name']} Region ---")

        # Generate grid points in this region
        r_min, r_max = region["r_range"]
        z_min, z_max = region["z_range"]

        r_vals = np.linspace(r_min, r_max, grid_resolution)
        z_vals = np.linspace(z_min, z_max, grid_resolution)

        # Create cache for this region
        ctx = MagneticFieldContext()
        cache = field.makeCache(ctx)

        times = []
        field_magnitudes = []

        start_time = time.perf_counter()

        for r in r_vals:
            for z in z_vals:
                # Use phi=0 for consistency
                x = r
                y = 0.0

                pos = acts.Vector3(x, y, z)
                eval_start = time.perf_counter()
                b_field = field.getField(pos, cache)
                eval_end = time.perf_counter()

                times.append(eval_end - eval_start)
                field_magnitudes.append(
                    np.sqrt(b_field[0] ** 2 + b_field[1] ** 2 + b_field[2] ** 2)
                )

        end_time = time.perf_counter()

        total_evaluations = len(times)
        total_time = end_time - start_time
        avg_time = np.mean(times)
        std_time = np.std(times)
        avg_field_magnitude = np.mean(field_magnitudes)

        print(f"  Total evaluations: {total_evaluations}")
        print(f"  Total time: " f"{total_time:.4f} seconds")
        print(
            f"  Average time per evaluation: "
            f"{avg_time*1e6:.2f} ± {std_time*1e6:.2f} microseconds"
        )
        print(f"  Average field magnitude: " f"{avg_field_magnitude:.4f} Tesla")
        print(f"  Evaluations per second: " f"{total_evaluations/total_time:.0f}")

        region_results.append(
            {
                "region": region["name"],
                "total_evaluations": total_evaluations,
                "avg_time_us": avg_time * 1e6,
                "std_time_us": std_time * 1e6,
                "avg_field_magnitude": avg_field_magnitude,
                "eval_per_sec": total_evaluations / total_time,
            }
        )

    return region_results


def plot_benchmark_results(batch_results, region_results, output_dir):
    """Create plots showing benchmark results"""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    # Plot 1: Batch size performance
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    batch_sizes = [r["batch_size"] for r in batch_results]
    avg_times = [r["avg_time_us"] for r in batch_results]
    eval_rates = [r["eval_per_sec"] for r in batch_results]

    ax1.semilogx(batch_sizes, avg_times, "o-", linewidth=2, markersize=8)
    ax1.set_xlabel("Batch Size")
    ax1.set_ylabel("Average Time per Evaluation (μs)")
    ax1.set_title("Evaluation Time vs Batch Size")
    ax1.grid(True, alpha=0.3)

    ax2.semilogx(
        batch_sizes, eval_rates, "s-", linewidth=2, markersize=8, color="orange"
    )
    ax2.set_xlabel("Batch Size")
    ax2.set_ylabel("Evaluations per Second")
    ax2.set_title("Evaluation Rate vs Batch Size")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / "batch_performance.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Plot 2: Regional performance
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    regions = [r["region"] for r in region_results]
    region_times = [r["avg_time_us"] for r in region_results]
    region_rates = [r["eval_per_sec"] for r in region_results]

    bars1 = ax1.bar(
        regions, region_times, color=["skyblue", "lightgreen", "lightcoral", "gold"]
    )
    ax1.set_ylabel("Average Time per Evaluation (μs)")
    ax1.set_title("Evaluation Time by Region")
    ax1.tick_params(axis="x", rotation=45)

    # Add value labels on bars
    for bar, time_val in zip(bars1, region_times):
        height = bar.get_height()
        ax1.text(
            bar.get_x() + bar.get_width() / 2.0,
            height + height * 0.01,
            f"{time_val:.1f}μs",
            ha="center",
            va="bottom",
        )

    bars2 = ax2.bar(
        regions, region_rates, color=["skyblue", "lightgreen", "lightcoral", "gold"]
    )
    ax2.set_ylabel("Evaluations per Second")
    ax2.set_title("Evaluation Rate by Region")
    ax2.tick_params(axis="x", rotation=45)

    # Add value labels on bars
    for bar, rate_val in zip(bars2, region_rates):
        height = bar.get_height()
        ax2.text(
            bar.get_x() + bar.get_width() / 2.0,
            height + height * 0.01,
            f"{rate_val:.0f}/s",
            ha="center",
            va="bottom",
        )

    plt.tight_layout()
    plt.savefig(output_dir / "regional_performance.png", dpi=150, bbox_inches="tight")
    plt.close()

    print(f"\nPlots saved to {output_dir}/")


def _eval_field_batch(field, points):
    """Evaluate B-field for an array of points (N,3) -> (N,3)."""
    ctx = MagneticFieldContext()
    cache = field.makeCache(ctx)
    out = np.empty_like(points, dtype=np.float64)
    for i in range(points.shape[0]):
        b_field = field.getField(acts.Vector3(points[i]), cache)
        out[i, 0] = b_field[0]
        out[i, 1] = b_field[1]
        out[i, 2] = b_field[2]
    return out


def plot_field_maps(
    field,
    output_dir,
    # XY slice params
    xy_z_plane=0.20,
    xy_xlim=(-10.0, 10.0),
    xy_ylim=(-10.0, 10.0),
    xy_nx=520,
    xy_ny=520,
    # ZX slice params
    zx_y_plane=0.10,
    zx_zlim=(-22.0, 22.0),
    zx_xlim=(-11.0, 11.0),
    zx_nz=560,
    zx_nx=560,
    # Viz params
    log_vmin=1e-4,
    log_vmax=4.1,
    quiver_stride_xy=28,  # ~ how many arrows across
    quiver_stride_zx=28,
):
    """
    Produce two figures:
      1) XY slice at fixed z=xy_z_plane
      2) ZX slice at fixed y=zx_y_plane
    Saved to output_dir as field_xy.png and field_zx.png
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    # ---------- XY slice ----------
    x = np.linspace(xy_xlim[0], xy_xlim[1], int(max(60, xy_nx)), dtype=np.float64)
    y = np.linspace(xy_ylim[0], xy_ylim[1], int(max(60, xy_ny)), dtype=np.float64)
    X, Y = np.meshgrid(x, y, indexing="xy")
    Z = np.full_like(X, float(xy_z_plane), dtype=np.float64)

    pts_xy = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    B_xy = _eval_field_batch(field, pts_xy).reshape(*X.shape, 3)
    Bx, By, Bz = B_xy[..., 0], B_xy[..., 1], B_xy[..., 2]
    Bmag = np.sqrt(Bx**2 + By**2 + Bz**2, dtype=np.float64)

    Bpos = np.clip(Bmag, log_vmin * 1e-2, None)  # avoid zeros on log scale
    norm = LogNorm(vmin=float(log_vmin), vmax=float(log_vmax))

    fig, ax = plt.subplots(figsize=(8.8, 8.8))
    im = ax.imshow(
        Bpos,
        extent=[x.min(), x.max(), y.min(), y.max()],
        origin="lower",
        aspect="equal",
        norm=norm,
        cmap="gnuplot2",
    )
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("|B| [T] (log scale)")

    # Quiver (unit-length for direction)
    step = max(1, X.shape[1] // quiver_stride_xy)
    Bsafe = np.where(Bmag > 0, Bmag, np.inf)
    Ux = Bx / Bsafe
    Uy = By / Bsafe
    ax.quiver(
        X[::step, ::step],
        Y[::step, ::step],
        Ux[::step, ::step],
        Uy[::step, ::step],
        scale=28,
        width=0.003,
    )

    ax.set_xlim(xy_xlim)
    ax.set_ylim(xy_ylim)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(f"|B| in z = {float(xy_z_plane):.2f} m plane")
    plt.tight_layout()
    (output_dir / "field_xy.png").unlink(missing_ok=True)
    plt.savefig(output_dir / "field_xy.png", dpi=150, bbox_inches="tight")
    plt.close()

    # ---------- ZX slice (Z horizontal, X vertical) ----------
    z = np.linspace(zx_zlim[0], zx_zlim[1], int(max(60, zx_nz)), dtype=np.float64)
    x = np.linspace(zx_xlim[0], zx_xlim[1], int(max(60, zx_nx)), dtype=np.float64)
    Z, Xg = np.meshgrid(z, x, indexing="xy")
    Y = np.full_like(Xg, float(zx_y_plane), dtype=np.float64)

    pts_zx = np.column_stack([Xg.ravel(), Y.ravel(), Z.ravel()])
    B_zx = _eval_field_batch(field, pts_zx).reshape(*Xg.shape, 3)
    Bx, By, Bz = B_zx[..., 0], B_zx[..., 1], B_zx[..., 2]
    Bmag = np.sqrt(Bx**2 + By**2 + Bz**2, dtype=np.float64)

    Bpos = np.clip(Bmag, log_vmin * 1e-2, None)
    norm = LogNorm(vmin=float(log_vmin), vmax=float(log_vmax))

    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(
        Bpos,
        extent=[z.min(), z.max(), x.min(), x.max()],
        origin="lower",
        aspect="equal",
        norm=norm,
        cmap="gnuplot2",
    )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.10)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label("|B| [T] (log scale)")

    # Quiver (projected in z–x plane)
    step = max(1, Z.shape[1] // quiver_stride_zx)
    Bsafe = np.where(Bmag > 0, Bmag, np.inf)
    Uz = Bz / Bsafe
    Ux = Bx / Bsafe
    ax.quiver(
        Z[::step, ::step],
        Xg[::step, ::step],
        Uz[::step, ::step],
        Ux[::step, ::step],
        scale=28,
        width=0.003,
    )

    ax.set_xlim(zx_zlim)
    ax.set_ylim(zx_xlim)
    ax.set_xlabel("z [m]")
    ax.set_ylabel("x [m]")
    ax.set_title(f"|B| in y = {float(zx_y_plane):.2f} m plane")
    plt.tight_layout()
    (output_dir / "field_zx.png").unlink(missing_ok=True)
    plt.savefig(output_dir / "field_zx.png", dpi=150, bbox_inches="tight")
    plt.close()

    print(
        f"\nSaved field maps to: {output_dir}/field_xy.png and {output_dir}/field_zx.png"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark toroidal magnetic field performance"
    )
    parser.add_argument(
        "--num-evaluations",
        type=int,
        default=10000,
        help="Number of evaluations for single benchmark (default: 10000)",
    )
    parser.add_argument(
        "--batch-sizes",
        type=int,
        nargs="+",
        default=[1, 10, 100, 1000, 10000],
        help="Batch sizes to test (default: 1 10 100 1000 10000)",
    )
    parser.add_argument(
        "--grid-resolution",
        type=int,
        default=50,
        help="Grid resolution for spatial benchmark (default: 50)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="toroidal_field_benchmark",
        help="Output directory for results (default: toroidal_field_benchmark)",
    )
    parser.add_argument(
        "--skip-plots", action="store_true", help="Skip generating performance plots"
    )
    parser.add_argument(
        "--plot-field",
        action="store_true",
        help="Also compute and save XY/ZX field maps",
    )

    # XY slice controls
    parser.add_argument(
        "--xy-z-plane",
        type=float,
        default=0.20,
        help="z (m) for XY slice (default: 0.20)",
    )
    parser.add_argument(
        "--xy-xlim",
        type=float,
        nargs=2,
        default=(-10.0, 10.0),
        help="x limits for XY slice (m) (default: -10 10)",
    )
    parser.add_argument(
        "--xy-ylim",
        type=float,
        nargs=2,
        default=(-10.0, 10.0),
        help="y limits for XY slice (m) (default: -10 10)",
    )
    parser.add_argument(
        "--xy-nx", type=int, default=520, help="grid Nx for XY slice (default: 520)"
    )
    parser.add_argument(
        "--xy-ny", type=int, default=520, help="grid Ny for XY slice (default: 520)"
    )

    # ZX slice controls
    parser.add_argument(
        "--zx-y-plane",
        type=float,
        default=0.10,
        help="y (m) for ZX slice (default: 0.10)",
    )
    parser.add_argument(
        "--zx-zlim",
        type=float,
        nargs=2,
        default=(-22.0, 22.0),
        help="z limits for ZX slice (m) (default: -22 22)",
    )
    parser.add_argument(
        "--zx-xlim",
        type=float,
        nargs=2,
        default=(-11.0, 11.0),
        help="x limits for ZX slice (m) (default: -11 11)",
    )
    parser.add_argument(
        "--zx-nz", type=int, default=560, help="grid Nz for ZX slice (default: 560)"
    )
    parser.add_argument(
        "--zx-nx", type=int, default=560, help="grid Nx for ZX slice (default: 560)"
    )

    # Color/arrow controls
    parser.add_argument(
        "--log-vmin",
        type=float,
        default=1e-4,
        help="LogNorm vmin for |B| (T) (default: 1e-4)",
    )
    parser.add_argument(
        "--log-vmax",
        type=float,
        default=4.1,
        help="LogNorm vmax for |B| (T) (default: 4.1)",
    )
    parser.add_argument(
        "--quiver-stride-xy",
        type=int,
        default=28,
        help="Quiver density for XY (default: 28)",
    )
    parser.add_argument(
        "--quiver-stride-zx",
        type=int,
        default=28,
        help="Quiver density for ZX (default: 28)",
    )

    args = parser.parse_args()

    print("=== Toroidal Magnetic Field Benchmark ===")
    print(f"ACTS version: {acts.version}")

    # Create toroidal field
    print("\nCreating toroidal field...")
    field = create_toroidal_field()
    print("✓ Toroidal field created successfully")

    # Run benchmarks
    try:
        # Single evaluation benchmark
        single_time, single_rate = benchmark_single_evaluations(
            field, args.num_evaluations
        )

        # Batch evaluation benchmark
        batch_results = benchmark_vectorized_evaluations(field, args.batch_sizes)

        # Spatial distribution benchmark
        region_results = benchmark_spatial_distribution(field, args.grid_resolution)

        # Print summary
        print("\n=== Benchmark Summary ===")
        print(
            f"Single evaluation average: {single_time*1e6:.2f} μs ({single_rate:.0f} eval/s)"
        )
        print(
            f"Best batch performance: {min(r['avg_time_us'] for r in batch_results):.2f} μs"
        )
        print(
            f"Fastest region: {min(region_results, key=lambda x: x['avg_time_us'])['region']}"
        )
        print(
            f"Slowest region: {max(region_results, key=lambda x: x['avg_time_us'])['region']}"
        )

        # Generate plots if requested
        if not args.skip_plots:
            try:
                plot_benchmark_results(batch_results, region_results, args.output_dir)
            except ImportError:
                print("\nWarning: matplotlib not available, skipping performance plots")

        # Field maps (XY, ZX) using ACTS field.getField
        if args.plot_field:
            try:
                plot_field_maps(
                    field,
                    output_dir=args.output_dir,
                    xy_z_plane=args.xy_z_plane,
                    xy_xlim=tuple(args.xy_xlim),
                    xy_ylim=tuple(args.xy_ylim),
                    xy_nx=args.xy_nx,
                    xy_ny=args.xy_ny,
                    zx_y_plane=args.zx_y_plane,
                    zx_zlim=tuple(args.zx_zlim),
                    zx_xlim=tuple(args.zx_xlim),
                    zx_nz=args.zx_nz,
                    zx_nx=args.zx_nx,
                    log_vmin=args.log_vmin,
                    log_vmax=args.log_vmax,
                    quiver_stride_xy=args.quiver_stride_xy,
                    quiver_stride_zx=args.quiver_stride_zx,
                )
            except ImportError:
                print(
                    "\nWarning: matplotlib (extras) not available, skipping field maps"
                )

        print("\n✓ Benchmark completed successfully!")

    except Exception as e:
        print(f"\n❌ Benchmark failed: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
