#!/usr/bin/env python3
"""
Compare ParaDiS subcycling and trapezoid results.

Expected folder structure:

subcycling_folder/
    time_Plastic_strain
    density

trapezoid_folder/
    time_Plastic_strain
    density

Usage:
    python compare_integrators.py /path/to/subcycling /path/to/trapezoid

Optional:
    python compare_integrators.py /path/to/subcycling /path/to/trapezoid \
        --output comparison_plots
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def load_plastic_strain(folder):
    """
    time_Plastic_strain format:
        column 1 = simulation time
        column 2 = plastic strain
    """
    path = Path(folder) / "time_Plastic_strain"
    data = np.loadtxt(path)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 2:
        raise ValueError(f"{path} must contain at least 2 columns.")

    time = data[:, 0]
    strain = data[:, 1]

    order = np.argsort(time)
    return time[order], strain[order]


def load_density(folder):
    """
    ParaDiS density format used here:
        column 2 = simulation time
        column 3 = total dislocation density

    Example:
        plastic_strain  time  density  ...
    """
    path = Path(folder) / "density"
    data = np.loadtxt(path)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 3:
        raise ValueError(f"{path} must contain at least 3 columns.")

    time = data[:, 1]
    density = data[:, 2]

    order = np.argsort(time)
    return time[order], density[order]


def remove_duplicate_times(time, values):
    """Keep the last value for duplicate time entries."""
    unique_time, reverse_index = np.unique(time[::-1], return_index=True)
    index = len(time) - 1 - reverse_index

    order = np.argsort(unique_time)
    return unique_time[order], values[index][order]


def common_grid(t1, t2):
    """
    Construct a comparison grid over the overlapping physical-time interval.

    Uses all available time points from both calculations within the overlap.
    """
    tmin = max(np.min(t1), np.min(t2))
    tmax = min(np.max(t1), np.max(t2))

    if tmax <= tmin:
        raise ValueError(
            "The two calculations have no overlapping physical-time interval."
        )

    a = t1[(t1 >= tmin) & (t1 <= tmax)]
    b = t2[(t2 >= tmin) & (t2 <= tmax)]

    return np.unique(np.concatenate((a, b)))


def relative_error(test, reference, floor=None):
    """
    Absolute relative difference:
        |test-reference| / |reference| * 100

    Very small reference values are returned as NaN to avoid meaningless
    percentages near zero.
    """
    reference_abs = np.abs(reference)

    if floor is None:
        finite_nonzero = reference_abs[reference_abs > 0.0]
        if len(finite_nonzero) == 0:
            floor = 0.0
        else:
            floor = max(
                1.0e-30,
                1.0e-8 * np.max(finite_nonzero)
            )

    error = np.full_like(reference, np.nan, dtype=float)
    mask = reference_abs > floor
    error[mask] = (
        np.abs(test[mask] - reference[mask])
        / reference_abs[mask]
        * 100.0
    )

    return error


def save_line_plot(
    x1, y1, x2, y2,
    label1, label2,
    xlabel, ylabel, title,
    output_file
):
    plt.figure(figsize=(8, 5.5))
    plt.plot(x1, y1, label=label1, linewidth=1.8)
    plt.plot(x2, y2, label=label2, linewidth=1.8)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def save_error_plot(time, error, ylabel, title, output_file):
    plt.figure(figsize=(8, 5.5))
    plt.plot(time, error, linewidth=1.8)
    plt.xlabel("Simulation time")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Compare ParaDiS subcycling and trapezoid outputs."
    )
    parser.add_argument(
        "subcycling_folder",
        help="Folder containing subcycling time_Plastic_strain and density files.",
    )
    parser.add_argument(
        "trapezoid_folder",
        help="Folder containing trapezoid time_Plastic_strain and density files.",
    )
    parser.add_argument(
        "--output",
        default="integrator_comparison",
        help="Output directory. Default: integrator_comparison",
    )

    args = parser.parse_args()

    sub_folder = Path(args.subcycling_folder)
    trap_folder = Path(args.trapezoid_folder)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------
    # Read data
    # ------------------------------------------------------------

    ts, eps_s = load_plastic_strain(sub_folder)
    tt, eps_t = load_plastic_strain(trap_folder)

    td_s, rho_s = load_density(sub_folder)
    td_t, rho_t = load_density(trap_folder)

    ts, eps_s = remove_duplicate_times(ts, eps_s)
    tt, eps_t = remove_duplicate_times(tt, eps_t)
    td_s, rho_s = remove_duplicate_times(td_s, rho_s)
    td_t, rho_t = remove_duplicate_times(td_t, rho_t)

    # ------------------------------------------------------------
    # Plot raw trajectories
    # ------------------------------------------------------------

    save_line_plot(
        ts, eps_s,
        tt, eps_t,
        "Subcycling", "Trapezoid",
        "Simulation time",
        "Plastic strain",
        "Plastic strain: subcycling vs trapezoid",
        output / "plastic_strain_comparison.png",
    )

    save_line_plot(
        td_s, rho_s,
        td_t, rho_t,
        "Subcycling", "Trapezoid",
        "Simulation time",
        r"Dislocation density (m$^{-2}$)",
        "Dislocation density: subcycling vs trapezoid",
        output / "density_comparison.png",
    )

    # ------------------------------------------------------------
    # Compare at common physical times
    # ------------------------------------------------------------

    strain_time = common_grid(ts, tt)

    strain_sub = np.interp(strain_time, ts, eps_s)
    strain_trap = np.interp(strain_time, tt, eps_t)

    strain_error = relative_error(
        strain_sub,
        strain_trap,
    )

    density_time = common_grid(td_s, td_t)

    density_sub = np.interp(density_time, td_s, rho_s)
    density_trap = np.interp(density_time, td_t, rho_t)

    density_error = relative_error(
        density_sub,
        density_trap,
    )

    # ------------------------------------------------------------
    # Error plots
    # ------------------------------------------------------------

    save_error_plot(
        strain_time,
        strain_error,
        "Relative difference (%)",
        "Plastic strain relative difference",
        output / "plastic_strain_relative_error.png",
    )

    save_error_plot(
        density_time,
        density_error,
        "Relative difference (%)",
        "Density relative difference",
        output / "density_relative_error.png",
    )

    # ------------------------------------------------------------
    # Save interpolated comparison data
    # ------------------------------------------------------------

    strain_table = np.column_stack(
        (
            strain_time,
            strain_sub,
            strain_trap,
            strain_error,
        )
    )

    np.savetxt(
        output / "plastic_strain_comparison.csv",
        strain_table,
        delimiter=",",
        header=(
            "time,"
            "subcycling_plastic_strain,"
            "trapezoid_plastic_strain,"
            "relative_difference_percent"
        ),
        comments="",
    )

    density_table = np.column_stack(
        (
            density_time,
            density_sub,
            density_trap,
            density_error,
        )
    )

    np.savetxt(
        output / "density_comparison.csv",
        density_table,
        delimiter=",",
        header=(
            "time,"
            "subcycling_density,"
            "trapezoid_density,"
            "relative_difference_percent"
        ),
        comments="",
    )

    # ------------------------------------------------------------
    # Print useful summary
    # ------------------------------------------------------------

    print()
    print("Comparison complete")
    print("-------------------")

    print(
        f"Subcycling strain time range : "
        f"{ts.min():.6e} -> {ts.max():.6e}"
    )
    print(
        f"Trapezoid strain time range  : "
        f"{tt.min():.6e} -> {tt.max():.6e}"
    )

    print(
        f"Common strain interval        : "
        f"{strain_time.min():.6e} -> {strain_time.max():.6e}"
    )

    finite_strain = strain_error[np.isfinite(strain_error)]
    if len(finite_strain):
        print(
            f"Mean strain difference        : "
            f"{np.mean(finite_strain):.3f} %"
        )
        print(
            f"Maximum strain difference     : "
            f"{np.max(finite_strain):.3f} %"
        )

    finite_density = density_error[np.isfinite(density_error)]
    if len(finite_density):
        print(
            f"Mean density difference       : "
            f"{np.mean(finite_density):.3f} %"
        )
        print(
            f"Maximum density difference    : "
            f"{np.max(finite_density):.3f} %"
        )

    print()
    print(f"Plots and CSV files written to: {output.resolve()}")
    print()
    print("Generated files:")
    print("  plastic_strain_comparison.png")
    print("  density_comparison.png")
    print("  plastic_strain_relative_error.png")
    print("  density_relative_error.png")
    print("  plastic_strain_comparison.csv")
    print("  density_comparison.csv")


if __name__ == "__main__":
    main()
