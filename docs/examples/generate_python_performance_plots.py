"""Generate illustrative figures for the Python performance documentation.

Run with: uv run --no-project --with matplotlib docs/examples/generate_python_performance_plots.py
The values are illustrative, not ACTS performance measurements.
"""

from math import exp, sqrt
from pathlib import Path
from random import Random

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGURES = Path(__file__).resolve().parents[1] / "figures" / "python"


def save(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = FIGURES / name
    fig.savefig(path, bbox_inches="tight")
    path.write_text(
        "\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n"
    )
    plt.close(fig)


def main():
    plt.rcParams.update({"font.size": 11, "svg.fonttype": "none"})
    rng = Random(42)

    pt = [0.5 * i for i in range(1, 21)]
    trials = 150
    efficiency = [
        sum(rng.random() < 0.94 * (1 - exp(-p / 1.3)) for _ in range(trials)) / trials
        for p in pt
    ]
    efficiency_error = [sqrt(value * (1 - value) / trials) for value in efficiency]
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.errorbar(
        pt,
        efficiency,
        yerr=efficiency_error,
        fmt="o",
        markersize=3,
        capsize=2,
        color="#16669e",
    )
    ax.set(xlabel=r"$p_T$ [GeV]", ylabel="Tracking efficiency", ylim=(0, 1.02))
    ax.grid(alpha=0.2)
    save(fig, "tracking_efficiency.svg")

    residual = [0.03 * i for i in range(-15, 16)]
    expected = [125 * exp(-0.5 * (x / 0.09) ** 2) for x in residual]
    counts = [max(0, round(rng.gauss(value, sqrt(value)))) for value in expected]
    count_error = [sqrt(count) for count in counts]
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.step(residual, counts, where="mid", color="#a4493b")
    ax.errorbar(
        residual, counts, yerr=count_error, fmt="none", capsize=2, color="#a4493b"
    )
    ax.set(xlabel=r"$d_0$ residual [mm]", ylabel="Tracks")
    ax.grid(alpha=0.2)
    save(fig, "track_residual.svg")


if __name__ == "__main__":
    main()
