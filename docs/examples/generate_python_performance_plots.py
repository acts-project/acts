"""Generate illustrative figures for the Python performance documentation.

Run with: uv run --no-project --with matplotlib docs/examples/generate_python_performance_plots.py
The values are illustrative, not ACTS performance measurements.
"""

from math import exp
from pathlib import Path

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

    pt = [0.5 * i for i in range(1, 21)]
    efficiency = [0.94 * (1 - exp(-p / 1.3)) for p in pt]
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.plot(pt, efficiency, marker="o", markersize=3, color="#16669e")
    ax.set(xlabel=r"$p_T$ [GeV]", ylabel="Tracking efficiency", ylim=(0, 1.02))
    ax.grid(alpha=0.2)
    save(fig, "tracking_efficiency.svg")

    residual = [0.02 * i for i in range(-20, 21)]
    counts = [125 * exp(-0.5 * (x / 0.09) ** 2) for x in residual]
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.step(residual, counts, where="mid", color="#a4493b")
    ax.set(xlabel=r"$d_0$ residual [mm]", ylabel="Tracks")
    ax.grid(alpha=0.2)
    save(fig, "track_residual.svg")


if __name__ == "__main__":
    main()
