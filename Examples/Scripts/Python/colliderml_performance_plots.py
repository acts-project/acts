#!/usr/bin/env python3
"""Performance plots for ColliderML truth-tracking.

Reads histograms.pkl produced by colliderml_truth_tracking.py and writes
performance_plots.pdf into the same directory.

Usage
-----
./run_in_env.sh python3 Examples/Scripts/Python/colliderml_performance_plots.py \\
    --output colliderml_output
"""

import argparse
import datetime
import pickle
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from atlasify import atlasify
from matplotlib.backends.backend_pdf import PdfPages


def _git_meta():
    def _run(cmd):
        try:
            return subprocess.check_output(cmd, text=True).strip()
        except Exception:
            return "unknown"

    branch = _run(["git", "rev-parse", "--abbrev-ref", "HEAD"])
    commit = _run(["git", "rev-parse", "--short", "HEAD"])
    return branch, commit


def _fix_label(s: str) -> str:
    """Translate ROOT-style axis labels to matplotlib."""
    return (
        s.replace("#eta", "η")
        .replace("#phi", "φ")
        .replace("#pi", "π")
        .replace("GeV/c", "GeV")
    )


def _eff_arrays(h):
    edges = h["edges"]
    centers = 0.5 * (edges[:-1] + edges[1:])
    acc, tot = h["accepted"], h["total"]
    mask = tot > 0
    with np.errstate(invalid="ignore", divide="ignore"):
        eff = np.where(mask, acc / tot, np.nan)
        err = np.where(mask, np.sqrt(acc * (tot - acc) / tot**3), np.nan)
    return centers, eff, err, _fix_label(h["label"])


def _prof_arrays(h):
    edges = h["edges"]
    centers = 0.5 * (edges[:-1] + edges[1:])
    counts, means, sum_dq = h["counts"], h["means"], h["sum_of_deltas_squared"]
    with np.errstate(invalid="ignore", divide="ignore"):
        var = np.where(counts > 1, sum_dq / (counts - 1), 0.0)
        err = np.where(counts > 0, np.sqrt(var / np.maximum(counts, 1)), np.nan)
    means = np.where(counts > 0, means, np.nan)
    return centers, means, err, _fix_label(h["label"])


def _plot_eff(ax, hists, key, label=None, color=None, xscale="linear"):
    if key not in hists:
        ax.text(0.5, 0.5, f"(no data: {key})", transform=ax.transAxes, ha="center")
        return
    x, y, ye, xlabel = _eff_arrays(hists[key])
    kw = dict(fmt="o-", markersize=3, capsize=2, lw=1, elinewidth=0.8)
    if label is not None:
        kw["label"] = label
    if color is not None:
        kw["color"] = color
    ax.errorbar(x, y, yerr=ye, **kw)
    ax.set_xscale(xscale)
    ax.set_xlabel(xlabel)
    ax.set_ylim(0, 1.15)
    ax.axhline(1.0, color="gray", ls="--", lw=0.8)


def _plot_prof(ax, hists, key, ylabel, ymin=None):
    if key not in hists:
        ax.text(0.5, 0.5, f"(no data: {key})", transform=ax.transAxes, ha="center")
        return
    x, y, ye, xlabel = _prof_arrays(hists[key])
    ax.errorbar(x, y, yerr=ye, fmt="o-", markersize=3, capsize=2, lw=1, elinewidth=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if ymin is not None:
        ax.set_ylim(bottom=ymin)


def _atlas(ax, subtext):
    atlasify(
        "Internal", subtext=subtext, axes=ax, font_size=9, sub_font_size=8, enlarge=1.0
    )


def _add_footer(fig, footer_text):
    fig.text(
        0.5,
        0.005,
        footer_text,
        ha="center",
        va="bottom",
        fontsize=7,
        color="gray",
        transform=fig.transFigure,
    )


def make_plots(hists, output_dir: Path, n_events: int):
    branch, commit = _git_meta()
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    footer = f"{ts}  |  {branch}  |  {commit}"
    subtext = f"ColliderML ttbar PU200, {n_events} particles"

    ts_file = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    output_path = output_dir / f"performance_plots_{ts_file}.pdf"
    with PdfPages(output_path) as pdf:
        d = pdf.infodict()
        d["Title"] = "ColliderML Truth-Tracking Performance"
        d["Subject"] = "Kalman filter, ODD, ttbar PU200"
        d["Keywords"] = "ACTS ColliderML tracking performance"
        d["CreationDate"] = datetime.datetime.now()

        # ------------------------------------------------------------------
        # Slide 1 — Core performance
        # ------------------------------------------------------------------
        fig, axes = plt.subplots(2, 2, figsize=(12, 6.75))
        fig.tight_layout(pad=2.5, rect=[0, 0.03, 1, 1])

        _plot_eff(axes[0, 0], hists, "trackeff_vs_eta")
        axes[0, 0].set_ylabel("Efficiency")
        axes[0, 0].set_title("Tracking efficiency vs η")
        _atlas(axes[0, 0], subtext)

        _plot_eff(axes[0, 1], hists, "trackeff_vs_pT", xscale="log")
        axes[0, 1].set_ylabel("Efficiency")
        axes[0, 1].set_title("Tracking efficiency vs $p_T$ [GeV]")
        _atlas(axes[0, 1], subtext)

        _plot_eff(axes[1, 0], hists, "fakeRatio_vs_eta")
        axes[1, 0].set_ylabel("Fake ratio")
        axes[1, 0].set_title("Fake ratio vs η")
        _atlas(axes[1, 0], subtext)

        _plot_eff(axes[1, 1], hists, "duplicationRatio_vs_eta")
        axes[1, 1].set_ylabel("Duplication ratio")
        axes[1, 1].set_title("Duplication ratio vs η")
        _atlas(axes[1, 1], subtext)

        _add_footer(fig, footer)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # ------------------------------------------------------------------
        # Slide 2 — Hit content
        # ------------------------------------------------------------------
        fig, axes = plt.subplots(2, 2, figsize=(12, 6.75))
        fig.tight_layout(pad=2.5, rect=[0, 0.03, 1, 1])

        hit_panels = [
            ("nMeasurements_vs_eta", "Measurements / track", "Measurements vs η"),
            ("nHoles_vs_eta", "Holes / track", "Holes vs η"),
            ("nOutliers_vs_eta", "Outliers / track", "Outliers vs η"),
            ("nSharedHits_vs_eta", "Shared hits / track", "Shared hits vs η"),
        ]
        for ax, (key, ylabel, title) in zip(axes.flat, hit_panels):
            _plot_prof(ax, hists, key, ylabel=ylabel, ymin=0)
            ax.set_title(title)
            _atlas(ax, subtext)

        _add_footer(fig, footer)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # ------------------------------------------------------------------
        # Slide 3 — Track quality
        # ------------------------------------------------------------------
        fig, axes = plt.subplots(2, 2, figsize=(12, 6.75))
        fig.tight_layout(pad=2.5, rect=[0, 0.03, 1, 1])

        quality_panels = [
            ("completeness_vs_eta", "Completeness", "Completeness vs η"),
            ("purity_vs_eta", "Purity", "Purity vs η"),
            ("completeness_vs_pT", "Completeness", "Completeness vs $p_T$ [GeV]"),
            ("purity_vs_pT", "Purity", "Purity vs $p_T$ [GeV]"),
        ]
        for ax, (key, ylabel, title) in zip(axes.flat, quality_panels):
            _plot_prof(ax, hists, key, ylabel=ylabel, ymin=0)
            ax.set_title(title)
            ax.axhline(1.0, color="gray", ls="--", lw=0.8)
            ax.set_ylim(0, 1.3)  # headroom for ATLAS label above data near 1.0
            _atlas(ax, subtext)

        _add_footer(fig, footer)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # ------------------------------------------------------------------
        # Slide 4 — Differential efficiency
        # ------------------------------------------------------------------
        fig, axes = plt.subplots(1, 2, figsize=(12, 6.75))
        fig.tight_layout(pad=2.5, rect=[0, 0.03, 1, 1])

        pt_labels = ["~1 GeV", "~10 GeV", "~100 GeV"]
        for i, (color, label) in enumerate(zip(["C0", "C1", "C2"], pt_labels)):
            _plot_eff(
                axes[0], hists, f"trackeff_vs_eta_ptRange_{i}", label=label, color=color
            )
        axes[0].set_ylabel("Efficiency")
        axes[0].set_title("Efficiency vs η — $p_T$ slices")
        axes[0].legend(fontsize=9, loc="lower center")
        _atlas(axes[0], subtext)

        eta_labels = ["|η| < 0.2", "|η| < 0.8", "1 < |η| < 2", "2 < |η| < 3"]
        for i, (color, label) in enumerate(zip(["C0", "C1", "C2", "C3"], eta_labels)):
            _plot_eff(
                axes[1],
                hists,
                f"trackeff_vs_pT_absEtaRange_{i}",
                label=label,
                color=color,
                xscale="log",
            )
        axes[1].set_ylabel("Efficiency")
        axes[1].set_title("Efficiency vs $p_T$ [GeV] — η slices")
        axes[1].legend(fontsize=9, loc="lower right")
        _atlas(axes[1], subtext)

        _add_footer(fig, footer)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"Saved performance plots → {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        required=True,
        help="Output directory containing histograms.pkl (from colliderml_truth_tracking.py)",
    )
    args = parser.parse_args()

    hist_path = args.output / "histograms.pkl"
    if not hist_path.exists():
        raise FileNotFoundError(f"histograms.pkl not found in {args.output}")

    with open(hist_path, "rb") as f:
        hists = pickle.load(f)

    # infer particle count from total entries in efficiency histogram
    n_events = "?"
    if "trackeff_vs_eta" in hists:
        try:
            n_events = int(hists["trackeff_vs_eta"]["total"].sum())
        except Exception:
            pass

    make_plots(hists, args.output, n_events)


if "__main__" == __name__:
    main()
