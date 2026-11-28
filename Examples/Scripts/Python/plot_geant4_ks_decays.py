"""Plot Geant4 K-short transport checks from geant4_ks_candidates.py CSV files.

The gun input has equal-sized samples at fixed, nonzero parent momenta.
The ODD input contains K-shorts in hadronic Z decays. All pion momenta are
transport truth at the daughter production vertex. Units are in the CSV header.
"""

import argparse
import csv
from itertools import cycle
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator
import numpy as np

# PDG 2025 strange- and light-meson tables, masses in MeV and lengths in mm.
KAON_MASS = 497.611
PION_MASS = 139.57039
CTAU = 26.844
INK = "#203349"
TEAL = "#00796d"
BLUE = "#365fa7"
ORANGE = "#ad521e"


def read_candidates(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def values(rows, field):
    return np.array([float(row[field]) for row in rows])


def flag(row, field):
    return row[field].lower() in {"true", "1"}


def figure(title, subtitle):
    fig = plt.figure(figsize=(11.8, 5.7))
    fig.text(0.075, 0.957, title, fontsize=19, weight="bold", va="top")
    fig.text(0.075, 0.888, subtitle, fontsize=10.5, va="top", color="#475769")
    return fig


def panel(ax, title):
    ax.set_title(title, loc="left", fontsize=12, weight="bold", pad=9)
    ax.grid(axis="y", alpha=0.17)
    ax.set_axisbelow(True)


def reconstruction(fits, odd):
    fig = figure(
        r"A $K_S^0$ mass peak from fitted pion tracks",
        f"{len(fits):,} truth-matched, opposite-charge CKF pairs  ·  OpenDataDetector in 2 T  ·  No mass or decay-radius cut",
    )
    gs = fig.add_gridspec(
        1,
        2,
        left=0.08,
        right=0.965,
        top=0.755,
        bottom=0.17,
        width_ratios=[1.18, 1],
        wspace=0.31,
    )
    mass, vertex = [fig.add_subplot(gs[i]) for i in range(2)]
    measured = values(fits, "m_fit_at_dca_mev")
    # A 1 MeV/c² bin width resolves the peak while retaining every candidate.
    edges = np.arange(np.floor(measured.min()) - 1, np.ceil(measured.max()) + 2)
    mass.hist(measured, bins=edges, color=BLUE, alpha=0.90, label="Fitted pion pairs")
    mass.axvline(KAON_MASS, color="#8b301e", lw=1.8, ls="--", label=r"PDG $K^0$ mass")
    lo, med, hi = np.quantile(measured, [0.16, 0.50, 0.84])
    panel(mass, "Dipion mass at two-track closest approach")
    mass.set(
        xlabel=r"Invariant mass  $M_{\pi^+\pi^-}$  [MeV/$c^2$]",
        ylabel=r"Pairs / (1 MeV/$c^2$)",
    )
    mass.yaxis.set_major_locator(MaxNLocator(integer=True))
    mass.legend(loc="upper right", frameon=False, fontsize=9)
    mass.text(
        0.97,
        0.73,
        f"Median {med:.2f}\nCentral 68%\n{lo:.2f}–{hi:.2f} MeV/$c^2$",
        transform=mass.transAxes,
        ha="right",
        va="top",
        fontsize=9,
    )

    truth = values(fits, "r_decay_mm")
    fitted = values(fits, "fit_r_mm")
    limit = 1.05 * max(truth.max(), fitted.max())
    vertex.plot([0, limit], [0, limit], color=INK, lw=1.3, ls="--")
    vertex.scatter(
        truth, fitted, s=24, color=TEAL, alpha=0.65, edgecolors="white", linewidth=0.35
    )
    panel(vertex, "Transverse decay radius")
    vertex.set(
        xlabel="Geant4 decay radius  [mm]",
        ylabel="Two-track decay radius  [mm]",
        xlim=(0, limit),
        ylim=(0, limit),
        aspect="equal",
    )
    vertex.text(
        0.04,
        0.94,
        "Median 3D vertex error\n"
        f"{np.median(values(fits, 'fit_xyz_err_mm')):.2f} mm",
        transform=vertex.transAxes,
        va="top",
        fontsize=9.5,
    )
    ndecays = sum(flag(row, "pair") for row in odd)
    fig.text(
        0.08,
        0.027,
        f"Both daughter tracks were uniquely matched for {len(fits):,} of {ndecays:,} simulated decays. "
        "The vertex uses the two fitted helices alone.",
        fontsize=9,
        color="#475769",
    )
    return fig


def kinematics(gun):
    pairs = [row for row in gun if flag(row, "pair")]
    fig = figure(
        r"Two-body kinematics of $K_S^0 \rightarrow \pi^+\pi^-$",
        f"{len(pairs):,} paired decays  ·  Pythia8 kaon gun  ·  Geant4 daughter momenta at creation",
    )
    grid = fig.add_gridspec(
        2,
        2,
        left=0.075,
        right=0.970,
        bottom=0.12,
        top=0.765,
        width_ratios=[1.1, 1],
        hspace=0.98,
        wspace=0.33,
    )
    locus = fig.add_subplot(grid[:, 0])
    mass = fig.add_subplot(grid[0, 1])
    angle = fig.add_subplot(grid[1, 1])
    parent_p = values(pairs, "gen_p_gev").round(6)
    alpha = values(pairs, "g4_alpha")
    qt = values(pairs, "g4_qt_mev")
    pstar = np.sqrt((KAON_MASS / 2) ** 2 - PION_MASS**2)
    cos = np.linspace(-1, 1, 600)
    styles = cycle([(TEAL, "o"), (BLUE, "s"), (ORANGE, "^")])
    for p, (color, marker) in zip(np.unique(parent_p), styles):
        beta = p / np.hypot(p, KAON_MASS / 1000)
        locus.plot(
            2 * pstar / (beta * KAON_MASS) * cos,
            pstar * np.sqrt(1 - cos**2),
            color=color,
            lw=1.6,
        )
        # Draw a subset so individual simulated decays remain distinguishable.
        selected = np.flatnonzero(parent_p == p)
        selected = selected[:: max(1, len(selected) // 55)]
        locus.scatter(
            alpha[selected],
            qt[selected],
            s=23,
            marker=marker,
            facecolor="white",
            edgecolor=color,
            linewidth=0.9,
            label=f"{p:g} GeV/$c$",
            zorder=3,
        )
    panel(locus, "Pion momenta in the laboratory")
    locus.legend(title="Kaon momentum", loc="lower center", frameon=False)
    locus.set(
        xlabel=r"Longitudinal asymmetry  $\alpha=(p_L^+-p_L^-)/(p_L^++p_L^-)$",
        ylabel=r"Pion momentum transverse to kaon  $q_T$  [MeV/$c$]",
        ylim=(0, 220),
    )
    locus.text(
        0.5,
        -0.19,
        "Sampled decays overlaid on the two-body prediction",
        transform=locus.transAxes,
        ha="center",
        fontsize=8.8,
        color="#475769",
    )

    residual = 1e6 * (values(pairs, "m_g4_mev") - values(pairs, "m_gen_mev"))
    width = max(np.max(np.abs(residual)) * 1.03, 1)
    mass.hist(residual, bins=np.linspace(-width, width, 55), color=BLUE)
    mass.axvline(0, color=INK, ls="--", lw=1)
    panel(mass, "Dipion invariant-mass closure")
    mass.set(
        xlabel=r"$M_{\pi\pi}^{\rm G4}-M_{\pi\pi}^{\rm HepMC}$  [eV/$c^2$]",
        ylabel="Decays / bin",
    )
    mass.yaxis.set_major_locator(MaxNLocator(3, integer=True))
    mass.text(
        0.02,
        0.91,
        f"Mean mass\n{values(pairs, 'm_g4_mev').mean():.5f} MeV/$c^2$",
        transform=mass.transAxes,
        va="top",
        fontsize=9,
    )

    edges = np.linspace(-1, 1, 21)
    counts, _ = np.histogram(values(pairs, "g4_costh"), bins=edges)
    generator, _ = np.histogram(values(pairs, "gen_costh"), bins=edges)
    expectation = len(pairs) / (len(edges) - 1)
    angle.stairs(
        generator, edges, color="#8da2c1", lw=1.15, label="HepMC", baseline=None
    )
    angle.errorbar(
        (edges[:-1] + edges[1:]) / 2,
        counts,
        np.sqrt(counts),
        fmt="o",
        ms=4,
        color=TEAL,
        capsize=2,
        label="Geant4",
    )
    angle.axhline(expectation, color=INK, ls="--", lw=1.3, label="Isotropic")
    panel(angle, "Positive-pion direction in the kaon rest frame")
    angle.set(
        xlabel=r"$\cos\,\theta^*$ relative to the kaon flight direction",
        ylabel="Decays / bin",
        xlim=(-1, 1),
    )
    angle.yaxis.set_major_locator(MaxNLocator(3, integer=True))
    angle.legend(loc="upper left", ncol=3, frameon=False, fontsize=8.5)
    low, high = angle.get_ylim()
    angle.set_ylim(low, high + (high - low) * 0.32)
    return fig


def flight(gun, odd):
    pairs = [row for row in gun if flag(row, "pair")]
    odd_pairs = [row for row in odd if flag(row, "pair")]
    if len(pairs) != len(gun):
        raise ValueError(
            "An unbiased gun lifetime plot requires all parents to have decayed"
        )
    fig = figure(
        "Kaon flight before decay",
        "Isolated kaons test the lifetime distribution. Hadronic Z events test the prescribed vertex in the OpenDataDetector.",
    )
    gs = fig.add_gridspec(
        1, 2, left=0.08, right=0.97, top=0.75, bottom=0.18, wspace=0.34
    )
    life, parity = [fig.add_subplot(gs[i]) for i in range(2)]
    proper = np.sort(values(pairs, "g4_L_mm") / values(pairs, "beta_gamma"))
    observed = 1 - np.arange(len(proper)) / len(proper)
    time = np.linspace(0, 5 * CTAU, 400)
    theory = np.exp(-time / CTAU)
    std = np.sqrt(theory * (1 - theory) / len(proper))
    life.fill_between(
        time,
        theory - 1.96 * std,
        theory + 1.96 * std,
        color="#d1dfed",
        label="95% pointwise band",
    )
    life.step(
        np.r_[0, proper[:-1]],
        observed,
        where="post",
        color=TEAL,
        lw=2.1,
        label="Geant4",
    )
    life.plot(time, theory, color=INK, lw=1.5, ls="--", label="PDG exponential")
    panel(life, "Survival versus proper flight length")
    life.set(
        xlabel=r"Proper flight length  $L\,m_K/|\mathbf{p}_K|$  [mm]",
        ylabel="Fraction surviving",
        yscale="log",
        xlim=(0, time[-1]),
        ylim=(0.004, 1.15),
    )
    life.legend(loc="upper right", frameon=False, fontsize=9)
    life.text(
        0.04,
        0.07,
        f"Mean {proper.mean():.2f} ± {proper.std(ddof=1)/np.sqrt(len(proper)):.2f} mm\nPDG  {CTAU:.3f} mm",
        transform=life.transAxes,
        fontsize=10,
        va="bottom",
    )

    dp = values(odd_pairs, "parent_dp_mev")
    dv = values(odd_pairs, "vertex_delta_mm")
    nuclear = values(odd_pairs, "direct_nuclear_children") > 0
    parity.scatter(
        dp[~nuclear],
        dv[~nuclear],
        s=16,
        color=BLUE,
        alpha=0.48,
        edgecolors="none",
        rasterized=True,
        label="Pion pair",
    )
    parity.scatter(
        dp[nuclear],
        dv[nuclear],
        s=40,
        marker="D",
        color=ORANGE,
        edgecolors="white",
        linewidth=0.45,
        zorder=3,
        label="Pair + nuclear products",
    )
    panel(parity, "Paired decays in the OpenDataDetector")
    parity.set(
        xlabel=r"$|\mathbf{p}_{\pi\pi}^{\rm G4}-\mathbf{p}_K^{\rm HepMC}|$  [MeV/$c$]",
        ylabel=r"3D decay-vertex shift from HepMC  [mm]",
        xscale="log",
        yscale="log",
    )
    parity.legend(loc="lower right", frameon=False, fontsize=9)
    parity.text(
        0.035,
        0.96,
        f"{len(odd_pairs):,} truth-matched pairs",
        transform=parity.transAxes,
        va="top",
        fontsize=10,
    )
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gun", type=Path, required=True)
    parser.add_argument("--odd", type=Path, required=True)
    parser.add_argument("--fits", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "text.color": INK,
            "axes.labelcolor": INK,
            "axes.edgecolor": "#84919c",
            "xtick.color": INK,
            "ytick.color": INK,
            "pdf.fonttype": 42,
        }
    )
    gun, odd, fits = (read_candidates(path) for path in (args.gun, args.odd, args.fits))
    with PdfPages(args.output / "ks-decay-validation.pdf") as pdf:
        for name, fig in [
            ("ks-reconstruction", reconstruction(fits, odd)),
            ("ks-kinematics", kinematics(gun)),
            ("ks-flight", flight(gun, odd)),
        ]:
            fig.savefig(
                args.output / f"{name}.png",
                dpi=155,
                bbox_inches="tight",
                pad_inches=0.18,
            )
            pdf.savefig(fig, bbox_inches="tight", pad_inches=0.18)
            plt.close(fig)


if __name__ == "__main__":
    main()
