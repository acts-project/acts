#!/usr/bin/env python3
"""Minimal truth-tracking demonstrator using ColliderML data.

Reads particles and tracker hits from a ColliderML parquet dataset,
converts them to ACTS EDM (SimParticleContainer, SimHitContainer,
MeasurementContainer), selects particles with pT > 1 GeV and >= 5 hits,
runs truth-smeared seeding, and fits tracks with the Kalman filter.

Usage
-----
./run_in_env.sh python3 Examples/Scripts/Python/colliderml_truth_tracking.py \\
    --input  colliderml-sample/CERN__ColliderML-Release-1 \\
    --digi-config build/_deps/odd-src/config/odd-digi-smearing-config.json \\
    --output colliderml_output
"""

import argparse
import pathlib

import acts
import acts.examples
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory
from acts.examples.simulation import (
    addDigiParticleSelection,
    ParticleSelectorConfig,
)
from acts.examples.reconstruction import (
    SeedingAlgorithm,
    addSeeding,
    addKalmanTracks,
)

u = acts.UnitConstants


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        "-i",
        type=pathlib.Path,
        required=True,
        help="ColliderML sample root directory "
        "(contains ttbar_pu200_{particles,tracker_hits}/)",
    )
    parser.add_argument(
        "--geo-map",
        type=pathlib.Path,
        default=None,
        help="ColliderML → ACTS geometry ID map CSV "
        "(from generate_colliderml_geo_map.py). "
        "Omit to use direct volume/layer/surface passthrough.",
    )
    parser.add_argument(
        "--digi-config",
        type=pathlib.Path,
        required=True,
        help="ACTS smearing digitisation config JSON "
        "(e.g. odd-digi-smearing-config.json)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "colliderml_output",
        help="Output directory (default: colliderml_output)",
    )
    parser.add_argument(
        "--events",
        "-n",
        type=int,
        default=10,
        help="Number of events to process (default: 10)",
    )
    parser.add_argument(
        "--odd-dir",
        type=pathlib.Path,
        default=None,
        help="ODD XML directory (default: auto-detect)",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel threads (default: 1)",
    )
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Detector + field
    # ------------------------------------------------------------------
    odd_dir = args.odd_dir or getOpenDataDetectorDirectory()
    odd = getOpenDataDetector(odd_dir=odd_dir)
    tgeo = odd.trackingGeometry()
    decorators = odd.contextDecorators()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    # ------------------------------------------------------------------
    # Sequencer
    # ------------------------------------------------------------------
    s = acts.examples.Sequencer(
        events=args.events,
        numThreads=args.jobs,
        logLevel=acts.logging.INFO,
        outputDir=str(args.output),
        # boost::histogram log-axis fills produce benign FPEs in ROOT writers
        failOnUnmaskedFpe=False,
    )

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    # ------------------------------------------------------------------
    # Imports
    # ------------------------------------------------------------------
    from acts.arrow import particleSchema, simHitSchema
    from acts.examples.arrow import (
        loadColliderMLGeoIdMap,
        ColliderMLInputConverter,
        ParquetReader,
    )
    from acts.examples.json import readDigiConfigFromJson
    from acts.examples.root import (
        RootTrackStatesWriter,
        RootTrackSummaryWriter,
        RootTrackFitterPerformanceWriter,
    )
    from acts.examples import PythonTrackFinderPerformanceWriter

    geo_id_map = loadColliderMLGeoIdMap(args.geo_map) if args.geo_map else {}
    digi_config = readDigiConfigFromJson(str(args.digi_config))

    # ------------------------------------------------------------------
    # ParquetReader: place both tables on the whiteboard per event
    # ------------------------------------------------------------------
    particles_dir = (
        args.input / "ttbar_pu200_particles" / "data" / "ttbar_pu200_particles"
    ).resolve()
    hits_dir = (
        args.input / "ttbar_pu200_tracker_hits" / "data" / "ttbar_pu200_tracker_hits"
    ).resolve()

    s.addReader(
        ParquetReader(
            level=acts.logging.INFO,
            inputDir=str(particles_dir.parent),  # unused but required
            collections={
                "cml_particles": str(particles_dir),
                "cml_hits": str(hits_dir),
            },
            expectedSchemas={
                "cml_particles": particleSchema(),
                "cml_hits": simHitSchema(),
            },
        )
    )

    # ------------------------------------------------------------------
    # ColliderMLInputConverter: particles + simhits + measurements
    # ------------------------------------------------------------------
    s.addAlgorithm(
        ColliderMLInputConverter(
            level=acts.logging.INFO,
            inputParticlesTable="cml_particles",
            inputHitsTable="cml_hits",
            outputParticles="particles",
            outputSimHits="simhits",
            outputMeasurements="measurements",
            outputMeasurementSubset="measurement_subset",
            outputMeasSimHitsMap="measurement_simhits_map",
            outputMeasParticlesMap="measurement_particles_map",
            outputParticleMeasurementsMap="particle_measurements_map",
            trackingGeometry=tgeo,
            digiConfig=digi_config,
            geoIdMap=geo_id_map,
        )
    )

    # ------------------------------------------------------------------
    # Particle selection: pT > 1 GeV, >= 5 measurements, charged only
    # addDigiParticleSelection reads "particles_simulated_selected"
    # ------------------------------------------------------------------
    s.addWhiteboardAlias("particles_simulated_selected", "particles")

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(1.0 * u.GeV, None),
            measurements=(5, None),
            removeNeutral=True,
        ),
    )
    # addDigiParticleSelection aliases output to "particles_selected"

    # ------------------------------------------------------------------
    # Truth-estimated seeding: space-point triplets matched to truth
    # particles.  Initial track parameters are estimated from the triplet
    # geometry (on measurement surfaces), not from the production vertex,
    # which avoids navigator failures for off-origin pileup vertices.
    # ------------------------------------------------------------------
    odd_seeding_cfg = (
        pathlib.Path(__file__).resolve().parent.parent.parent.parent
        / "build/_deps/odd-src/config/odd-seeding-config.json"
    )
    addSeeding(
        s,
        trackingGeometry=tgeo,
        field=field,
        rnd=rnd,
        seedingAlgorithm=SeedingAlgorithm.TruthEstimated,
        selectedParticles="particles_selected",
        geoSelectionConfigFile=odd_seeding_cfg,
        initialSigmas=[
            1 * u.mm,
            1 * u.mm,
            1 * u.degree,
            1 * u.degree,
            0 / u.GeV,
            1 * u.ns,
        ],
        initialSigmaQoverPt=0.1 / u.GeV,
        initialSigmaPtRel=0.1,
        initialVarInflation=[1e0, 1e0, 1e0, 1e0, 1e0, 1e0],
        logLevel=acts.logging.INFO,
    )

    # ------------------------------------------------------------------
    # Kalman filter track fitting
    # ------------------------------------------------------------------
    addKalmanTracks(
        s,
        trackingGeometry=tgeo,
        field=field,
        logLevel=acts.logging.INFO,
    )

    # ------------------------------------------------------------------
    # Track selection
    # ------------------------------------------------------------------
    s.addAlgorithm(
        acts.examples.TrackSelectorAlgorithm(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTracks="selected_tracks",
            selectorConfig=acts.TrackSelector.Config(
                minMeasurements=5,
            ),
        )
    )
    s.addWhiteboardAlias("tracks", "selected_tracks")

    # ------------------------------------------------------------------
    # ROOT output
    # ------------------------------------------------------------------
    s.addWriter(
        RootTrackStatesWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            inputSimHits="simhits",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(args.output / "trackstates_kf.root"),
        )
    )
    s.addWriter(
        RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(args.output / "tracksummary_kf.root"),
        )
    )
    s.addWriter(
        RootTrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(args.output / "performance_kf.root"),
        )
    )

    # ------------------------------------------------------------------
    # Python performance writer (ROOT-free, exposes histograms() after run)
    # ------------------------------------------------------------------
    perf_cfg = PythonTrackFinderPerformanceWriter.Config()
    perf_cfg.inputTracks = "tracks"
    perf_cfg.inputParticles = "particles_selected"
    perf_cfg.inputTrackParticleMatching = "track_particle_matching"
    perf_cfg.inputParticleTrackMatching = "particle_track_matching"
    perf_cfg.inputParticleMeasurementsMap = "particle_measurements_map"
    perf_writer = PythonTrackFinderPerformanceWriter(perf_cfg, acts.logging.INFO)
    s.addWriter(perf_writer)

    s.run()

    _make_plots(perf_writer.histograms(), args.output)


def _make_plots(hists, output_dir):
    """Save four-page performance PDF to output_dir/performance_plots.pdf."""
    import numpy as np
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "figure.dpi": 120,
        }
    )

    # ------------------------------------------------------------------
    # Low-level helpers
    # ------------------------------------------------------------------

    def _eff_arrays(eff1):
        edges = np.asarray(eff1.total.axis(0).edges)
        centers = 0.5 * (edges[:-1] + edges[1:])
        acc = np.asarray(eff1.accepted.values())
        tot = np.asarray(eff1.total.values())
        mask = tot > 0
        eff = np.where(mask, acc / tot, np.nan)
        # Wilson / binomial error
        err = np.where(mask, np.sqrt(acc * (tot - acc) / tot**3), np.nan)
        return centers, eff, err, eff1.total.axis(0).label

    def _prof_arrays(prof1):
        bh = prof1.histogram
        edges = np.asarray(bh.axis(0).edges)
        centers = 0.5 * (edges[:-1] + edges[1:])
        counts = np.asarray(bh.counts())
        means = np.asarray(bh.means())
        sum_dq = np.asarray(bh.sum_of_deltas_squared())
        with np.errstate(invalid="ignore", divide="ignore"):
            var = np.where(counts > 1, sum_dq / (counts - 1), 0.0)
            err = np.where(counts > 0, np.sqrt(var / np.maximum(counts, 1)), np.nan)
        means = np.where(counts > 0, means, np.nan)
        return centers, means, err, bh.axis(0).label

    def _plot_eff(ax, key, label=None, color=None, xscale="linear"):
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
        ax.grid(True, alpha=0.3)

    def _plot_prof(ax, key, ylabel=None, label=None, color=None, ymin=None):
        if key not in hists:
            ax.text(0.5, 0.5, f"(no data: {key})", transform=ax.transAxes, ha="center")
            return
        x, y, ye, xlabel = _prof_arrays(hists[key])
        kw = dict(fmt="o-", markersize=3, capsize=2, lw=1, elinewidth=0.8)
        if label is not None:
            kw["label"] = label
        if color is not None:
            kw["color"] = color
        ax.errorbar(x, y, yerr=ye, **kw)
        ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        if ymin is not None:
            ax.set_ylim(bottom=ymin)
        ax.grid(True, alpha=0.3)

    output_path = output_dir / "performance_plots.pdf"
    with PdfPages(output_path) as pdf:
        # ------------------------------------------------------------------
        # Page 1 — Core tracking performance (efficiency, fake, duplication)
        # ------------------------------------------------------------------
        fig, axes = plt.subplots(2, 2, figsize=(14, 9))
        fig.suptitle(
            "Tracking Performance — Kalman Filter, ColliderML ttbar PU200",
            fontsize=14,
        )

        _plot_eff(axes[0, 0], "trackeff_vs_eta")
        axes[0, 0].set_ylabel("Efficiency")
        axes[0, 0].set_title("Tracking efficiency vs η")

        _plot_eff(axes[0, 1], "trackeff_vs_pT", xscale="log")
        axes[0, 1].set_ylabel("Efficiency")
        axes[0, 1].set_title("Tracking efficiency vs pT")

        _plot_eff(axes[1, 0], "fakeRatio_vs_eta")
        axes[1, 0].set_ylabel("Fake ratio")
        axes[1, 0].set_title("Fake ratio vs η")

        _plot_eff(axes[1, 1], "duplicationRatio_vs_eta")
        axes[1, 1].set_ylabel("Duplication ratio")
        axes[1, 1].set_title("Duplication ratio vs η")

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # ------------------------------------------------------------------
        # Page 2 — Hit content (measurements, holes, outliers, shared hits)
        # ------------------------------------------------------------------
        fig, axes = plt.subplots(2, 3, figsize=(16, 9))
        fig.suptitle("Hit Content per Track", fontsize=14)

        hit_panels = [
            ("nMeasurements_vs_eta", "nMeasurements"),
            ("nHoles_vs_eta", "nHoles"),
            ("nOutliers_vs_eta", "nOutliers"),
            ("nMeasurements_vs_pT", "nMeasurements"),
            ("nHoles_vs_pT", "nHoles"),
            ("nSharedHits_vs_eta", "nSharedHits"),
        ]
        for ax, (key, ylabel) in zip(axes.flat, hit_panels):
            _plot_prof(ax, key, ylabel=ylabel, ymin=0)
            ax.set_title(key.replace("_", " "))

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # ------------------------------------------------------------------
        # Page 3 — Track quality (completeness and purity)
        # ------------------------------------------------------------------
        fig, axes = plt.subplots(2, 2, figsize=(14, 9))
        fig.suptitle("Track Quality", fontsize=14)

        quality_panels = [
            ("completeness_vs_eta", "Completeness"),
            ("purity_vs_eta", "Purity"),
            ("completeness_vs_pT", "Completeness"),
            ("purity_vs_pT", "Purity"),
        ]
        for ax, (key, ylabel) in zip(axes.flat, quality_panels):
            _plot_prof(ax, key, ylabel=ylabel, ymin=0)
            ax.set_title(key.replace("_", " "))
            ax.axhline(1.0, color="gray", ls="--", lw=0.8)
            ax.set_ylim(0, 1.05)

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # ------------------------------------------------------------------
        # Page 4 — Differential efficiency (pT slices × η slices)
        # ------------------------------------------------------------------
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle("Differential Tracking Efficiency", fontsize=14)

        pt_labels = ["~1 GeV", "~10 GeV", "~100 GeV"]
        for i, (color, label) in enumerate(zip(["C0", "C1", "C2"], pt_labels)):
            _plot_eff(
                axes[0],
                f"trackeff_vs_eta_ptRange_{i}",
                label=f"pT ∈ {label}",
                color=color,
            )
        axes[0].set_ylabel("Efficiency")
        axes[0].set_title("Efficiency vs η in pT slices")
        axes[0].legend(fontsize=10)

        eta_labels = ["|η| < 0.2", "|η| < 0.8", "1 < |η| < 2", "2 < |η| < 3"]
        for i, (color, label) in enumerate(zip(["C0", "C1", "C2", "C3"], eta_labels)):
            _plot_eff(
                axes[1],
                f"trackeff_vs_pT_absEtaRange_{i}",
                label=label,
                color=color,
                xscale="log",
            )
        axes[1].set_ylabel("Efficiency")
        axes[1].set_title("Efficiency vs pT in |η| slices")
        axes[1].legend(fontsize=10)

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    print(f"Saved performance plots → {output_path}")


if __name__ == "__main__":
    main()
