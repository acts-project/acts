#!/usr/bin/env python3

from pathlib import Path
from typing import Optional

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

_srcdir = Path(__file__).resolve().parent.parent.parent.parent


def runColliderMLTruthTracking(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    outputDir: Path,
    inputDir: Path,
    geoIdMapPath: Optional[Path] = None,
    geoIdMapSourcePrefix: str = "gen1",
    geoIdMapTargetPrefix: str = "gen3",
    decorators=[],
    events: int = 10,
    numThreads: int = 1,
    sample: str = "ttbar_pu200",
    s: Optional[acts.examples.Sequencer] = None,
):
    """Set up a ColliderML truth-tracking sequencer and return it with the performance writer.

    Returns
    -------
    (Sequencer, PythonTrackFinderPerformanceWriter)
        Call s.run() on the sequencer, then access perf_writer.histograms().
    """
    from acts.examples.arrow import (
        ColliderMLRelease1InputConverter,
        ParquetReader,
    )
    from acts.examples.root import (
        RootTrackStatesWriter,
        RootTrackSummaryWriter,
        RootTrackFitterPerformanceWriter,
    )
    from acts.examples import PythonTrackFinderPerformanceWriter

    outputDir = Path(outputDir)
    outputDir.mkdir(parents=True, exist_ok=True)

    s = s or acts.examples.Sequencer(
        events=events,
        numThreads=numThreads,
        logLevel=acts.logging.INFO,
        outputDir=str(outputDir),
        failOnUnmaskedFpe=False,
    )

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    particles_dir = (
        inputDir / f"{sample}_particles" / "data" / f"{sample}_particles"
    ).resolve()
    hits_dir = (
        inputDir / f"{sample}_tracker_hits" / "data" / f"{sample}_tracker_hits"
    ).resolve()

    s.addReader(
        ParquetReader(
            level=acts.logging.INFO,
            inputDir=str(particles_dir.parent),
            collections={
                "cml_particles": str(particles_dir),
                "cml_hits": str(hits_dir),
            },
            expectedSchemas={
                "cml_particles": ColliderMLRelease1InputConverter.particleSchema(),
                "cml_hits": ColliderMLRelease1InputConverter.hitSchema(),
            },
        )
    )

    converter_kwargs = dict(
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
        trackingGeometry=trackingGeometry,
    )
    if geoIdMapPath is not None:
        converter_kwargs["geoIdMapPath"] = geoIdMapPath
        converter_kwargs["geoIdMapSourcePrefix"] = geoIdMapSourcePrefix
        converter_kwargs["geoIdMapTargetPrefix"] = geoIdMapTargetPrefix

    s.addAlgorithm(ColliderMLRelease1InputConverter(**converter_kwargs))

    s.addWhiteboardAlias("particles_simulated_selected", "particles")
    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(1.0 * u.GeV, None),
            measurements=(5, None),
            removeNeutral=True,
        ),
    )

    addSeeding(
        s,
        trackingGeometry=trackingGeometry,
        field=field,
        rnd=rnd,
        seedingAlgorithm=SeedingAlgorithm.TruthEstimated,
        selectedParticles="particles_selected",
        geoSelectionConfigFile=_srcdir / "Examples/Configs/odd-seeding-config.json",
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

    addKalmanTracks(
        s,
        trackingGeometry=trackingGeometry,
        field=field,
        logLevel=acts.logging.INFO,
    )

    s.addAlgorithm(
        acts.examples.TrackSelectorAlgorithm(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTracks="selected_tracks",
            selectorConfig=acts.TrackSelector.Config(
                minMeasurements=7,
            ),
        )
    )
    s.addWhiteboardAlias("tracks", "selected_tracks")

    s.addWriter(
        RootTrackStatesWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            inputSimHits="simhits",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDir / "trackstates_kf.root"),
        )
    )
    s.addWriter(
        RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "tracksummary_kf.root"),
        )
    )
    s.addWriter(
        RootTrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "performance_kf.root"),
        )
    )

    # Both writers use truth_seeded_particles as the denominator so efficiency
    # is evaluated only on particles that the seeding algorithm actually found.
    # This cleanly separates seeding coverage (determined upstream by the
    # geoSelection config) from KF quality.

    # Proto-track level: how many seeded particles survive as seed-tracks.
    perf_proto_cfg = PythonTrackFinderPerformanceWriter.Config()
    perf_proto_cfg.inputTracks = "seed-tracks"
    perf_proto_cfg.inputParticles = "truth_seeded_particles"
    perf_proto_cfg.inputTrackParticleMatching = "seed_particle_matching"
    perf_proto_cfg.inputParticleTrackMatching = "particle_seed_matching"
    perf_proto_cfg.inputParticleMeasurementsMap = "particle_measurements_map"
    perf_proto_writer = PythonTrackFinderPerformanceWriter(
        perf_proto_cfg, acts.logging.INFO
    )
    s.addWriter(perf_proto_writer)

    # KF track level: how many seeded particles survive KF fitting + selection.
    perf_cfg = PythonTrackFinderPerformanceWriter.Config()
    perf_cfg.inputTracks = "tracks"
    perf_cfg.inputParticles = "truth_seeded_particles"
    perf_cfg.inputTrackParticleMatching = "track_particle_matching"
    perf_cfg.inputParticleTrackMatching = "particle_track_matching"
    perf_cfg.inputParticleMeasurementsMap = "particle_measurements_map"
    perf_writer = PythonTrackFinderPerformanceWriter(perf_cfg, acts.logging.INFO)
    s.addWriter(perf_writer)

    return s, perf_proto_writer, perf_writer


def _serialize_hists(hists):
    """Convert PythonTrackFinderPerformanceWriter histograms to a picklable dict.

    The writer exposes three C++ histogram wrapper types, none of which are
    picklable. We extract numpy arrays keyed by type tag:
      "efficiency"  — Efficiency1/2  (has .accepted / .total BoostHistogram)
      "profile"     — ProfileHistogram1 (.histogram is BoostProfileHistogram)
      "histogram"   — Histogram1/2/3   (.histogram is plain BoostHistogram)
    """
    import numpy as np

    out = {}
    for key, h in hists.items():
        if hasattr(h, "accepted"):
            # Efficiency1 / Efficiency2
            out[key] = {
                "type": "efficiency",
                "edges": np.asarray(h.total.axis(0).edges),
                "accepted": np.asarray(h.accepted.values()),
                "total": np.asarray(h.total.values()),
                "label": h.total.axis(0).label,
            }
        elif hasattr(h, "histogram"):
            bh = h.histogram
            if hasattr(bh, "counts"):
                # ProfileHistogram1 — BoostProfileHistogram
                out[key] = {
                    "type": "profile",
                    "edges": np.asarray(bh.axis(0).edges),
                    "counts": np.asarray(bh.counts()),
                    "means": np.asarray(bh.means()),
                    "sum_of_deltas_squared": np.asarray(bh.sum_of_deltas_squared()),
                    "label": bh.axis(0).label,
                }
            else:
                # Histogram1/2/3 — plain BoostHistogram (from FakePlotTool)
                out[key] = {
                    "type": "histogram",
                    "edges": np.asarray(bh.axis(0).edges),
                    "values": np.asarray(bh.values()),
                    "label": bh.axis(0).label,
                }
    return out


if __name__ == "__main__":
    import argparse
    import pickle

    parser = argparse.ArgumentParser(
        description="ColliderML truth-tracking Kalman filter demo on ttbar PU200."
    )
    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        required=True,
        help="ColliderML sample root (contains ttbar_pu200_{particles,tracker_hits}/)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path.cwd() / "colliderml_output",
        help="Output directory (default: colliderml_output)",
    )
    parser.add_argument(
        "--events",
        "-n",
        type=int,
        default=10,
        help="Number of events (default: 10)",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel threads (default: 1)",
    )
    args = parser.parse_args()

    detector = getOpenDataDetector()
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    s, perf_proto_writer, perf_writer = runColliderMLTruthTracking(
        trackingGeometry=trackingGeometry,
        field=field,
        outputDir=args.output,
        inputDir=args.input,
        decorators=decorators,
        events=args.events,
        numThreads=args.jobs,
    )
    s.run()

    hist_path = args.output / "histograms.pkl"
    with open(hist_path, "wb") as f:
        pickle.dump(_serialize_hists(perf_writer.histograms()), f)
    print(f"Saved histograms → {hist_path}")

    proto_hist_path = args.output / "histograms_proto.pkl"
    with open(proto_hist_path, "wb") as f:
        pickle.dump(_serialize_hists(perf_proto_writer.histograms()), f)
    print(f"Saved proto-track histograms → {proto_hist_path}")
