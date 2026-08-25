"""ODD KF alignment workflow — Kalman reconstruction stage (nominal / misaligned / aligned)."""

from __future__ import annotations

import sys
from pathlib import Path

import acts
import acts.examples
import acts.examples.alignment
from acts.examples.alignment import AlignmentDecorator

from kfalignment_config import (
    DEFAULT_INITIAL_VAR_INFLATION,
    DEFAULT_NUM_EVENTS,
    DEFAULT_NUM_THREADS,
)
from kfalignment_transform_io import (
    geo_id_map_from_aligned,
    geo_id_map_from_misalignment,
)

u = acts.UnitConstants


def runReconstruction(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    inputDir: Path,
    outputDir: Path,
    numEvents: int = DEFAULT_NUM_EVENTS,
    reverseFilteringMomThreshold=0 * u.GeV,
    reverseFilteringCovarianceScaling=1,
    s: acts.examples.Sequencer = None,
    mode: str = "nominal",
    alignmentDir: Path = None,
):
    """Run Kalman reconstruction for ``mode``: nominal / misaligned / aligned."""
    from acts.examples.simulation import (
        ParticleSelectorConfig,
        addDigiParticleSelection,
    )
    from acts.examples.reconstruction import (
        addSeeding,
        SeedingAlgorithm,
        addKalmanTracks,
    )

    s = s or acts.examples.Sequencer(
        events=numEvents, numThreads=DEFAULT_NUM_THREADS, logLevel=acts.logging.INFO
    )

    rnd = acts.examples.RandomNumbers(seed=42)
    inputDir = Path(inputDir)
    outputDir = Path(outputDir)

    logger = acts.getDefaultLogger("Reconstruction", acts.logging.INFO)

    if alignmentDir is None:
        alignmentDir = Path(__file__).resolve().parent / "alignment_output"

    if mode == "nominal":
        use_alignment_decorator = False
    elif mode in ("misaligned", "aligned"):
        use_alignment_decorator = True
    else:
        raise ValueError(
            f"Unknown mode: {mode}. Must be 'nominal', 'misaligned', or 'aligned'"
        )

    if use_alignment_decorator:
        if mode == "misaligned":
            applied_path = alignmentDir / "misalignment_applied.txt"
            index_path = alignmentDir / "misalignment_index_map.txt"
            if not applied_path.exists() or not index_path.exists():
                print(f"Error: misalignment files not found in {alignmentDir}")
                print(f"  expected: {applied_path.name} and {index_path.name}")
                print("Please run 'misalignment' first to generate the required files.")
                sys.exit(1)
            geoIdMap = geo_id_map_from_misalignment(
                trackingGeometry, applied_path, index_path
            )
            print(
                f"Reconstruction (misaligned): {len(geoIdMap)} elements from "
                f"{applied_path.name} + {index_path.name}"
            )
        else:
            mis_applied = alignmentDir / "misalignment_applied.txt"
            mis_index = alignmentDir / "misalignment_index_map.txt"
            ali_result = alignmentDir / "aligned_result.txt"
            ali_index = alignmentDir / "aligned_result_index_map.txt"
            if not ali_index.exists():
                ali_index = mis_index
            missing = [
                p.name for p in (mis_applied, mis_index, ali_result) if not p.exists()
            ]
            if missing:
                print(
                    f"Error: aligned files incomplete in {alignmentDir}: "
                    f"missing {missing}"
                )
                print("Please run 'misalignment' then 'alignment' to generate them.")
                sys.exit(1)
            geoIdMap = geo_id_map_from_aligned(
                trackingGeometry,
                mis_applied,
                mis_index,
                ali_result,
                ali_index,
            )
            print(
                f"Reconstruction (aligned): {len(geoIdMap)} elements from "
                f"{ali_result.name} + {mis_applied.name}"
            )

        mutableStore = acts.examples.alignment.MutableGeoIdAlignmentStore(geoIdMap)

        cfg = AlignmentDecorator.Config()
        cfg.iovStores = [((0, 1_000_000), mutableStore)]
        cfg.garbageCollection = False
        alignDeco = AlignmentDecorator(cfg, acts.logging.INFO)

        s.addContextDecorator(alignDeco)

    logger.info("Reading simulation CSV inputs from {}", str(inputDir))

    s.addReader(
        acts.examples.CsvSimHitReader(
            level=acts.logging.INFO,
            inputDir=str(inputDir),
            inputStem="hits",
            outputSimHits="simhits",
        )
    )

    s.addReader(
        acts.examples.CsvMeasurementReader(
            level=acts.logging.INFO,
            inputDir=str(inputDir),
            outputMeasurements="measurements",
            outputMeasurementSimHitsMap="measurement_simhits_map",
            outputMeasurementParticlesMap="measurement_particles_map",
            outputParticleMeasurementsMap="particle_measurements_map",
            inputSimHits="simhits",
        )
    )

    s.addReader(
        acts.examples.CsvParticleReader(
            level=acts.logging.INFO,
            inputDir=str(inputDir),
            inputStem="particles_simulated",
            outputParticles="particles_simulated",
        )
    )

    s.addWhiteboardAlias("particles", "particles_simulated")
    s.addWhiteboardAlias("particles_simulated_selected", "particles_simulated")
    s.addWhiteboardAlias("particles_selected", "particles_simulated")

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(2 * u.GeV, None),
            measurements=(5, None),
            removeNeutral=True,
            removeSecondaries=True,
        ),
    )

    addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        inputParticles="particles_selected",
        initialVarInflation=DEFAULT_INITIAL_VAR_INFLATION,
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        particleHypothesis=acts.ParticleHypothesis.muon,
    )

    addKalmanTracks(
        s,
        trackingGeometry,
        field,
        reverseFilteringMomThreshold,
        reverseFilteringCovarianceScaling,
    )

    s.addAlgorithm(
        acts.examples.TrackSelectorAlgorithm(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTracks="selected-tracks",
            selectorConfig=acts.TrackSelector.Config(
                minMeasurements=5,
            ),
        )
    )
    s.addWhiteboardAlias("tracks", "selected-tracks")

    s.addWriter(
        acts.examples.root.RootTrackStatesWriter(
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
        acts.examples.root.RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "tracksummary_kf.root"),
        )
    )

    s.addWriter(
        acts.examples.root.RootTrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "performance_kf.root"),
        )
    )

    return s
