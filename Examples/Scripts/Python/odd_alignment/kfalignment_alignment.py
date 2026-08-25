"""ODD KF alignment workflow — alignment calibration stage."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import acts
import acts.examples
import acts.examples.alignment
from acts.examples.alignment import AlignmentDecorator

from kfalignment_config import (
    DEFAULT_INITIAL_VAR_INFLATION,
    DEFAULT_MISALIGNMENT_CONFIG,
    DEFAULT_NUM_EVENTS,
    DEFAULT_NUM_THREADS,
    MisalignmentConfig,
)
from kfalignment_transform_io import (
    load_misalignment_for_alignment,
    read_alignment_delta_res,
)

u = acts.UnitConstants


def runAlignment(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    inputDir: Path,
    outputDir: Path,
    numEvents: int = DEFAULT_NUM_EVENTS,
    s: acts.examples.Sequencer = None,
    misalignmentDir: Path = None,
    misalignmentConfig: Optional[MisalignmentConfig] = None,
):
    """Align misaligned geometry; write ``aligned_result.txt`` (+ index map)."""
    from acts.examples.simulation import (
        ParticleSelectorConfig,
        addDigiParticleSelection,
    )
    from acts.examples.reconstruction import (
        addSeeding,
        SeedingAlgorithm,
    )

    s = s or acts.examples.Sequencer(
        events=numEvents, numThreads=DEFAULT_NUM_THREADS, logLevel=acts.logging.INFO
    )

    rnd = acts.examples.RandomNumbers(seed=42)
    inputDir = Path(inputDir)
    outputDir = Path(outputDir)
    outputDir.mkdir(parents=True, exist_ok=True)
    misalignmentDir = Path(misalignmentDir) if misalignmentDir else outputDir

    logger = acts.getDefaultLogger("Alignment", acts.logging.INFO)

    cfg = (
        misalignmentConfig
        if misalignmentConfig is not None
        else DEFAULT_MISALIGNMENT_CONFIG
    )

    applied_path = misalignmentDir / "misalignment_applied.txt"
    index_path = misalignmentDir / "misalignment_index_map.txt"
    if not applied_path.exists() or not index_path.exists():
        raise FileNotFoundError(
            "Alignment requires misalignment files:\n"
            f"  {applied_path}\n  {index_path}\n"
            "Run `python full_chain_odd_KFalignment.py misalignment` first "
            "(or pass --misalignment-dir)."
        )

    geoIdMap, alignment_placements, surface_to_geoid = load_misalignment_for_alignment(
        trackingGeometry, applied_path, index_path
    )
    tx, ty, tz = cfg.tx, cfg.ty, cfg.tz
    rx, ry, rz = cfg.rx, cfg.ry, cfg.rz

    print(
        f"Alignment: loaded misalignment ({len(geoIdMap)} elements) from "
        f"{applied_path}; fit DoF tx={tx} ty={ty} tz={tz} rx={rx} ry={ry} rz={rz}"
    )

    mutableStore = acts.examples.alignment.MutableGeoIdAlignmentStore(geoIdMap)

    deco_cfg = AlignmentDecorator.Config()
    deco_cfg.iovStores = [((0, 1_000_000), mutableStore)]
    deco_cfg.garbageCollection = False
    alignDeco = AlignmentDecorator(deco_cfg, acts.logging.INFO)
    s.addContextDecorator(alignDeco)

    logger.info("Reading simulation results from CSV files in {}", str(inputDir))
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

    Center0 = 1 << 0
    Center1 = 1 << 1
    Center2 = 1 << 2
    Rotation0 = 1 << 3
    Rotation1 = 1 << 4
    Rotation2 = 1 << 5

    alignment_dof = 0
    if tx:
        alignment_dof |= Center0
    if ty:
        alignment_dof |= Center1
    if tz:
        alignment_dof |= Center2
    if rx:
        alignment_dof |= Rotation0
    if ry:
        alignment_dof |= Rotation1
    if rz:
        alignment_dof |= Rotation2

    det_elements = alignment_placements
    if len(det_elements) != len(geoIdMap):
        print(
            f"WARNING: misaligned ({len(geoIdMap)}) vs alignment ({len(det_elements)}) count mismatch"
        )

    aal_cfg = acts.examples.alignment.AlignmentAlgorithmConfig()
    aal_cfg.inputMeasurements = "measurements"
    aal_cfg.inputProtoTracks = "truth_particle_tracks"
    aal_cfg.inputInitialTrackParameters = "estimatedparameters"
    aal_cfg.outputAlignmentParameters = "alignmentParameters"
    aligned_result_file = outputDir / "aligned_result.txt"
    aligned_index_file = outputDir / "aligned_result_index_map.txt"
    aal_cfg.outputAlignmentFile = str(aligned_result_file)
    aal_cfg.outputAlignmentIndexFile = str(aligned_index_file)
    aal_cfg.trackingGeometry = trackingGeometry
    aal_cfg.chi2ONdfCutOff = 0.1
    aal_cfg.deltaChi2ONdfCutOff = (5, 0.0001)

    aal_cfg.alignedTransformUpdater = (
        acts.examples.alignment.makeAlignedTransformUpdater(mutableStore)
    )

    aal_cfg.alignedDetElements = det_elements
    aal_cfg.align = acts.examples.alignment.makeAlignmentFunction(
        trackingGeometry, field, acts.logging.INFO
    )
    aal_cfg.maxNumIterations = 100
    aal_cfg.maxNumTracks = 100000
    aal_cfg.iterationState = {
        i: int(alignment_dof) for i in range(aal_cfg.maxNumIterations)
    }

    alignment_algo = acts.examples.alignment.AlignmentAlgorithm(
        aal_cfg, acts.logging.INFO
    )
    s.addAlgorithm(alignment_algo)

    s.run()

    aligned_deltas = (
        read_alignment_delta_res(aligned_result_file)
        if aligned_result_file.exists()
        else {}
    )
    n_nonzero = sum(
        1 for deltas in aligned_deltas.values() if any(abs(v) > 0.0 for v in deltas)
    )
    if n_nonzero == 0 and surface_to_geoid:
        print(
            "WARNING: aligned_result.txt has no non-zero deltas — "
            "aligned reconstruction will match misaligned. Check alignment fit."
        )
    print(
        f"Alignment done: {len(surface_to_geoid)} surfaces -> "
        f"{aligned_result_file.name}, {aligned_index_file.name} "
        f"({n_nonzero} with non-zero deltas; written by C++ finalize)"
    )

    return s
