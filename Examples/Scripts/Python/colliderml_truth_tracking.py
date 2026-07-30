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
    particlesDir: Path,
    hitsDir: Path,
    tracksDir: Optional[Path] = None,
    geoIdMapPath: Optional[Path] = None,
    geoIdMapSourcePrefix: str = "gen1",
    geoIdMapTargetPrefix: str = "gen3",
    decorators=[],
    events: Optional[int] = None,
    numThreads: int = 1,
    s: Optional[acts.examples.Sequencer] = None,
):
    """Set up a ColliderML truth-tracking sequencer and return it with the performance writer.

    If `tracksDir` is given, the dataset's own published tracks are also read
    and converted to a `ConstTrackContainer` on the whiteboard under
    "colliderml_tracks" -- distinct from "tracks", which holds this script's
    own truth-seeded KF tracks.

    Returns
    -------
    (Sequencer, PythonTrackFinderPerformanceWriter)
        Call s.run() on the sequencer, then access perf_writer.histograms().
    """
    from acts.examples.arrow import (
        ColliderMLRelease1InputConverter,
        ParquetReader,
    )

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

    readerCollections = {
        "cml_particles": str(particlesDir),
        "cml_hits": str(hitsDir),
    }
    readerSchemas = {
        "cml_particles": ColliderMLRelease1InputConverter.particleSchema(),
        "cml_hits": ColliderMLRelease1InputConverter.hitSchema(),
    }
    if tracksDir is not None:
        readerCollections["cml_tracks"] = str(tracksDir)
        readerSchemas["cml_tracks"] = ColliderMLRelease1InputConverter.tracksSchema()

    s.addReader(
        ParquetReader(
            level=acts.logging.INFO,
            collections=readerCollections,
            expectedSchemas=readerSchemas,
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
    if tracksDir is not None:
        converter_kwargs["inputTracksTable"] = "cml_tracks"
        converter_kwargs["outputTracks"] = "colliderml_tracks"

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

    return s


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="ColliderML truth-tracking Kalman filter demo on ttbar PU200."
    )
    parser.add_argument(
        "--particlesDir",
        "-p",
        type=Path,
        required=True,
        help="ColliderML particles directory",
    )
    parser.add_argument(
        "--hitsDir",
        "-m",
        type=Path,
        required=True,
        help="ColliderML hits directory",
    )
    parser.add_argument(
        "--tracksDir",
        "-t",
        type=Path,
        default=None,
        help="ColliderML published-tracks directory. Optional: not all samples ship these.",
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
        default=None,
        help="Number of events (default: all available events in the dataset)",
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

    s = runColliderMLTruthTracking(
        trackingGeometry=trackingGeometry,
        field=field,
        outputDir=args.output,
        particlesDir=args.particlesDir,
        hitsDir=args.hitsDir,
        tracksDir=args.tracksDir,
        decorators=decorators,
        events=args.events,
        numThreads=args.jobs,
    )
    s.run()
