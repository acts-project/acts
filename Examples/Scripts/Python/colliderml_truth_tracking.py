#!/usr/bin/env python3

from pathlib import Path
from typing import Optional

import acts
import acts.examples
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory
from acts.examples.dataset import addColliderML
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
    channel: str = "ttbar",
    pileup: str = "pu0",
    dataDir: Optional[Path] = None,
    particlesDir: Optional[Path] = None,
    hitsDir: Optional[Path] = None,
    tracksDir: Optional[Path] = None,
    readTracks: bool = True,
    geoIdMapPath: Optional[Path] = None,
    geoIdMapSourcePrefix: str = "gen1",
    geoIdMapTargetPrefix: str = "gen3",
    decorators=[],
    events: Optional[int] = None,
    numThreads: int = 1,
    s: Optional[acts.examples.Sequencer] = None,
):
    """Set up a ColliderML truth-tracking sequencer and return it with the performance writer.

    By default, the sample is located automatically from `channel`/`pileup`
    via `acts.examples.dataset.getColliderMLObjectDirectory` (which in turn
    honors `$COLLIDERML_DATA_DIR`, falling back to `~/.cache/colliderml` --
    the same cache used by the `colliderml` Python library). Pass
    `particlesDir`/`hitsDir`/`tracksDir` explicitly to bypass this resolution.

    If `readTracks` is True (default), the dataset's own published tracks are
    also read and converted to a `ConstTrackContainer` on the whiteboard under
    "colliderml_tracks" -- distinct from "tracks", which holds this script's
    own truth-seeded KF tracks.

    Returns
    -------
    (Sequencer, PythonPatternRecognitionPerformanceWriter)
        Call s.run() on the sequencer, then access perf_writer.histograms().
    """
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

    addColliderML(
        s,
        trackingGeometry=trackingGeometry,
        channel=channel,
        pileup=pileup,
        dataDir=dataDir,
        particlesDir=particlesDir,
        hitsDir=hitsDir,
        tracksDir=tracksDir,
        readTracks=readTracks,
        geoIdMapPath=geoIdMapPath,
        geoIdMapSourcePrefix=geoIdMapSourcePrefix,
        geoIdMapTargetPrefix=geoIdMapTargetPrefix,
        logLevel=acts.logging.INFO,
    )

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
        "--channel",
        type=str,
        default="ttbar",
        help="ColliderML physics channel (default: ttbar)",
    )
    parser.add_argument(
        "--pileup",
        type=str,
        default="pu0",
        help="ColliderML pileup token, e.g. pu0 or pu200 (default: pu0)",
    )
    parser.add_argument(
        "--dataDir",
        type=Path,
        default=None,
        help="ColliderML cache directory (default: $COLLIDERML_DATA_DIR or ~/.cache/colliderml)",
    )
    parser.add_argument(
        "--particlesDir",
        "-p",
        type=Path,
        default=None,
        help="Explicit ColliderML particles directory (overrides --channel/--pileup/--dataDir)",
    )
    parser.add_argument(
        "--hitsDir",
        "-m",
        type=Path,
        default=None,
        help="Explicit ColliderML hits directory (overrides --channel/--pileup/--dataDir)",
    )
    parser.add_argument(
        "--tracksDir",
        "-t",
        type=Path,
        default=None,
        help="Explicit ColliderML published-tracks directory (overrides --channel/--pileup/--dataDir)",
    )
    parser.add_argument(
        "--no-tracks",
        dest="readTracks",
        action="store_false",
        help="Don't read the dataset's own published tracks. Not all samples ship these.",
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
        channel=args.channel,
        pileup=args.pileup,
        dataDir=args.dataDir,
        particlesDir=args.particlesDir,
        hitsDir=args.hitsDir,
        tracksDir=args.tracksDir,
        readTracks=args.readTracks,
        decorators=decorators,
        events=args.events,
        numThreads=args.jobs,
    )
    s.run()
