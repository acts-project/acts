#!/usr/bin/env python3
from pathlib import Path
from typing import Optional, Union

import acts.examples
import acts
from acts import UnitConstants as u


def addExaTrkx(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    geometrySelection: Union[Path, str],
    onnxModelDir: Union[Path, str],
    outputDirRoot: Optional[Union[Path, str]] = None,
) -> acts.examples.Sequencer:

    # Run the particle selection
    # The pre-selection will select truth particles satisfying provided criteria
    # from all particles read in by particle reader for further processing. It
    # has no impact on the truth hits themselves
    s.addAlgorithm(
        acts.examples.TruthSeedSelector(
            level=acts.logging.INFO,
            ptMin=500 * u.MeV,
            nHitsMin=9,
            inputParticles="particles_initial",
            inputMeasurementParticlesMap="measurement_particle_map",
            outputParticles="particles_seed_selected",
        )
    )

    # Create space points
    s.addAlgorithm(
        acts.examples.SpacePointMaker(
            level=acts.logging.INFO,
            inputSourceLinks="source_links",
            inputMeasurements="measurements",
            outputSpacePoints="spacepoints",
            trackingGeometry=trackingGeometry,
            geometrySelection=acts.examples.readJsonGeometryList(
                str(geometrySelection)
            ),
        )
    )

    # Setup the track finding algorithm with ExaTrkX
    # It takes all the source links created from truth hit smearing, seeds from
    # truth particle smearing and source link selection config
    exaTrkxFinding = acts.examples.ExaTrkXTrackFinding(
        inputMLModuleDir=str(onnxModelDir),
        spacepointFeatures=3,
        embeddingDim=8,
        rVal=1.6,
        knnVal=500,
        filterCut=0.21,
    )

    s.addAlgorithm(
        acts.examples.TrackFindingMLBasedAlgorithm(
            level=acts.logging.INFO,
            inputSpacePoints="spacepoints",
            outputProtoTracks="protoTracks",
            trackFinderML=exaTrkxFinding,
        )
    )

    # Write truth track finding / seeding performance
    if outputDirRoot is not None:
        s.addWriter(
            acts.examples.TrackFinderPerformanceWriter(
                level=acts.logging.INFO,
                inputProtoTracks="protoTracks",
                inputParticles="particles_initial",  # the original selected particles after digitization
                inputMeasurementParticlesMap="measurement_particle_map",
                filePath=str(Path(outputDirRoot) / "performance_seeding_trees.root"),
            )
        )

    return s


if "__main__" == __name__:
    import os
    from digitization import configureDigitization

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    detector, trackingGeometry, decorators = acts.examples.GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    inputParticlePath = Path("particles.root")
    if not inputParticlePath.exists():
        inputParticlePath = None

    srcdir = Path(__file__).resolve().parent.parent.parent.parent
    geometrySelection = (
        srcdir
        / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
    )
    assert geometrySelection.exists()

    onnxdir = Path(os.getcwd()) / "onnx_models"
    assert onnxdir.exists()

    s = acts.examples.Sequencer(events=100, numThreads=-1)
    s.config.logLevel = acts.logging.INFO

    rnd = acts.examples.RandomNumbers()
    outputDir = Path(os.getcwd())

    s = configureDigitization(
        trackingGeometry,
        field,
        outputDir,
        inputParticlePath,
        outputRoot=True,
        outputCsv=True,
        s=s,
    )

    s = addExaTrkx(s, trackingGeometry, geometrySelection, onnxdir, outputDir)

    s.run()
