#!/usr/bin/env python3
from pathlib import Path
from typing import Optional

from acts.examples import Sequencer, GenericDetector, RootParticleReader

import acts.examples

import acts

from acts import UnitConstants as u


def runExaTrkX(
    trackingGeometry,
    decorators,
    geometrySelection: Path,
    digiConfigFile: Path,
    field,
    outputDir: Path,
    outputCsv=True,
    inputParticlePath: Optional[Path] = None,
    s=None,
):
    s = s or Sequencer(events=1, numThreads=1)

    logger = acts.logging.getLogger("ExaTrkXExample")

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    if inputParticlePath is None:
        logger.info("Generating particles using particle gun")

        evGen = acts.examples.EventGenerator(
            level=acts.logging.INFO,
            generators=[
                acts.examples.EventGenerator.Generator(
                    multiplicity=acts.examples.FixedMultiplicityGenerator(n=2),
                    vertex=acts.examples.GaussianVertexGenerator(
                        stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
                    ),
                    particles=acts.examples.ParametricParticleGenerator(
                        p=(1 * u.GeV, 10 * u.GeV),
                        eta=(-2, 2),
                        phi=(0, 360 * u.degree),
                        randomizeCharge=True,
                        numParticles=4,
                    ),
                )
            ],
            outputParticles="particles_input",
            randomNumbers=rnd,
        )
        s.addReader(evGen)
        inputParticles = evGen.config.outputParticles
    else:
        logger.info("Reading particles from %s", inputParticlePath.resolve())
        assert inputParticlePath.exists()
        inputParticles = "particles_read"
        s.addReader(
            RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                particleCollection=inputParticles,
                orderedEvents=False,
            )
        )

    # Selector
    selector = acts.examples.ParticleSelector(
        level=acts.logging.INFO,
        inputParticles=inputParticles,
        outputParticles="particles_selected",
    )
    s.addAlgorithm(selector)

    # Simulation
    simAlg = acts.examples.FatrasSimulation(
        level=acts.logging.INFO,
        inputParticles=selector.config.outputParticles,
        outputParticlesInitial="particles_initial",
        outputParticlesFinal="particles_final",
        outputSimHits="simhits",
        randomNumbers=rnd,
        trackingGeometry=trackingGeometry,
        magneticField=field,
        generateHitsOnSensitive=True,
    )
    s.addAlgorithm(simAlg)

    # Run the sim hits smearing
    digiCfg = acts.examples.DigitizationConfig(
        acts.examples.readDigiConfigFromJson(str(digiConfigFile)),
        trackingGeometry=trackingGeometry,
        randomNumbers=rnd,
        inputSimHits=simAlg.config.outputSimHits,
    )
    digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, acts.logging.INFO)
    s.addAlgorithm(digiAlg)

    # Run the particle selection
    # The pre-selection will select truth particles satisfying provided criteria
    # from all particles read in by particle reader for further processing. It
    # has no impact on the truth hits themselves
    selAlg = acts.examples.TruthSeedSelector(
        level=acts.logging.INFO,
        ptMin=500 * u.MeV,
        nHitsMin=9,
        inputParticles=simAlg.config.outputParticlesInitial,
        inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
        outputParticles="particles_seed_selected",
    )
    s.addAlgorithm(selAlg)

    inputParticles = selAlg.config.outputParticles

    # Create space points
    spAlg = acts.examples.SpacePointMaker(
        level=acts.logging.INFO,
        inputSourceLinks=digiAlg.config.outputSourceLinks,
        inputMeasurements=digiAlg.config.outputMeasurements,
        outputSpacePoints="spacepoints",
        trackingGeometry=trackingGeometry,
        geometrySelection=acts.examples.readJsonGeometryList(str(geometrySelection)),
    )
    s.addAlgorithm(spAlg)

    # Setup the track finding algorithm with ExaTrkX
    # It takes all the source links created from truth hit smearing, seeds from
    # truth particle smearing and source link selection config
    exaTrkxFinding = acts.examples.ExaTrkXTrackFinding(
        inputMLModuleDir="/home/benjamin/Documents/acts_project/gnn_integration/run/onnx_models",
        spacepointFeatures=3,
        embeddingDim=8,
        rVal=1.6,
        knnVal=500,
        filterCut=0.21,
    )

    trackFinderAlg = acts.examples.TrackFindingAlgorithmExaTrkX(
        level=acts.logging.INFO,
        inputSpacePoints="spacepoints",
        outputProtoTracks="protoTracks",
        trackFinderML=exaTrkxFinding,
    )
    s.addAlgorithm(trackFinderAlg)

    # Write truth track finding / seeding performance
    trackFinderPerformanceWriter = acts.examples.TrackFinderPerformanceWriter(
        level=acts.logging.INFO,
        inputProtoTracks="protoTracks",
        inputParticles=inputParticles,  # the original selected particles after digitization
        inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
        filePath=str(outputDir / "performance_seeding_trees.root"),
    )
    s.addWriter(trackFinderPerformanceWriter)

    return s


if "__main__" == __name__:
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    inputParticlePath = Path("particles.root")
    if not inputParticlePath.exists():
        inputParticlePath = None

    runExaTrkX(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=srcdir
        / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json",
        digiConfigFile=srcdir
        / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json",
        outputDir=Path.cwd(),
    ).run()
