#!/usr/bin/env python3
from pathlib import Path
from typing import Optional, Union
import enum
import warnings

from acts.examples import (
    Sequencer,
    ParticleSelector,
    ParticleSmearing,
    TruthVertexFinder,
    VertexFitterAlgorithm,
    IterativeVertexFinderAlgorithm,
    RootParticleReader,
    AdaptiveMultiVertexFinderAlgorithm,
    RootVertexPerformanceWriter,
    RootTrajectorySummaryReader,
    TrackSelector,
    GenericDetector,
)

import acts

from acts import UnitConstants as u

import pythia8


class VertexFinder(enum.Enum):
    Truth = (1,)
    AMVF = (2,)
    Iterative = (3,)


def addVertexFitting(
    s,
    field,
    outputDirRoot: Optional[Union[Path, str]] = None,
    associatedParticles: str = "particles_input",
    trackParameters: str = "estimatedparameters",
    vertexFinder: VertexFinder = VertexFinder.Truth,
    logLevel: Optional[acts.logging.Level] = None,
):
    """This function steers the vertex fitting

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addVertexFitting)
    field : magnetic field
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    associatedParticles : str, "associatedTruthParticles"
        RootVertexPerformanceWriter.inputAssociatedTruthParticles
    vertexFinder : VertexFinder, Truth
        vertexFinder algorithm: one of Truth, AMVF, Iterative
    logLevel : acts.logging.Level, None
        logging level to override setting given in `s`
    """

    def customLogLevel(custom: acts.logging.Level = acts.logging.INFO):
        """override logging level"""
        if logLevel is None:
            return s.config.logLevel
        return acts.logging.Level(max(custom.value, logLevel.value))

    if int(customLogLevel()) <= int(acts.logging.DEBUG):
        acts.examples.dump_args_calls(locals())

    inputParticles = "particles_input"
    outputVertices = "fittedVertices"
    selectedParticles = "particles_selected"

    outputTime = ""
    if vertexFinder == VertexFinder.Truth:
        findVertices = TruthVertexFinder(
            level=customLogLevel(acts.logging.VERBOSE),
            inputParticles=selectedParticles,
            outputProtoVertices="protovertices",
            excludeSecondaries=True,
        )
        s.addAlgorithm(findVertices)
        fitVertices = VertexFitterAlgorithm(
            level=customLogLevel(acts.logging.VERBOSE),
            bField=field,
            inputTrackParameters=trackParameters,
            inputProtoVertices=findVertices.config.outputProtoVertices,
            outputVertices=outputVertices,
        )
        s.addAlgorithm(fitVertices)

    elif vertexFinder == VertexFinder.Iterative:
        findVertices = IterativeVertexFinderAlgorithm(
            level=customLogLevel(),
            bField=field,
            inputTrackParameters=trackParameters,
            outputProtoVertices="protovertices",
            outputVertices=outputVertices,
        )
        s.addAlgorithm(findVertices)
    elif vertexFinder == VertexFinder.AMVF:
        outputTime = "outputTime"
        findVertices = AdaptiveMultiVertexFinderAlgorithm(
            level=customLogLevel(),
            bField=field,
            inputTrackParameters=trackParameters,
            outputProtoVertices="protovertices",
            outputVertices=outputVertices,
            outputTime=outputTime,
        )

        s.addAlgorithm(findVertices)
    else:
        raise RuntimeError("Invalid finder argument")

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()
        if associatedParticles == selectedParticles:
            warnings.warn(
                "Using RootVertexPerformanceWriter with smeared particles is not necessarily supported. "
                "Please get in touch with us"
            )
        s.addWriter(
            RootVertexPerformanceWriter(
                level=customLogLevel(),
                inputAllTruthParticles=inputParticles,
                inputSelectedTruthParticles=selectedParticles,
                inputAssociatedTruthParticles=associatedParticles,
                inputFittedTracks=trackParameters,
                inputVertices=outputVertices,
                inputTime=outputTime,
                treeName="vertexing",
                filePath=str(outputDirRoot / "performance_vertexing.root"),
            )
        )

    return s


def runVertexFitting(
    field,
    outputDir: Path,
    outputRoot: bool = True,
    inputParticlePath: Optional[Path] = None,
    inputTrackSummary: Path = None,
    vertexFinder: VertexFinder = VertexFinder.Truth,
    s=None,
):
    s = s or Sequencer(events=100, numThreads=-1)

    logger = acts.logging.getLogger("VertexFittingExample")

    rnd = acts.examples.RandomNumbers(seed=42)

    inputParticles = "particles_input"
    if inputParticlePath is None:
        logger.info("Generating particles using Pythia8")
        pythia8.addPythia8(s, rnd)
    else:
        logger.info("Reading particles from %s", inputParticlePath.resolve())
        assert inputParticlePath.exists()
        s.addReader(
            RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                particleCollection=inputParticles,
                orderedEvents=False,
            )
        )

    selectedParticles = "particles_selected"
    ptclSelector = ParticleSelector(
        level=acts.logging.INFO,
        inputParticles=inputParticles,
        outputParticles=selectedParticles,
        removeNeutral=True,
        absEtaMax=2.5,
        rhoMax=4.0 * u.mm,
        ptMin=500 * u.MeV,
    )
    s.addAlgorithm(ptclSelector)

    trackParameters = "trackparameters"
    if inputTrackSummary is None or inputParticlePath is None:
        logger.info("Using smeared particles")

        ptclSmearing = ParticleSmearing(
            level=acts.logging.INFO,
            inputParticles=selectedParticles,
            outputTrackParameters=trackParameters,
            randomNumbers=rnd,
        )
        s.addAlgorithm(ptclSmearing)
        associatedParticles = selectedParticles
    else:
        logger.info("Reading track summary from %s", inputTrackSummary.resolve())
        assert inputTrackSummary.exists()
        associatedParticles = "associatedTruthParticles"
        trackSummaryReader = RootTrajectorySummaryReader(
            level=acts.logging.VERBOSE,
            outputTracks="fittedTrackParameters",
            outputParticles=associatedParticles,
            filePath=str(inputTrackSummary.resolve()),
            orderedEvents=False,
        )
        s.addReader(trackSummaryReader)

        s.addAlgorithm(
            TrackSelector(
                level=acts.logging.INFO,
                inputTrackParameters=trackSummaryReader.config.outputTracks,
                outputTrackParameters=trackParameters,
                outputTrackIndices="outputTrackIndices",
                removeNeutral=True,
                absEtaMax=2.5,
                loc0Max=4.0 * u.mm,  # rho max
                ptMin=500 * u.MeV,
            )
        )

    logger.info("Using vertex finder: %s", vertexFinder.name)

    return addVertexFitting(
        s,
        field,
        outputDirRoot=outputDir if outputRoot else None,
        associatedParticles=associatedParticles,
        trackParameters=trackParameters,
        vertexFinder=vertexFinder,
    )


if "__main__" == __name__:
    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    inputParticlePath = Path("particles.root")
    if not inputParticlePath.exists():
        inputParticlePath = None

    inputTrackSummary = None
    for p in ("tracksummary_fitter.root", "tracksummary_ckf.root"):
        p = Path(p)
        if p.exists():
            inputTrackSummary = p
            break

    runVertexFitting(
        field,
        vertexFinder=VertexFinder.Truth,
        inputParticlePath=inputParticlePath,
        inputTrackSummary=inputTrackSummary,
        outputDir=Path.cwd(),
    ).run()
