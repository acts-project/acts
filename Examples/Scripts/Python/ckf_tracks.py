#!/usr/bin/env python3
from pathlib import Path
from typing import Optional, Union
from collections import namedtuple

from acts.examples import Sequencer, GenericDetector, RootParticleReader

import acts

from acts import UnitConstants as u

CKFPerformanceConfig = namedtuple(
    "CKFPerformanceConfig",
    ["truthMatchProbMin", "nMeasurementsMin", "ptMin"],
    defaults=[None] * 3,
)


@acts.examples.NamedTypeArgs(
    CKFPerformanceConfigArg=CKFPerformanceConfig,
)
def addCKFTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    CKFPerformanceConfigArg: CKFPerformanceConfig = CKFPerformanceConfig(),
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    selectedParticles: str = "truth_seeds_selected",
) -> acts.examples.Sequencer:
    """This function steers the seeding

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    CKFPerformanceConfigArg : CKFPerformanceConfig(truthMatchProbMin, nMeasurementsMin, ptMin)
        CKFPerformanceWriter configuration.
        Defaults specified in Examples/Io/Performance/ActsExamples/Io/Performance/CKFPerformanceWriter.hpp
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    selectedParticles : str, "truth_seeds_selected"
        CKFPerformanceWriter truth input
    """

    if int(s.config.logLevel) <= int(acts.logging.DEBUG):
        acts.examples.dump_args_calls(locals())

    # Setup the track finding algorithm with CKF
    # It takes all the source links created from truth hit smearing, seeds from
    # truth particle smearing and source link selection config
    trackFinder = acts.examples.TrackFindingAlgorithm(
        level=s.config.logLevel,
        measurementSelectorCfg=acts.MeasurementSelector.Config(
            [(acts.GeometryIdentifier(), ([], [15.0], [10]))]
        ),
        inputMeasurements="measurements",
        inputSourceLinks="sourcelinks",
        inputInitialTrackParameters="estimatedparameters",
        outputTrajectories="trajectories",
        findTracks=acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
            trackingGeometry, field
        ),
    )
    s.addAlgorithm(trackFinder)

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

        # write track states from CKF
        trackStatesWriter = acts.examples.RootTrajectoryStatesWriter(
            level=s.config.logLevel,
            inputTrajectories=trackFinder.config.outputTrajectories,
            # @note The full particles collection is used here to avoid lots of warnings
            # since the unselected CKF track might have a majority particle not in the
            # filtered particle collection. This could be avoided when a seperate track
            # selection algorithm is used.
            inputParticles="particles_selected",
            inputSimHits="simhits",
            inputMeasurementParticlesMap="measurement_particles_map",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDirRoot / "trackstates_ckf.root"),
            treeName="trackstates",
        )
        s.addWriter(trackStatesWriter)

        # write track summary from CKF
        trackSummaryWriter = acts.examples.RootTrajectorySummaryWriter(
            level=s.config.logLevel,
            inputTrajectories=trackFinder.config.outputTrajectories,
            # @note The full particles collection is used here to avoid lots of warnings
            # since the unselected CKF track might have a majority particle not in the
            # filtered particle collection. This could be avoided when a seperate track
            # selection algorithm is used.
            inputParticles="particles_selected",
            inputMeasurementParticlesMap="measurement_particles_map",
            filePath=str(outputDirRoot / "tracksummary_ckf.root"),
            treeName="tracksummary",
        )
        s.addWriter(trackSummaryWriter)

        # Write CKF performance data
        ckfPerfWriter = acts.examples.CKFPerformanceWriter(
            level=s.config.logLevel,
            inputParticles=selectedParticles,
            inputTrajectories=trackFinder.config.outputTrajectories,
            inputMeasurementParticlesMap="measurement_particles_map",
            **acts.examples.defaultKWArgs(
                # The bottom seed could be the first, second or third hits on the truth track
                nMeasurementsMin=CKFPerformanceConfigArg.nMeasurementsMin,
                ptMin=CKFPerformanceConfigArg.ptMin,
                truthMatchProbMin=CKFPerformanceConfigArg.truthMatchProbMin,
            ),
            filePath=str(outputDirRoot / "performance_ckf.root"),
        )
        s.addWriter(ckfPerfWriter)

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()
        acts.logging.getLogger("CKFExample").info("Writing CSV files")
        csvMTJWriter = acts.examples.CsvMultiTrajectoryWriter(
            level=s.config.logLevel,
            inputTrajectories=trackFinder.config.outputTrajectories,
            inputMeasurementParticlesMap="measurement_particles_map",
            outputDir=str(outputDirCsv),
        )
        s.addWriter(csvMTJWriter)

    return s


def runCKFTracks(
    trackingGeometry,
    decorators,
    geometrySelection: Path,
    digiConfigFile: Path,
    field,
    outputDir: Path,
    truthSmearedSeeded=False,
    truthEstimatedSeeded=False,
    outputCsv=True,
    inputParticlePath: Optional[Path] = None,
    s=None,
):

    from particle_gun import addParticleGun, EtaConfig, PhiConfig, ParticleConfig
    from fatras import addFatras
    from digitization import addDigitization
    from seeding import (
        addSeeding,
        TruthSeedRanges,
        ParticleSmearingSigmas,
        SeedfinderConfigArg,
        SeedingAlgorithm,
        TrackParamsEstimationConfig,
    )

    s = s or acts.examples.Sequencer(
        events=100, numThreads=-1, logLevel=acts.logging.INFO
    )
    for d in decorators:
        s.addContextDecorator(d)
    rnd = acts.examples.RandomNumbers(seed=42)
    outputDir = Path(outputDir)

    if inputParticlePath is None:
        s = addParticleGun(
            s,
            EtaConfig(-2.0, 2.0),
            ParticleConfig(4, acts.PdgParticle.eMuon, True),
            PhiConfig(0.0, 360.0 * u.degree),
            multiplicity=2,
            rnd=rnd,
        )
    else:
        acts.logging.getLogger("CKFExample").info(
            "Reading particles from %s", inputParticlePath.resolve()
        )
        assert inputParticlePath.exists()
        s.addReader(
            RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                particleCollection="particles_input",
                orderedEvents=False,
            )
        )

    s = addFatras(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
    )

    s = addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
    )

    s = addSeeding(
        s,
        trackingGeometry,
        field,
        TruthSeedRanges(pt=(500.0 * u.MeV, None), nHits=(9, None)),
        ParticleSmearingSigmas(pRel=0.01),  # only used by SeedingAlgorithm.TruthSmeared
        SeedfinderConfigArg(
            r=(None, 200 * u.mm),  # rMin=default, 33mm
            deltaR=(1 * u.mm, 60 * u.mm),
            collisionRegion=(-250 * u.mm, 250 * u.mm),
            z=(-2000 * u.mm, 2000 * u.mm),
            maxSeedsPerSpM=1,
            sigmaScattering=50,
            radLengthPerSeed=0.1,
            minPt=500 * u.MeV,
            bFieldInZ=1.99724 * u.T,
            impactMax=3 * u.mm,
        ),
        TrackParamsEstimationConfig(deltaR=(10.0 * u.mm, None)),
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared
        if truthSmearedSeeded
        else SeedingAlgorithm.TruthEstimated
        if truthEstimatedSeeded
        else SeedingAlgorithm.Default,
        geoSelectionConfigFile=geometrySelection,
        outputDirRoot=outputDir,
        rnd=rnd,  # only used by SeedingAlgorithm.TruthSmeared
    )

    s = addCKFTracks(
        s,
        trackingGeometry,
        field,
        CKFPerformanceConfig(ptMin=400.0 * u.MeV, nMeasurementsMin=6),
        outputDirRoot=outputDir,
        outputDirCsv=outputDir / "csv" if outputCsv else None,
    )

    return s


if "__main__" == __name__:
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    inputParticlePath = Path("particles.root")
    if not inputParticlePath.exists():
        inputParticlePath = None

    runCKFTracks(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=srcdir
        / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json",
        digiConfigFile=srcdir
        / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json",
        outputCsv=True,
        truthSmearedSeeded=False,
        truthEstimatedSeeded=False,
        inputParticlePath=inputParticlePath,
        outputDir=Path.cwd(),
    ).run()
