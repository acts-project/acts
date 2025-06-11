#!/usr/bin/env python3
import pathlib, acts, acts.examples, acts.examples.itk
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addPythia8,
    ParticleSelectorConfig,
    addGenParticleSelection,
    addFatras,
    addDigitization,
    addDigiParticleSelection,
)
from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    addCKFTracks,
    CkfConfig,
    TrackSelectorConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
    addHoughVertexFinding,
)

u = acts.UnitConstants
geo_dir = pathlib.Path("acts-itk")
outputDir = pathlib.Path.cwd() / "itk_output"
# acts.examples.dump_args_calls(locals())  # show acts.examples python binding calls

detector = acts.examples.itk.buildITkGeometry(geo_dir, logLevel=acts.logging.INFO)
trackingGeometry = detector.trackingGeometry()
# decorators = detector.contextDecorators()

field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))
rnd = acts.examples.RandomNumbers(seed=42)

# this will produce a problematic event, please avoid it
# s = acts.examples.Sequencer(events=1, numThreads=1, outputDir=str(outputDir), skip=34935)

s = acts.examples.Sequencer(events=1, numThreads=1, outputDir=str(outputDir), skip=1)

addPythia8(
    s,
    hardProcess=["Top:qqbar2ttbar=on"],
    npileup=0,
    # beam=acts.PdgParticle.eLead,
    cmsEnergy=5.0 * u.TeV,
    vtxGen=acts.examples.GaussianVertexGenerator(
        # stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
        stddev=acts.Vector4(0.5 * u.mm, 0.5 * u.mm, 50. * u.mm, 0.0 * u.ns),
        mean=acts.Vector4(0, 0, 0, 0),
    ),
    rnd=rnd,
    outputDirRoot=outputDir,
    # logLevel=acts.logging.VERBOSE,
)

addGenParticleSelection(
    s,
    ParticleSelectorConfig(
        rho=(0.0 * u.mm, 28.0 * u.mm),
        absZ=(0.0 * u.mm, 1.0 * u.m),
        eta=(-4.0, 4.0),
        pt=(150 * u.MeV, None),
    ),
    # logLevel=acts.logging.VERBOSE,
)

addFatras(
    s,
    trackingGeometry,
    field,
    rnd=rnd,
    outputDirRoot=outputDir,
    # logLevel=acts.logging.VERBOSE,
)

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=geo_dir / "itk-hgtd/itk-smearing-config.json",
    outputDirRoot=outputDir,
    rnd=rnd,
    # logLevel=acts.logging.VERBOSE,
)

addDigiParticleSelection(
    s,
    ParticleSelectorConfig(
        pt=(1.0 * u.GeV, None),
        eta=(-4.0, 4.0),
        measurements=(9, None),
        removeNeutral=True,
    ),
    # logLevel=acts.logging.VERBOSE,
)

addSeeding(
    s,
    trackingGeometry,
    field,
    seedingAlgorithm=SeedingAlgorithm.Default,
    *acts.examples.itk.itkSeedingAlgConfig(
        acts.examples.itk.InputSpacePointsType.PixelSpacePoints
    ),
    initialSigmas=[
        1 * u.mm,
        1 * u.mm,
        1 * u.degree,
        1 * u.degree,
        0.1 * u.e / u.GeV,
        1 * u.ns,
    ],
    initialSigmaPtRel=0.1,
    initialVarInflation=[1.0] * 6,
    geoSelectionConfigFile=geo_dir / "itk-hgtd/geoSelection-ITk.json",
    outputDirRoot=outputDir,
    logLevel=acts.logging.INFO,
)

addHoughVertexFinding(
    s,
    # logLevel=acts.logging.VERBOSE,
    outputDirRoot=outputDir,
    outputVertices="fittedHoughVertices"
)

addCKFTracks(
    s,
    trackingGeometry,
    field,
    trackSelectorConfig=(
        # fmt: off
        TrackSelectorConfig(absEta=(None, 2.0), pt=(0.9 * u.GeV, None), nMeasurementsMin=9, maxHoles=2, maxOutliers=2, maxSharedHits=2),
        TrackSelectorConfig(absEta=(None, 2.6), pt=(0.4 * u.GeV, None), nMeasurementsMin=8, maxHoles=2, maxOutliers=2, maxSharedHits=2),
        TrackSelectorConfig(absEta=(None, 4.0), pt=(0.4 * u.GeV, None), nMeasurementsMin=7, maxHoles=2, maxOutliers=2, maxSharedHits=2),
        # fmt: on
    ),
    ckfConfig=CkfConfig(
        seedDeduplication=True,
        stayOnSeed=True,
        numMeasurementsCutOff=3,
        # ITk volumes from Noemi's plot
        pixelVolumes=[8, 9, 10, 13, 14, 15, 16, 18, 19, 20],
        stripVolumes=[22, 23, 24],
        maxPixelHoles=1,
        maxStripHoles=2,
    ),
    outputDirRoot=outputDir,
    logLevel=acts.logging.INFO,
)

addAmbiguityResolution(
    s,
    AmbiguityResolutionConfig(
        maximumSharedHits=3,
        maximumIterations=10000,
        nMeasurementsMin=6,
    ),
    outputDirRoot=outputDir,
    logLevel=acts.logging.INFO,
)

addVertexFitting(
    s,
    field,
    # logLevel=acts.logging.VERBOSE,
    vertexFinder=VertexFinder.AMVF,
    outputDirRoot=outputDir,
)



s.run()
