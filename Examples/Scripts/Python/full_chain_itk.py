#!/usr/bin/env python3
import pathlib, acts, acts.examples
from acts.examples.itk import buildITkGeometry

u = acts.UnitConstants
geo_dir = pathlib.Path("acts-itk")
outputDir = pathlib.Path.cwd()

# acts.examples.dump_args_calls(locals())  # show acts.examples python binding calls
detector, trackingGeometry, decorators = buildITkGeometry(geo_dir)
field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))
rnd = acts.examples.RandomNumbers(seed=42)

from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    SeedingAlgorithm,
    ParticleSmearingSigmas,
    addCKFTracks,
    CKFPerformanceConfig,
)

from acts.examples.itk import itkSeedingAlgConfig

s = acts.examples.Sequencer(events=100, numThreads=-1)
addParticleGun(
    s,
    MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
    EtaConfig(-4.0, 4.0, uniform=True),
    ParticleConfig(1, acts.PdgParticle.eMuon, randomizeCharge=True),
    rnd=rnd,
)
# # Uncomment addPythia8 and ParticleSelectorConfig, instead of addParticleGun, to generate ttbar with mu=200 pile-up.
# addPythia8(
#     s,
#     hardProcess=["Top:qqbar2ttbar=on"],
#     vtxGen=acts.examples.GaussianVertexGenerator(
#         stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
#         mean=acts.Vector4(0, 0, 0, 0),
#     ),
#     rnd=rnd,
#     outputDirRoot=outputDir,
# )
addFatras(
    s,
    trackingGeometry,
    field,
    # ParticleSelectorConfig(eta=(-4.0, 4.0), pt=(150 * u.MeV, None), removeNeutral=True),
    outputDirRoot=outputDir,
    rnd=rnd,
)
addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=geo_dir / "itk-hgtd/itk-smearing-config.json",
    outputDirRoot=outputDir,
    rnd=rnd,
)
addSeeding(
    s,
    trackingGeometry,
    field,
    TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-4.0, 4.0), nHits=(9, None)),
    # SeedingAlgorithm.TruthEstimated,
    # SeedingAlgorithm.TruthSmeared, ParticleSmearingSigmas(pRel=0.01), rnd=rnd,
    *itkSeedingAlgConfig("PixelSpacePoints"),
    geoSelectionConfigFile=geo_dir / "itk-hgtd/geoSelection-ITk.json",
    outputDirRoot=outputDir,
)
addCKFTracks(
    s,
    trackingGeometry,
    field,
    CKFPerformanceConfig(ptMin=1.0 * u.GeV, nMeasurementsMin=6),
    outputDirRoot=outputDir,
)

s.run()
