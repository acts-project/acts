#!/usr/bin/python
import pathlib, acts, acts.examples
from pathlib import Path
from typing import Optional, Union
from acts import UnitConstants as u
from acts.examples.geant4 import TelescopeG4DetectorConstructionFactory 
from acts.examples import TelescopeDetector

from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addFatras,
    addGeant4,
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addTelescopeSeeding,
    addCKFTracks,
    TrackSelectorConfig,
    CkfConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
)

teleG4Config=TelescopeDetector.Config();
teleG4Config.bounds=[200, 200]
teleG4Config.positions=[30, 60, 90, 105, 120, 150, 180]
teleG4Config.thickness = [80, 80, 80, 1, 80, 80, 80]
teleG4Config.binValue=0

u = acts.UnitConstants

detector, trackingGeometry, decorators = acts.examples.TelescopeDetector.create(
    bounds=[200, 200],
    positions=[30, 60, 90, 105, 120, 150, 180],
    thickness=[0.08, 0.08, 0.08, 0.001, 0.08, 0.08, 0.08],
    binValue=0,
)

# ParticleGun Configurationi
# multiplicity
'''
mul = int(input("multiplicity = "))
if mul == 1:
    outputDir = Path.cwd() / "wot_multiplicity_1"
elif mul == 2:
    outputDir = Path.cwd() / "wot_multiplicity_2"
elif mul == 3:
    outputDir = Path.cwd() / "wot_multiplicity_3"
elif mul == 4:
    outputDir = Path.cwd() / "wot_multiplicity_4"
elif mul == 5:
    outputDir = Path.cwd() / "wot_multiplicity_5"
else:
    print("multiplicity error")
'''

# position stddev
'''
stddev_p = float(input("value of Y-Z stddev"))
if stddev_p == 0:
    outputDir = Path.cwd() / "wot_pos_stddev_0"
elif stddev_p == 20:
    outputDir = Path.cwd() / "wot_pos_stddev_20"
elif stddev_p == 40:
    outputDir = Path.cwd() / "wot_pos_stddev_40"
elif stddev_p == 60:
    outputDir = Path.cwd() / "wot_pos_stddev_60"
elif stddev_p == 80:
    outputDir = Path.cwd() / "wot_pos_stddev_80"
elif stddev_p == 100:
    outputDir = Path.cwd() / "wot_pos_stddev_100"
elif stddev_p == 120:
    outputDir = Path.cwd() / "wot_pos_stddev_120"
elif stddev_p == 140:
    outputDir = Path.cwd() / "wot_pos_stddev_140"
elif stddev_p == 160:
    outputDir = Path.cwd() / "wot_pos_stddev_160"
elif stddev_p == 180:
    outputDir = Path.cwd() / "wot_pos_stddev_180"
elif stddev_p == 200:
    outputDir = Path.cwd() / "wot_pos_stddev_200"
else:
    print("pos_stddev error")
'''

# time stddev

stddev_t = float(input("value of time stddev"))
if stddev_t == 0:
    outputDir = Path.cwd() / "wot_time_stddev_0"
elif stddev_t == 2:
    outputDir = Path.cwd() / "wot_time_stddev_2"
elif stddev_t == 4:
    outputDir = Path.cwd() / "wot_time_stddev_4"
elif stddev_t == 6:
    outputDir = Path.cwd() / "wot_time_stddev_6"
elif stddev_t == 8:
    outputDir = Path.cwd() / "wot_time_stddev_8"
elif stddev_t == 10:
    outputDir = Path.cwd() / "wot_time_stddev_10"
else:
    print("time_stddev error")


field = acts.ConstantBField(acts.Vector3(0, 0, 0)) # u.T
if not outputDir.exists():
    outputDir.mkdir()

rnd = acts.examples.RandomNumbers(seed=42)
s = acts.examples.Sequencer(events=100, numThreads=1, outputDir=str(outputDir))

addParticleGun(
    s,
    MomentumConfig(4 * u.GeV, 4 * u.GeV, transverse=True),
    EtaConfig(-0.0125, 0.0125, uniform=True),
    PhiConfig(-0.0 * u.degree, 0.71 * u.degree),
    ParticleConfig(1, acts.PdgParticle.eMuon, randomizeCharge=False),
    multiplicity=3,
    #multiplicity=mul,
    rnd=rnd,
    #vtxGen=acts.examples.GaussianVertexGenerator(mean=acts.Vector4(0, 0, 0, 0), stddev=acts.Vector4(0.*u.mm, 50.*u.um, 50.*u.um, 4*u.ns)),
    #vtxGen=acts.examples.GaussianVertexGenerator(mean=acts.Vector4(0, 0, 0, 0), stddev=acts.Vector4(0.*u.mm, stddev_p*u.um, stddev_p*u.um, 4*u.ns)),
    vtxGen=acts.examples.GaussianVertexGenerator(mean=acts.Vector4(0, 0, 0, 0), stddev=acts.Vector4(5.*u.mm, 5.*u.mm, 5.*u.um, stddev_t*u.ns)),
)


addFatras(
    s,
    trackingGeometry,
    field,
    rnd=rnd,
    preSelectParticles=ParticleSelectorConfig(
        rho=(0.0 * u.mm, 300.0 * u.mm),
        absZ=(0.0 * u.mm, 200.0 * u.mm),
        eta=(-3, 3),
        pt=(1 * u.GeV, None),
        removeNeutral=True,
    ),
    outputDirRoot=outputDir,
)

'''
addGeant4(
    s,
    detector=None,
    trackingGeometry=trackingGeometry,
    field=field,
    rnd=rnd,
    volumeMappings = ["Layer #0 Phys"],
    g4DetectorConstructionFactory=TelescopeG4DetectorConstructionFactory(teleG4Config),
    preSelectParticles=ParticleSelectorConfig(
            rho=(0.0, 300 * u.mm),
            absZ=(-200.0, 200.0 * u.m),
            eta=(-3.0, 3.0),
            pt=(1 * u.GeV, None),
            removeNeutral=True,
    ),
    outputDirRoot=outputDir,
    killVolume=trackingGeometry.worldVolume,
    killAfterTime=1000 * u.ns,
)
'''

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=Path("../Examples/Algorithms/Digitization/share/default-smearing-config-telescope.json"),
    outputDirRoot=outputDir,
    rnd=rnd,
    # logLevel=acts.logging.VERBOSE
)


addTelescopeSeeding(
    s,
    trackingGeometry,
    initialSigmas={1, 1, 1, 1, 0.1, 1},
    initialVarInflation={1, 1, 1, 1, 1, 1},
    outputDirRoot=outputDir,
    # logLevel=acts.logging.VERBOSE,
)

addCKFTracks(                                                                                                                                                   s,                                                                                                                                                          trackingGeometry,                                                                                                                                           field,                                                                                                                                                      CkfConfig(                                                                                                                                                      chi2CutOff=15.0,                                                                                                                                            numMeasurementsCutOff=10,                                                                                                                               ),                                                                                                                                                          outputDirRoot=outputDir,                                                                                                                                    logLevel=acts.logging.DEBUG
)

s.run()
