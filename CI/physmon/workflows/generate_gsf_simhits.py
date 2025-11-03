#!/usr/bin/env python3

from pathlib import Path

import acts
from acts.examples.simulation import (
    addParticleGun,
    ParticleConfig,
    EtaConfig,
    PhiConfig,
    MomentumConfig,
    addGeant4,
)
from acts.examples.odd import getOpenDataDetector

u = acts.UnitConstants

# Hardcoded configuration matching runTruthTrackingGsf
srcdir = Path(__file__).resolve().parent.parent.parent.parent
outputDir = srcdir / "CI/physmon/gsf_simhits"

# Setup detector and field
matDeco = acts.IMaterialDecorator.fromFile(
    srcdir / "thirdparty/OpenDataDetector/data/odd-material-maps.root",
    level=acts.logging.INFO,
)
detector = getOpenDataDetector(matDeco)
trackingGeometry = detector.trackingGeometry()
field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

# Create sequencer (single-threaded for Geant4)
s = acts.examples.Sequencer(
    events=10000,
    numThreads=1,
    logLevel=acts.logging.INFO,
)

# Generate electrons with same config as runTruthTrackingGsf
addParticleGun(
    s,
    ParticleConfig(num=1, pdg=acts.PdgParticle.eElectron, randomizeCharge=True),
    EtaConfig(-3.0, 3.0, uniform=True),
    MomentumConfig(1.0 * u.GeV, 100.0 * u.GeV, transverse=True),
    PhiConfig(0.0, 360.0 * u.degree),
    vtxGen=acts.examples.GaussianVertexGenerator(
        mean=acts.Vector4(0, 0, 0, 0),
        stddev=acts.Vector4(0, 0, 0, 0),
    ),
    multiplicity=1,
    rnd=rnd,
)

# Run Geant4 simulation
addGeant4(
    s,
    detector,
    trackingGeometry,
    field,
    rnd=rnd,
    killVolume=trackingGeometry.highestTrackingVolume,
    killAfterTime=25 * u.ns,
    outputDirRoot=outputDir,
)

# Run simulation
outputDir.mkdir(parents=True, exist_ok=True)
s.run()

print(f"Generated simulation files in {outputDir}")
print(f"  - particles_simulation.root")
print(f"  - hits.root")
