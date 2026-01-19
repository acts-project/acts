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
from physmon_common import makeSetup

u = acts.UnitConstants
setup = makeSetup()

# Create output directory for simulation files
outputDir = setup.outdir / "simulation"
outputDir.mkdir(exist_ok=True)

# Single-threaded sequencer (Geant4 requirement)
s = acts.examples.Sequencer(
    events=10000,
    numThreads=1,  # Geant4 must run single-threaded
    logLevel=acts.logging.INFO,
)

# Add context decorators
for d in setup.decorators:
    s.addContextDecorator(d)

rnd = acts.examples.RandomNumbers(seed=42)

# Particle gun: electrons, matching truth_tracking_gsf.py configuration
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

# Geant4 simulation
addGeant4(
    s,
    setup.detector,
    setup.trackingGeometry,
    setup.field,
    rnd,
    killVolume=setup.trackingGeometry.highestTrackingVolume,
    killAfterTime=25 * u.ns,
    killSecondaries=True,
    outputDirRoot=outputDir,  # Enable ROOT output
)

# Run simulation
s.run()

print(f"Simulation complete. Output files:")
print(f"  - {outputDir / 'particles_simulation.root'}")
print(f"  - {outputDir / 'hits.root'}")
