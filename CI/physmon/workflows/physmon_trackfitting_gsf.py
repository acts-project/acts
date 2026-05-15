#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil

import acts
from truth_tracking_gsf import runTruthTrackingGsf

from physmon_common import makeSetup

setup = makeSetup()

# Paths to pre-simulated data (from physmon_simulation_gsf.py)
simDir = (setup.outdir / "../simulation_gsf/simulation").resolve()
particlesPath = simDir / "particles_simulation.root"
simhitsPath = simDir / "hits.root"

# Verify simulation files exist
assert particlesPath.exists(), f"Simulation not found: {particlesPath}"
assert simhitsPath.exists(), f"SimHits not found: {simhitsPath}"

with tempfile.TemporaryDirectory() as temp:
    s = acts.examples.Sequencer(
        events=10000,
        numThreads=-1,  # Multi-threaded now!
        logLevel=acts.logging.INFO,
    )

    tp = Path(temp)
    runTruthTrackingGsf(
        trackingGeometry=setup.trackingGeometry,
        field=setup.field,
        digiConfigFile=setup.digiConfig,
        outputDir=tp,
        inputParticlePath=particlesPath,
        inputSimHitsPath=simhitsPath,
        decorators=setup.decorators,
        s=s,
    )

    s.run()

    perf_file = tp / "performance_gsf.root"
    assert perf_file.exists(), "Performance file not found"
    shutil.copy(perf_file, setup.outdir / "performance_trackfitting.root")
