#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil

import acts
from truth_tracking_kalman import runTruthTrackingKalman
from truth_tracking_gx2f import runTruthTrackingGx2f

from physmon_common import makeSetup

setup = makeSetup()

digiConfigFile = setup.digiConfig
nSkip = 0
nEvents = 100000
numThreads = -1

with tempfile.TemporaryDirectory() as temp:
    s = acts.examples.Sequencer(
        skip=nSkip,
        events=nEvents,
        numThreads=numThreads,
        logLevel=acts.logging.INFO,
        trackFpes=True,
    )

    tp = Path(temp)
    runTruthTrackingKalman(
        trackingGeometry=setup.trackingGeometry,
        field=setup.field,
        digiConfigFile=digiConfigFile,
        outputDir=tp,
        s=s,
    )

    s.run()
    del s

    perf_file = tp / "performance_kf.root"
    assert perf_file.exists(), "Performance file not found"
    shutil.copy(perf_file, setup.outdir / "performance_trackfitting_kf.root")

with tempfile.TemporaryDirectory() as temp:
    s = acts.examples.Sequencer(
        skip=nSkip,
        events=nEvents,
        numThreads=numThreads,
        logLevel=acts.logging.INFO,
        trackFpes=True,
    )

    tp = Path(temp)
    runTruthTrackingGx2f(
        trackingGeometry=setup.trackingGeometry,
        field=setup.field,
        digiConfigFile=digiConfigFile,
        outputDir=tp,
        s=s,
    )

    s.run()
    del s

    perf_file = tp / "performance_gx2f.root"
    assert perf_file.exists(), "Performance file not found"
    shutil.copy(perf_file, setup.outdir / "performance_trackfitting_gx2f.root")
