#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil

import acts
from truth_tracking_gsf import runTruthTrackingGsf

from physmon_common import makeSetup

setup = makeSetup()

with tempfile.TemporaryDirectory() as temp:
    s = acts.examples.Sequencer(
        events=10000,
        numThreads=-1,
        logLevel=acts.logging.INFO,
    )

    tp = Path(temp)
    runTruthTrackingGsf(
        setup.trackingGeometry,
        setup.field,
        setup.digiConfig,
        outputDir=tp,
        s=s,
    )

    s.run()

    perf_file = tp / "performance_gsf.root"
    assert perf_file.exists(), "Performance file not found"
    shutil.copy(perf_file, setup.outdir / "performance_trackfitting.root")
