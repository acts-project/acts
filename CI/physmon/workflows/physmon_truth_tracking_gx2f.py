#!/usr/bin/env python3
import tempfile
from pathlib import Path
import shutil

import acts
from truth_tracking_gx2f import runTruthTrackingGx2f

from physmon_common import makeSetup

setup = makeSetup()

with tempfile.TemporaryDirectory() as temp:
    s = acts.examples.Sequencer(
        events=10000,
        numThreads=-1,
        logLevel=acts.logging.INFO,
    )

    tp = Path(temp)
    runTruthTrackingGx2f(
        trackingGeometry=setup.trackingGeometry,
        field=setup.field,
        outputDir=tp,
        digiConfigFile=setup.digiConfig,
        s=s,
    )

    s.run()
    del s

    perf_file = tp / "performance_gx2f.root"
    assert perf_file.exists(), "Performance file not found"
    shutil.copy(perf_file, setup.outdir / "performance_gx2f.root")
