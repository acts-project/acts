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
        events=500,
        numThreads=-1,
        logLevel=acts.logging.INFO,
        fpeMasks=acts.examples.Sequencer.FpeMask.fromFile(
            Path(__file__).parent.parent / "fpe_masks.yml"
        ),
    )

    tp = Path(temp)
    runTruthTrackingGsf(
        setup.trackingGeometry,
        setup.digiConfig,
        setup.field,
        outputDir=tp,
        s=s,
    )

    s.run()
    del s

    perf_file = tp / "performance_gsf.root"
    assert perf_file.exists(), "Performance file not found"
    shutil.copy(perf_file, setup.outdir / "performance_gsf.root")
