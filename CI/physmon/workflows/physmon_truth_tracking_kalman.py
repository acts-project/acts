#!/usr/bin/env python3
import tempfile
from pathlib import Path
import shutil

import acts
from truth_tracking_kalman import runTruthTrackingKalman

from physmon_common import makeSetup

setup = makeSetup()

with tempfile.TemporaryDirectory() as temp:
    s = acts.examples.Sequencer(
        events=10000,
        numThreads=-1,
        logLevel=acts.logging.INFO,
        fpeMasks=acts.examples.Sequencer.FpeMask.fromFile(
            Path(__file__).parent.parent / "fpe_masks.yml"
        ),
    )

    tp = Path(temp)
    runTruthTrackingKalman(
        setup.trackingGeometry,
        setup.field,
        setup.digiConfig,
        outputDir=tp,
        s=s,
    )

    s.run()
    del s

    perf_file = tp / "performance_track_fitter.root"
    assert perf_file.exists(), "Performance file not found"
    shutil.copy(perf_file, setup.outdir / "performance_truth_tracking.root")
