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
        fpeMasks=[
            (
                "Fatras/include/ActsFatras/Physics/ElectroMagnetic/BetheHeitler.hpp:66",
                acts.FpeType.FLTUND,
                1,
            ),
            (
                "Fatras/include/ActsFatras/Kernel/detail/SimulationActor.hpp:178",
                acts.FpeType.FLTUND,
                1,
            ),
            (
                "Core/include/Acts/TrackFitting/detail/GsfUtils.hpp:187",
                acts.FpeType.FLTUND,
                1,
            ),
            (
                "Acts/Utilities/GaussianMixtureReduction.hpp:198",
                acts.FpeType.FLTUND,
                1,
            ),
        ],
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
