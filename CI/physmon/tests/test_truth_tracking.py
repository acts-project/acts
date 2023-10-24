from pathlib import Path

import acts.examples

from truth_tracking_kalman import runTruthTrackingKalman
from truth_tracking_gsf import runTruthTrackingGsf

from helpers import failure_threshold


def test_truth_tracking_kalman(physmon: "Physmon"):
    s = acts.examples.Sequencer(
        events=100,
        numThreads=-1,
        logLevel=acts.logging.INFO,
        fpeMasks=acts.examples.Sequencer.FpeMask.fromFile(
            Path(__file__).parent.parent / "fpe_masks.yml"
        ),
    )

    runTruthTrackingKalman(
        physmon.trackingGeometry,
        physmon.field,
        digiConfigFile=physmon.digiConfig,
        outputDir=physmon.tmp_path,
        s=s,
    )

    s.run()
    del s

    physmon.add_output_file(
        "performance_track_fitter.root", rename="performance_truth_tracking.root"
    )
    physmon.histogram_comparison(
        "performance_truth_tracking.root", title="Truth tracking KF"
    )


def test_truth_tracking_gsf(physmon: "Physmon"):
    s = acts.examples.Sequencer(
        events=500,
        numThreads=-1,
        logLevel=acts.logging.INFO,
        fpeMasks=acts.examples.Sequencer.FpeMask.fromFile(
            Path(__file__).parent.parent / "fpe_masks.yml"
        ),
    )

    runTruthTrackingGsf(
        physmon.trackingGeometry,
        physmon.digiConfig,
        physmon.field,
        outputDir=physmon.tmp_path,
        s=s,
    )

    with failure_threshold(acts.logging.FATAL):
        s.run()
    del s

    physmon.add_output_file("performance_gsf.root", rename="performance_gsf.root")
    physmon.histogram_comparison("performance_gsf.root", title="Truth tracking GSF")
