from pathlib import Path

import acts.examples

from truth_tracking_kalman import runTruthTrackingKalman


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
