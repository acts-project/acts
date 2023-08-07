from pathlib import Path

import acts.examples

from truth_tracking_kalman import runTruthTrackingKalman


def test_truth_tracking_kalman(output_path: Path, tmp_path: Path, setup):
    print(setup)
    s = acts.examples.Sequencer(
        events=100,
        numThreads=-1,
        logLevel=acts.logging.INFO,
        fpeMasks=acts.examples.Sequencer.FpeMask.fromFile(
            Path(__file__).parent.parent / "fpe_masks.yml"
        ),
    )

    runTruthTrackingKalman(
        setup.trackingGeometry,
        setup.field,
        digiConfigFile=setup.digiConfig,
        outputDir=tmp_path,
        s=s,
    )

    s.run()
    del s

    perf_file = tmp_path / "performance_track_fitter.root"
    assert perf_file.exists(), "Performance file not found"
    #  shutil.copy(perf_file, setup.outdir / "performance_truth_tracking.root")
