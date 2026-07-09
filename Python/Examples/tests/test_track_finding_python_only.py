from pathlib import Path

import pytest

import acts
from acts import UnitConstants as u
from acts.examples import Sequencer

_srcdir = Path(__file__).resolve().parent.parent.parent.parent


@pytest.mark.pypi
def test_track_finding_python_only(tmp_path, generic_detector_config):
    """Pytest wrapper for track_finding_python_only.py, the script used as
    the cibuildwheel wheel smoke test. Exercises particle gun -> Fatras ->
    digitization -> pure-Python proto-tracking/fitting -> truth matching
    -> performance histograms, using only what's in the wheel (no ROOT).
    """
    with generic_detector_config.detector:
        from track_finding_python_only import runTrackFindingPythonOnly

        field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
        geoSelectionConfigFile = (
            _srcdir / "Examples/Configs/generic-pixel-sstrips-lstrips-spacepoints.json"
        )

        s, perfWriter = runTrackFindingPythonOnly(
            trackingGeometry=generic_detector_config.trackingGeometry,
            field=field,
            digiConfigFile=generic_detector_config.digiConfigFile,
            geoSelectionConfigFile=geoSelectionConfigFile,
            outputDir=tmp_path,
            decorators=generic_detector_config.decorators,
            s=Sequencer(events=1, numThreads=1, logLevel=acts.logging.INFO),
        )
        s.run()

        histograms = perfWriter.histograms()
        assert len(histograms) > 0, "no performance histograms were produced"
