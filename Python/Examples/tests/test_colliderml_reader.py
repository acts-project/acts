import os
import subprocess
import sys
from pathlib import Path

import filelock
import pytest

import acts
import acts.examples
from acts import UnitConstants as u

from helpers import arrowEnabled, dd4hepEnabled

pytestmark = [
    pytest.mark.skipif(not arrowEnabled, reason="Arrow/Parquet bindings not built"),
    pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up"),
]

_srcdir = Path(__file__).resolve().parent.parent.parent.parent

# Maximum number of events to process in the test — keep small so it runs fast
# even on the full (1000-event) downloaded shard.
_N_EVENTS = 3


@pytest.fixture(scope="session")
def colliderml_pu0_sample(tmp_path_factory):
    """Locate or download a ttbar PU0 ColliderML sample.

    If COLLIDERML_DATA_DIR is set, use that directory directly (expected
    to contain CERN__ColliderML-Release-1/). Otherwise download via the
    colliderml CLI (requires network access and uv in PATH).

    Returns the path to the CERN__ColliderML-Release-1 directory.
    """
    data_dir_env = os.environ.get("COLLIDERML_DATA_DIR")
    if data_dir_env:
        base = Path(data_dir_env)
    else:
        base = tmp_path_factory.getbasetemp().parent / "colliderml_pu0"
        with filelock.FileLock(str(base) + ".lock"):
            release_dir = base / "CERN__ColliderML-Release-1"
            if not (release_dir / "ttbar_pu0_particles").exists():
                base.mkdir(parents=True, exist_ok=True)
                # Clear PYTHONHOME/PYTHONPATH so uv can manage its own
                # interpreter without interference from LCG or other venvs.
                dl_env = os.environ.copy()
                dl_env.pop("PYTHONHOME", None)
                dl_env.pop("PYTHONPATH", None)
                subprocess.check_call(
                    [
                        "uv",
                        "run",
                        "--with",
                        "colliderml",
                        "--no-project",
                        "colliderml",
                        "download",
                        "--pileup",
                        "pu0",
                        "--channels",
                        "ttbar",
                        "--objects",
                        "particles,tracker_hits",
                        "--max-events",
                        str(_N_EVENTS),
                        "--out",
                        str(base),
                    ],
                    env=dl_env,
                )

    release_dir = base / "CERN__ColliderML-Release-1"
    if not release_dir.exists():
        pytest.skip(f"ColliderML PU0 data not found at {release_dir}")
    return release_dir


def test_colliderml_truth_tracking_pu0(
    tmp_path, assert_root_hash, colliderml_pu0_sample
):
    """Read ttbar PU0 ColliderML data and run truth-seeded Kalman tracking.

    Verifies that the ColliderML reader pipeline (ParquetReader +
    ColliderMLInputConverter + TruthEstimated seeding + KF) runs without
    error and produces ROOT output files with the expected content hashes.
    """
    from acts.examples.odd import getOpenDataDetector

    from colliderml_truth_tracking import runColliderMLTruthTracking

    odd_dir = acts.examples.odd.getOpenDataDetectorDirectory()
    matDeco = acts.IMaterialDecorator.fromFile(
        odd_dir / "data/odd-material-maps.root", level=acts.logging.WARNING
    )
    detector = getOpenDataDetector(matDeco)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    with detector:
        s, _perf_proto, _perf_kf = runColliderMLTruthTracking(
            trackingGeometry=trackingGeometry,
            field=field,
            outputDir=tmp_path,
            inputDir=colliderml_pu0_sample,
            geoIdMapPath=_srcdir / "Examples/Configs/colliderml_geo_map.parquet",
            digiConfigFile=_srcdir
            / "Examples/Configs/odd-digi-smearing-config-notime.json",
            decorators=decorators,
            events=_N_EVENTS,
            numThreads=1,
            sample="ttbar_pu0",
        )
        s.run()

    root_files = [
        "trackstates_kf.root",
        "tracksummary_kf.root",
        "performance_kf.root",
    ]

    for fn in root_files:
        fp = tmp_path / fn
        assert fp.exists(), f"{fn} was not produced"
        assert fp.stat().st_size > 1024, f"{fn} is suspiciously small"
        assert_root_hash(fn, fp)
