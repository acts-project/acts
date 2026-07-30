import os
import sys
from pathlib import Path

import pytest

import acts
import acts.examples
from acts import UnitConstants as u

from helpers import arrowEnabled, dd4hepEnabled, AssertCollectionExistsAlg

_srcdir = Path(__file__).resolve().parent.parent.parent.parent

sys.path.insert(0, str(_srcdir / "Examples/Scripts/Python"))
from colliderml_truth_tracking import runColliderMLTruthTracking

_COLLIDERML_DATA_ENV_VAR = "COLLIDERML_DATA"
_SAMPLE = "ttbar_pu0"

_colliderml_data_dir = os.environ.get(_COLLIDERML_DATA_ENV_VAR)
if _colliderml_data_dir is not None:
    _dataDir = Path(_colliderml_data_dir)
    _particlesDir = _dataDir / f"{_SAMPLE}_particles" / "data" / f"{_SAMPLE}_particles"
    _hitsDir = _dataDir / f"{_SAMPLE}_tracker_hits" / "data" / f"{_SAMPLE}_tracker_hits"
    _tracksDir = _dataDir / f"{_SAMPLE}_tracks" / "data" / f"{_SAMPLE}_tracks"
    _colliderml_data_available = (
        _particlesDir.exists() and _hitsDir.exists() and _tracksDir.exists()
    )
else:
    _particlesDir = _hitsDir = _tracksDir = None
    _colliderml_data_available = False

pytestmark = [
    pytest.mark.skipif(not arrowEnabled, reason="Arrow/Parquet bindings not built"),
    pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up"),
    pytest.mark.skipif(
        not _colliderml_data_available,
        reason=(
            f"ColliderML CI sample not found; set {_COLLIDERML_DATA_ENV_VAR} to a "
            f"directory containing the '{_SAMPLE}' sample"
        ),
    ),
]


@pytest.mark.root
@pytest.mark.parametrize("reco_geo", ["gen1", "gen3"])
def test_colliderml_truth_tracking(tmp_path, reco_geo, assert_root_hash):
    """Read a small real ColliderML (ttbar, PU0) sample and run truth-seeded
    KF tracking on it, for both gen1 (no geo-id map needed) and gen3 (geo-id
    map generated on the fly) reconstruction geometries. Also checks that the
    dataset's own published tracks are read and converted.

    gen3's seeding config hardcodes gen1 volume/layer IDs, so it currently
    produces zero seeds on real data; the ROOT performance check is skipped
    for gen3 until that's fixed.
    """
    from acts.examples.odd import getOpenDataDetector
    from generate_geoid_map import generate_geoid_map

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    if reco_geo == "gen3":
        # Build both geometries to generate the mapping on the fly
        detector_gen1 = getOpenDataDetector()
        geo_gen1 = detector_gen1.trackingGeometry()

        detector_gen3 = getOpenDataDetector(gen3=True)
        trackingGeometry = detector_gen3.trackingGeometry()
        decorators = detector_gen3.contextDecorators()

        geoid_map_path = tmp_path / "geoid_map.csv"
        generate_geoid_map(
            geo_gen1,
            trackingGeometry,
            output_path=geoid_map_path,
            prefix_a="gen1",
            prefix_b="gen3",
        )

        ctx = None
    else:
        odd_dir = acts.examples.odd.getOpenDataDetectorDirectory()
        matDeco = acts.IMaterialDecorator.fromFile(
            odd_dir / "data/odd-material-maps.root", level=acts.logging.WARNING
        )
        detector = getOpenDataDetector(matDeco)
        trackingGeometry = detector.trackingGeometry()
        decorators = detector.contextDecorators()
        ctx = detector

    check_alg = AssertCollectionExistsAlg(
        "colliderml_tracks", name="check_colliderml_tracks", level=acts.logging.WARNING
    )

    perfFile = tmp_path / "performance_colliderml_truth_tracking.root"
    perfWriter = None
    if reco_geo == "gen1":
        RootTrackFinderPerformanceWriter = acts.examples._tryImportRoot(
            "RootTrackFinderPerformanceWriter"
        )
        perfWriter = RootTrackFinderPerformanceWriter(
            level=acts.logging.WARNING,
            inputTracks="tracks",
            inputParticles="particles",
            inputTrackParticleMatching="track_particle_matching",
            inputParticleTrackMatching="particle_track_matching",
            inputParticleMeasurementsMap="particle_measurements_map",
            filePath=str(perfFile),
        )

    def _run():
        s = runColliderMLTruthTracking(
            trackingGeometry=trackingGeometry,
            field=field,
            outputDir=tmp_path,
            particlesDir=_particlesDir,
            hitsDir=_hitsDir,
            tracksDir=_tracksDir,
            geoIdMapPath=geoid_map_path if reco_geo == "gen3" else None,
            decorators=decorators,
            numThreads=1,
        )
        s.addAlgorithm(check_alg)
        if perfWriter is not None:
            s.addWriter(perfWriter)
        s.run()

    if ctx is not None:
        with ctx:
            _run()
    else:
        _run()

    assert check_alg.events_seen > 0

    if reco_geo == "gen1":
        assert_root_hash("performance_colliderml_truth_tracking.root", perfFile)
