import contextlib
import sys
from pathlib import Path

import pytest

import acts
import acts.examples
from acts import UnitConstants as u
from acts.examples.dataset import getColliderMLObjectDirectory

from helpers import arrowEnabled, dd4hepEnabled, AssertCollectionExistsAlg

_srcdir = Path(__file__).resolve().parent.parent.parent.parent

sys.path.insert(0, str(_srcdir / "Examples/Scripts/Python"))
from colliderml_truth_tracking import runColliderMLTruthTracking

_CHANNEL = "ttbar"
_PILEUP = "pu0"

try:
    for _object in ("particles", "tracker_hits", "tracks"):
        getColliderMLObjectDirectory(_object, channel=_CHANNEL, pileup=_PILEUP)
    _colliderml_data_available = True
except RuntimeError:
    _colliderml_data_available = False

pytestmark = [
    pytest.mark.skipif(not arrowEnabled, reason="Arrow/Parquet bindings not built"),
    pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up"),
    pytest.mark.skipif(
        not _colliderml_data_available,
        reason=(
            f"ColliderML CI sample not found; set COLLIDERML_DATA_DIR to a "
            f"directory containing the '{_CHANNEL}_{_PILEUP}' sample"
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

    The ROOT performance check now runs for gen3 too, but there is no stored
    reference hash for it yet in root_file_hashes.txt -- this is expected to
    fail with "Hash lookup key ... not found" until someone with access to
    the ColliderML dataset runs it, confirms the tracking performance output
    is sane, and adds the reference. Left failing on purpose rather than
    skipped, so this isn't silently forgotten.
    """
    from acts.examples.odd import getOpenDataDetector
    from generate_geoid_map import generate_geoid_map

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    if reco_geo == "gen3":
        # Build both geometries to generate the mapping on the fly
        detector_gen1 = getOpenDataDetector(gen3=False)
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
        detector = getOpenDataDetector(matDeco, gen3=False)
        trackingGeometry = detector.trackingGeometry()
        decorators = detector.contextDecorators()
        ctx = detector

    check_alg = AssertCollectionExistsAlg(
        "colliderml_tracks", name="check_colliderml_tracks", level=acts.logging.WARNING
    )

    perfFile = tmp_path / "performance_colliderml_truth_tracking.root"
    RootPatternRecognitionPerformanceWriter = acts.examples._tryImportRoot(
        "RootPatternRecognitionPerformanceWriter"
    )
    perfWriter = RootPatternRecognitionPerformanceWriter(
        level=acts.logging.WARNING,
        inputTracks="tracks",
        inputParticles="particles",
        inputTrackParticleMatching="track_particle_matching",
        inputParticleTrackMatching="particle_track_matching",
        inputParticleMeasurementsMap="particle_measurements_map",
        filePath=str(perfFile),
    )

    with ctx if ctx is not None else contextlib.nullcontext():
        s = runColliderMLTruthTracking(
            trackingGeometry=trackingGeometry,
            field=field,
            outputDir=tmp_path,
            channel=_CHANNEL,
            pileup=_PILEUP,
            geoIdMapPath=geoid_map_path if reco_geo == "gen3" else None,
            decorators=decorators,
            numThreads=1,
        )
        s.addAlgorithm(check_alg)
        s.addWriter(perfWriter)
        s.run()

    assert check_alg.events_seen > 0

    assert_root_hash("performance_colliderml_truth_tracking.root", perfFile)
