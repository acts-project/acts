from pathlib import Path
import os
import json
import functools
import subprocess

import pytest

from helpers import (
    rootEnabled,
    dd4hepEnabled,
    hepmc3Enabled,
    AssertCollectionExistsAlg,
    isCI,
)

pytestmark = pytest.mark.skipif(not rootEnabled, reason="ROOT not set up")


import acts
from acts.examples import (
    Sequencer,
    GenericDetector,
    AlignedDetector,
    RootParticleWriter,
)

from common import getOpenDataDetector

u = acts.UnitConstants


@pytest.fixture
def field():
    return acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))


@pytest.fixture
def seq():
    return Sequencer(events=10, numThreads=1)


def assert_csv_output(csv_path, stem):
    __tracebackhide__ = True
    # print(list(csv_path.iterdir()))
    assert len([f for f in csv_path.iterdir() if f.name.endswith(stem + ".csv")]) > 0
    assert all([f.stat().st_size > 100 for f in csv_path.iterdir()])


def assert_entries(root_file, tree_name, exp):
    __tracebackhide__ = True
    import ROOT

    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)

    rf = ROOT.TFile.Open(str(root_file))
    keys = [k.GetName() for k in rf.GetListOfKeys()]
    assert tree_name in keys
    assert rf.Get(tree_name).GetEntries() == exp, f"{root_file}:{tree_name}"


def test_propagation(tmp_path, trk_geo, field, seq):
    from propagation import runPropagation

    obj = tmp_path / "obj"
    obj.mkdir()

    root_files = [("propagation_steps.root", "propagation_steps", 10000)]

    for fn, _, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    assert len(list(obj.iterdir())) == 0

    runPropagation(trk_geo, field, str(tmp_path), s=seq).run()

    for fn, tn, ee in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 2 ** 10 * 50
        assert_entries(fp, tn, ee)

    assert len(list(obj.iterdir())) > 0


def test_particle_gun(tmp_path):
    from particle_gun import runParticleGun

    s = Sequencer(events=20, numThreads=1)

    csv_dir = tmp_path / "csv"
    root_file = tmp_path / "particles.root"

    assert not csv_dir.exists()
    assert not root_file.exists()

    runParticleGun(str(tmp_path), s=s).run()

    assert csv_dir.exists()
    assert root_file.exists()

    assert len([f for f in csv_dir.iterdir() if f.name.endswith("particles.csv")]) > 0
    assert all([f.stat().st_size > 100 for f in csv_dir.iterdir()])

    assert root_file.stat().st_size > 200
    assert_entries(root_file, "particles", 20)


@pytest.mark.parametrize(
    "geoFactory,nobj",
    [
        (GenericDetector.create, 450),
        pytest.param(
            getOpenDataDetector,
            540,
            marks=pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up"),
        ),
        (functools.partial(AlignedDetector.create, iovSize=1), 450),
    ],
)
def test_geometry_example(geoFactory, nobj, tmp_path):
    detector, trackingGeometry, decorators = geoFactory()

    from geometry import runGeometry

    json_dir = tmp_path / "json"
    csv_dir = tmp_path / "csv"
    obj_dir = tmp_path / "obj"

    for d in (json_dir, csv_dir, obj_dir):
        d.mkdir()

    events = 5

    kwargs = dict(
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        events=events,
        outputDir=str(tmp_path),
    )

    runGeometry(outputJson=True, **kwargs)
    runGeometry(outputJson=False, **kwargs)

    assert len(list(obj_dir.iterdir())) == nobj
    assert all(f.stat().st_size > 200 for f in obj_dir.iterdir())

    assert len(list(csv_dir.iterdir())) == 3 * events
    assert all(f.stat().st_size > 200 for f in csv_dir.iterdir())

    detector_files = [csv_dir / f"event{i:>09}-detectors.csv" for i in range(events)]
    for detector_file in detector_files:
        assert detector_file.exists()
        assert detector_file.stat().st_size > 200

    contents = [f.read_text() for f in detector_files]
    ref = contents[0]
    for c in contents[1:]:
        if isinstance(detector, AlignedDetector):
            assert c != ref, "Detector writeout is expected to be different"
        else:
            assert c == ref, "Detector writeout is expected to be identical"

    if not isinstance(detector, AlignedDetector):
        for f in [json_dir / f"event{i:>09}-detector.json" for i in range(events)]:
            assert detector_file.exists()
            with f.open() as fh:
                data = json.load(fh)
                assert data
        material_file = tmp_path / "geometry-map.json"
        assert material_file.exists()
        assert material_file.stat().st_size > 200


def test_ckf_tracks_example_full_seeding(tmp_path):
    # the example as written is only compatible with the generic detector
    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    events = 10
    s = Sequencer(events=events, numThreads=1)  # Digitization is not thread-safe

    root_files = [
        ("performance_ckf.root", None, None),
        ("performance_seeding_trees.root", "track_finder_tracks", 368),
        ("performance_seeding_trees.root", "track_finder_particles", 80),
        ("trackstates_ckf.root", "trackstates", 368),
        ("tracksummary_ckf.root", "tracksummary", 10),
    ]

    csv = tmp_path / "csv"

    assert not csv.exists()
    for rf, _, _ in root_files:
        assert not (tmp_path / rf).exists()

    from ckf_tracks import runCKFTracks

    runCKFTracks(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=Path(
            "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
        ),
        digiConfigFile=Path(
            "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
        ),
        outputCsv=True,
        outputDir=tmp_path,
        truthSmearedSeeded=False,
        truthEstimatedSeeded=False,
        s=s,
    )
    s.run()

    del s  # files are closed in destructors, not great

    assert csv.exists()
    for rf, tn, nume in root_files:
        rp = tmp_path / rf
        assert rp.exists()
        if tn is not None and nume is not None:
            assert_entries(rp, tn, nume)

    assert len([f for f in csv.iterdir() if f.name.endswith("CKFtracks.csv")]) == events
    assert all([f.stat().st_size > 300 for f in csv.iterdir()])


def test_ckf_tracks_example_truth_estimate(tmp_path):
    # the example as written is only compatible with the generic detector
    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    events = 10
    s = Sequencer(events=events, numThreads=1)  # Digitization is not thread-safe

    root_files = [
        ("performance_ckf.root", None, None),
        ("performance_seeding_trees.root", "track_finder_tracks", 80),
        ("performance_seeding_trees.root", "track_finder_particles", 80),
        ("trackstates_ckf.root", "trackstates", 80),
        ("tracksummary_ckf.root", "tracksummary", 10),
    ]

    csv = tmp_path / "csv"

    assert not csv.exists()
    for rf, _, _ in root_files:
        assert not (tmp_path / rf).exists()

    from ckf_tracks import runCKFTracks

    runCKFTracks(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=Path(
            "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
        ),
        digiConfigFile=Path(
            "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
        ),
        outputCsv=True,
        outputDir=tmp_path,
        truthSmearedSeeded=False,
        truthEstimatedSeeded=True,
        s=s,
    )
    s.run()

    del s  # files are closed in destructors, not great

    assert csv.exists()
    for rf, tn, nume in root_files:
        rp = tmp_path / rf
        assert rp.exists()
        if tn is not None and nume is not None:
            assert_entries(rp, tn, nume)

    assert len([f for f in csv.iterdir() if f.name.endswith("CKFtracks.csv")]) == events
    assert all([f.stat().st_size > 100 for f in csv.iterdir()])


def test_ckf_tracks_example_truth_smeared(tmp_path):
    # the example as written is only compatible with the generic detector
    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    events = 10
    s = Sequencer(events=events, numThreads=1)  # Digitization is not thread-safe

    root_files = [
        ("performance_ckf.root", None, None),
        ("trackstates_ckf.root", "trackstates", 80),
        ("tracksummary_ckf.root", "tracksummary", 10),
    ]

    csv = tmp_path / "csv"

    assert not csv.exists()
    for rf, _, _ in root_files:
        assert not (tmp_path / rf).exists()

    from ckf_tracks import runCKFTracks

    runCKFTracks(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=Path(
            "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
        ),
        digiConfigFile=Path(
            "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
        ),
        outputCsv=True,
        outputDir=tmp_path,
        truthSmearedSeeded=True,
        truthEstimatedSeeded=False,
        s=s,
    )
    s.run()

    del s  # files are closed in destructors, not great

    assert csv.exists()
    for rf, tn, nume in root_files:
        rp = tmp_path / rf
        assert rp.exists()
        if tn is not None and nume is not None:
            assert_entries(rp, tn, nume)

    assert len([f for f in csv.iterdir() if f.name.endswith("CKFtracks.csv")]) == events
    assert all([f.stat().st_size > 300 for f in csv.iterdir()])
