from pathlib import Path
import os
import json
import functools
import subprocess

import pytest

from helpers import (
    geant4Enabled,
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


def test_fatras(trk_geo, tmp_path, field):
    from fatras import runFatras

    csv = tmp_path / "csv"
    csv.mkdir()

    nevents = 10

    root_files = [
        ("fatras_particles_final.root", "particles", nevents),
        ("fatras_particles_initial.root", "particles", nevents),
        ("hits.root", "hits", 115),
    ]

    assert len(list(csv.iterdir())) == 0
    for rf, _, _ in root_files:
        assert not (tmp_path / rf).exists()

    seq = Sequencer(events=nevents)
    runFatras(trk_geo, field, str(tmp_path), s=seq).run()

    del seq

    assert_csv_output(csv, "particles_final")
    assert_csv_output(csv, "particles_initial")
    assert_csv_output(csv, "hits")
    for f, tn, exp_entries in root_files:
        rfp = tmp_path / f
        assert rfp.exists()
        assert rfp.stat().st_size > 2 ** 10 * 10

        assert_entries(rfp, tn, exp_entries)


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


@pytest.mark.slow
@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_material_recording(tmp_path, material_recording):

    # Not quite sure why this isn't 200
    root_files = [("geant4_material_tracks.root", "material-tracks", 198)]

    for fn, tn, ee in root_files:
        fp = material_recording / fn
        assert fp.exists()
        assert fp.stat().st_size > 2 ** 10 * 50
        assert_entries(fp, tn, ee)


@pytest.mark.slow
@pytest.mark.skipif(not hepmc3Enabled, reason="HepMC3 plugin not available")
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
def test_event_recording(tmp_path):

    script = (
        Path(__file__).parent.parent.parent.parent
        / "Examples"
        / "Scripts"
        / "Python"
        / "event_recording.py"
    )
    assert script.exists()

    env = os.environ.copy()
    env["NEVENTS"] = "1"
    subprocess.check_call([str(script)], cwd=tmp_path, env=env)

    from acts.examples.hepmc3 import HepMC3AsciiReader

    out_path = tmp_path / "hepmc3"
    # out_path.mkdir()

    assert len([f for f in out_path.iterdir() if f.name.endswith("events.hepmc3")]) > 0
    assert all([f.stat().st_size > 100 for f in out_path.iterdir()])

    s = Sequencer(numThreads=1)

    s.addReader(
        HepMC3AsciiReader(
            level=acts.logging.INFO,
            inputDir=str(out_path),
            inputStem="events",
            outputEvents="hepmc-events",
        )
    )

    alg = AssertCollectionExistsAlg(
        "hepmc-events", name="check_alg", level=acts.logging.INFO
    )
    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 1


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


@pytest.mark.slow
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_material_mapping(material_recording, tmp_path):
    map_file = tmp_path / "material-maps_tracks.root"
    assert not map_file.exists()

    s = Sequencer(numThreads=1)

    detector, trackingGeometry, decorators = getOpenDataDetector()

    from material_mapping import runMaterialMapping

    runMaterialMapping(
        trackingGeometry,
        decorators,
        outputDir=str(tmp_path),
        inputDir=material_recording,
        s=s,
    )

    s.run()

    # MaterialMapping alg only writes on destruct.
    # See https://github.com/acts-project/acts/issues/881
    del s

    mat_file = tmp_path / "material-map.json"

    assert mat_file.exists()
    assert mat_file.stat().st_size > 10

    with mat_file.open() as fh:
        assert json.load(fh)

    assert map_file.exists()
    assert_entries(map_file, "material-tracks", 198)

    val_file = tmp_path / "propagation-material.root"
    assert not val_file.exists()

    # test the validation as well

    # we need to destroy the ODD to reload with material
    # del trackingGeometry
    # del detector

    detector, trackingGeometry, decorators = getOpenDataDetector(
        mdecorator=acts.IMaterialDecorator.fromFile(mat_file)
    )

    from material_validation import runMaterialValidation

    s = Sequencer(events=10, numThreads=1)

    field = acts.NullBField()

    runMaterialValidation(
        trackingGeometry, decorators, field, outputDir=str(tmp_path), s=s
    )

    s.run()

    assert val_file.exists()
    assert_entries(val_file, "material-tracks", 10000)


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
