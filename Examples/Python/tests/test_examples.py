from pathlib import Path
import os
import json
import functools
import tarfile
import urllib.request
import subprocess
import sys
import re
import collections

import pytest

from helpers import (
    geant4Enabled,
    rootEnabled,
    dd4hepEnabled,
    hepmc3Enabled,
    pythia8Enabled,
    exatrkxEnabled,
    onnxEnabled,
    AssertCollectionExistsAlg,
    isCI,
    doHashChecks,
    failure_threshold,
)

pytestmark = pytest.mark.skipif(not rootEnabled, reason="ROOT not set up")


import acts
from acts.examples import (
    Sequencer,
    GenericDetector,
    AlignedDetector,
    RootParticleWriter,
)

from acts.examples.odd import getOpenDataDetector
from common import getOpenDataDetectorDirectory

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


def assert_entries(root_file, tree_name, exp=None, non_zero=False):
    __tracebackhide__ = True
    import ROOT

    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)

    rf = ROOT.TFile.Open(str(root_file))
    keys = [k.GetName() for k in rf.GetListOfKeys()]
    assert tree_name in keys
    print("Entries:", rf.Get(tree_name).GetEntries())
    if non_zero:
        assert rf.Get(tree_name).GetEntries() > 0, f"{root_file}:{tree_name}"
    if exp is not None:
        assert rf.Get(tree_name).GetEntries() == exp, f"{root_file}:{tree_name}"


def assert_has_entries(root_file, tree_name):
    __tracebackhide__ = True
    assert_entries(root_file, tree_name, non_zero=True)


@pytest.mark.slow
@pytest.mark.skipif(not pythia8Enabled, reason="Pythia8 not set up")
def test_pythia8(tmp_path, seq, assert_root_hash):
    from pythia8 import runPythia8

    (tmp_path / "csv").mkdir()

    assert not (tmp_path / "pythia8_particles.root").exists()
    assert len(list((tmp_path / "csv").iterdir())) == 0

    events = seq.config.events

    runPythia8(str(tmp_path), outputRoot=True, outputCsv=True, s=seq).run()

    del seq

    fp = tmp_path / "pythia8_particles.root"
    assert fp.exists()
    assert fp.stat().st_size > 2**10 * 50
    assert_entries(fp, "particles", events)
    assert_root_hash(fp.name, fp)

    assert len(list((tmp_path / "csv").iterdir())) > 0
    assert_csv_output(tmp_path / "csv", "particles")


def test_fatras(trk_geo, tmp_path, field, assert_root_hash):
    from fatras import runFatras

    csv = tmp_path / "csv"
    csv.mkdir()

    nevents = 10

    root_files = [
        (
            "particles_final.root",
            "particles",
        ),
        (
            "particles_initial.root",
            "particles",
        ),
        (
            "hits.root",
            "hits",
        ),
    ]

    assert len(list(csv.iterdir())) == 0
    for rf, _ in root_files:
        assert not (tmp_path / rf).exists()

    seq = Sequencer(events=nevents)
    runFatras(trk_geo, field, str(tmp_path), s=seq).run()

    del seq

    assert_csv_output(csv, "particles_final")
    assert_csv_output(csv, "particles_initial")
    assert_csv_output(csv, "hits")
    for f, tn in root_files:
        rfp = tmp_path / f
        assert rfp.exists()
        assert rfp.stat().st_size > 2**10 * 10

        assert_has_entries(rfp, tn)
        assert_root_hash(f, rfp)


@pytest.mark.slow
@pytest.mark.odd
@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_geant4(tmp_path, assert_root_hash):
    # This test literally only ensures that the geant 4 example can run without erroring out
    getOpenDataDetector(
        getOpenDataDetectorDirectory()
    )  # just to make sure it can build

    csv = tmp_path / "csv"
    csv.mkdir()

    root_files = [
        "particles_final.root",
        "particles_initial.root",
        "hits.root",
    ]

    assert len(list(csv.iterdir())) == 0
    for rf in root_files:
        assert not (tmp_path / rf).exists()

    script = (
        Path(__file__).parent.parent.parent.parent
        / "Examples"
        / "Scripts"
        / "Python"
        / "geant4.py"
    )
    assert script.exists()
    env = os.environ.copy()
    env["ACTS_LOG_FAILURE_THRESHOLD"] = "WARNING"
    try:
        subprocess.check_call(
            [sys.executable, str(script)],
            cwd=tmp_path,
            env=env,
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as e:
        print(e.output.decode("utf-8"))
        raise

    assert_csv_output(csv, "particles_final")
    assert_csv_output(csv, "particles_initial")
    assert_csv_output(csv, "hits")
    for f in root_files:
        rfp = tmp_path / f
        assert rfp.exists()
        assert rfp.stat().st_size > 2**10 * 10

        assert_root_hash(f, rfp)


def test_seeding(tmp_path, trk_geo, field, assert_root_hash):
    from seeding import runSeeding

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))

    csv = tmp_path / "csv"
    csv.mkdir()

    seq = Sequencer(events=10, numThreads=1)

    root_files = [
        (
            "estimatedparams.root",
            "estimatedparams",
        ),
        (
            "performance_seeding.root",
            None,
        ),
        (
            "particles.root",
            "particles",
        ),
        (
            "particles_final.root",
            "particles",
        ),
        (
            "particles_initial.root",
            "particles",
        ),
    ]

    for fn, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    assert len(list(csv.iterdir())) == 0

    runSeeding(trk_geo, field, outputDir=str(tmp_path), s=seq).run()

    del seq

    for fn, tn in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 100

        if tn is not None:
            assert_has_entries(fp, tn)
            assert_root_hash(fn, fp)

    assert_csv_output(csv, "particles")
    assert_csv_output(csv, "particles_final")
    assert_csv_output(csv, "particles_initial")


def test_seeding_orthogonal(tmp_path, trk_geo, field, assert_root_hash):
    from seeding import runSeeding, SeedingAlgorithm

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))

    csv = tmp_path / "csv"
    csv.mkdir()

    seq = Sequencer(events=10, numThreads=1)

    root_files = [
        (
            "estimatedparams.root",
            "estimatedparams",
        ),
        (
            "performance_seeding.root",
            None,
        ),
        (
            "particles.root",
            "particles",
        ),
        (
            "particles_final.root",
            "particles",
        ),
        (
            "particles_initial.root",
            "particles",
        ),
    ]

    for fn, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    assert len(list(csv.iterdir())) == 0

    runSeeding(
        trk_geo,
        field,
        outputDir=str(tmp_path),
        s=seq,
        seedingAlgorithm=SeedingAlgorithm.Orthogonal,
    ).run()

    del seq

    for fn, tn in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 100

        if tn is not None:
            assert_has_entries(fp, tn)
            assert_root_hash(fn, fp)

    assert_csv_output(csv, "particles")
    assert_csv_output(csv, "particles_final")
    assert_csv_output(csv, "particles_initial")


def test_itk_seeding(tmp_path, trk_geo, field, assert_root_hash):
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))

    csv = tmp_path / "csv"
    csv.mkdir()

    seq = Sequencer(events=10, numThreads=1)

    root_files = [
        (
            "estimatedparams.root",
            "estimatedparams",
        ),
        (
            "performance_seeding.root",
            None,
        ),
        (
            "particles.root",
            "particles",
        ),
        (
            "particles_final.root",
            "particles",
        ),
        (
            "particles_initial.root",
            "particles",
        ),
    ]

    for fn, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    assert len(list(csv.iterdir())) == 0

    rnd = acts.examples.RandomNumbers(seed=42)

    from acts.examples.simulation import (
        addParticleGun,
        EtaConfig,
        MomentumConfig,
        ParticleConfig,
        addFatras,
        addDigitization,
    )

    addParticleGun(
        seq,
        MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, True),
        EtaConfig(-4.0, 4.0, True),
        ParticleConfig(1, acts.PdgParticle.eMuon, True),
        outputDirCsv=tmp_path / "csv",
        outputDirRoot=str(tmp_path),
        rnd=rnd,
    )

    addFatras(
        seq,
        trk_geo,
        field,
        outputDirCsv=tmp_path / "csv",
        outputDirRoot=str(tmp_path),
        rnd=rnd,
    )

    srcdir = Path(__file__).resolve().parent.parent.parent.parent
    addDigitization(
        seq,
        trk_geo,
        field,
        digiConfigFile=srcdir
        / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json",
        rnd=rnd,
    )

    from acts.examples.reconstruction import (
        addSeeding,
        TruthSeedRanges,
    )
    from acts.examples.reconstruction import (
        addSeeding,
        TruthSeedRanges,
    )
    from acts.examples.itk import itkSeedingAlgConfig, InputSpacePointsType

    addSeeding(
        seq,
        trk_geo,
        field,
        TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-4, 4), nHits=(9, None)),
        *itkSeedingAlgConfig(InputSpacePointsType.PixelSpacePoints),
        acts.logging.VERBOSE,
        geoSelectionConfigFile=srcdir
        / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json",
        inputParticles="particles_final",  # use this to reproduce the original root_file_hashes.txt - remove to fix
        outputDirRoot=str(tmp_path),
    )

    seq.run()

    del seq

    for fn, tn in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 100

        if tn is not None:
            assert_has_entries(fp, tn)
            assert_root_hash(fn, fp)

    assert_csv_output(csv, "particles")
    assert_csv_output(csv, "particles_final")
    assert_csv_output(csv, "particles_initial")


@pytest.mark.slow
def test_propagation(tmp_path, trk_geo, field, seq, assert_root_hash):
    from propagation import runPropagation

    obj = tmp_path / "obj"
    obj.mkdir()

    root_files = [
        (
            "propagation_steps.root",
            "propagation_steps",
            10000,
        )
    ]

    for fn, _, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    assert len(list(obj.iterdir())) == 0

    runPropagation(trk_geo, field, str(tmp_path), s=seq).run()

    for fn, tn, ee in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 2**10 * 50
        assert_entries(fp, tn, ee)
        assert_root_hash(fn, fp)

    assert len(list(obj.iterdir())) > 0


@pytest.mark.slow
@pytest.mark.odd
@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_material_recording(tmp_path, material_recording, assert_root_hash):
    root_files = [
        (
            "geant4_material_tracks.root",
            "material-tracks",
            200,
        )
    ]

    for fn, tn, ee in root_files:
        fp = material_recording / fn
        assert fp.exists()
        assert fp.stat().st_size > 2**10 * 50
        assert_entries(fp, tn, ee)
        assert_root_hash(fn, fp)


@pytest.mark.slow
@pytest.mark.odd
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
    env["ACTS_LOG_FAILURE_THRESHOLD"] = "WARNING"
    try:
        subprocess.check_call(
            [sys.executable, str(script)],
            cwd=tmp_path,
            env=env,
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as e:
        print(e.output.decode("utf-8"))
        raise

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


@pytest.mark.parametrize("revFiltMomThresh", [0 * u.GeV, 1 * u.TeV])
def test_truth_tracking_kalman(
    tmp_path, assert_root_hash, revFiltMomThresh, detector_config
):
    from truth_tracking_kalman import runTruthTrackingKalman

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    seq = Sequencer(events=10, numThreads=1)

    root_files = [
        ("trackstates_fitter.root", "trackstates", 19),
        ("tracksummary_fitter.root", "tracksummary", 10),
        ("performance_track_finder.root", "track_finder_tracks", 19),
        ("performance_track_fitter.root", None, -1),
    ]

    for fn, _, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    runTruthTrackingKalman(
        trackingGeometry=detector_config.trackingGeometry,
        field=field,
        digiConfigFile=detector_config.digiConfigFile,
        outputDir=tmp_path,
        reverseFilteringMomThreshold=revFiltMomThresh,
        s=seq,
    )

    seq.run()

    del seq

    for fn, tn, ee in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 1024
        if tn is not None:
            assert_has_entries(fp, tn)
            assert_root_hash(fn, fp)


def test_truth_tracking_gsf(tmp_path, assert_root_hash, detector_config):
    from truth_tracking_gsf import runTruthTrackingGsf

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    seq = Sequencer(
        events=10,
        numThreads=1,
        fpeMasks=[
            (
                "Fatras/include/ActsFatras/Kernel/detail/SimulationActor.hpp:177",
                acts.FpeType.FLTUND,
                1,
            ),
            (
                "Core/include/Acts/TrackFitting/detail/GsfUtils.hpp:197",
                acts.FpeType.FLTUND,
                1,
            ),
        ],
    )

    root_files = [
        ("trackstates_gsf.root", "trackstates"),
        ("tracksummary_gsf.root", "tracksummary"),
    ]

    for fn, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    runTruthTrackingGsf(
        trackingGeometry=detector_config.trackingGeometry,
        decorators=detector_config.decorators,
        field=field,
        digiConfigFile=detector_config.digiConfigFile,
        outputDir=tmp_path,
        s=seq,
    )

    # See https://github.com/acts-project/acts/issues/1300
    with failure_threshold(acts.logging.FATAL):
        seq.run()

    del seq

    for fn, tn in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 1024
        if tn is not None:
            assert_root_hash(fn, fp)


def test_particle_gun(tmp_path, assert_root_hash):
    from particle_gun import runParticleGun

    s = Sequencer(events=20, numThreads=-1)

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
    assert_root_hash(root_file.name, root_file)


@pytest.mark.slow
@pytest.mark.odd
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_material_mapping(material_recording, tmp_path, assert_root_hash):
    map_file = tmp_path / "material-map_tracks.root"
    assert not map_file.exists()

    s = Sequencer(numThreads=1)

    detector, trackingGeometry, decorators = getOpenDataDetector(
        getOpenDataDetectorDirectory()
    )

    from material_mapping import runMaterialMapping

    runMaterialMapping(
        trackingGeometry,
        decorators,
        outputDir=str(tmp_path),
        inputDir=material_recording,
        mappingStep=1,
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
    assert_entries(map_file, "material-tracks", 200)
    assert_root_hash(map_file.name, map_file)

    val_file = tmp_path / "propagation-material.root"
    assert not val_file.exists()

    # test the validation as well

    # we need to destroy the ODD to reload with material
    del trackingGeometry
    del detector

    detector, trackingGeometry, decorators = getOpenDataDetector(
        getOpenDataDetectorDirectory(),
        mdecorator=acts.IMaterialDecorator.fromFile(mat_file),
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
    assert_root_hash(val_file.name, val_file)


@pytest.mark.slow
@pytest.mark.odd
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_volume_material_mapping(material_recording, tmp_path, assert_root_hash):
    map_file = tmp_path / "material-map-volume_tracks.root"
    assert not map_file.exists()

    s = Sequencer(numThreads=1)

    geo_map = Path(__file__).parent / "geometry-volume-map.json"
    assert geo_map.exists()
    assert geo_map.stat().st_size > 10
    with geo_map.open() as fh:
        assert json.load(fh)

    detector, trackingGeometry, decorators = getOpenDataDetector(
        getOpenDataDetectorDirectory(),
        mdecorator=acts.IMaterialDecorator.fromFile(geo_map),
    )

    from material_mapping import runMaterialMapping

    runMaterialMapping(
        trackingGeometry,
        decorators,
        mapName="material-map-volume",
        outputDir=str(tmp_path),
        inputDir=material_recording,
        mappingStep=1,
        s=s,
    )

    s.run()

    # MaterialMapping alg only writes on destruct.
    # See https://github.com/acts-project/acts/issues/881
    del s

    mat_file = tmp_path / "material-map-volume.json"

    assert mat_file.exists()
    assert mat_file.stat().st_size > 10

    with mat_file.open() as fh:
        assert json.load(fh)

    assert map_file.exists()
    assert_entries(map_file, "material-tracks", 200)
    assert_root_hash(map_file.name, map_file)

    val_file = tmp_path / "propagation-volume-material.root"
    assert not val_file.exists()

    # test the validation as well

    # we need to destroy the ODD to reload with material
    del trackingGeometry
    del detector

    detector, trackingGeometry, decorators = getOpenDataDetector(
        getOpenDataDetectorDirectory(),
        mdecorator=acts.IMaterialDecorator.fromFile(mat_file),
    )

    from material_validation import runMaterialValidation

    s = Sequencer(events=10, numThreads=1)

    field = acts.NullBField()

    runMaterialValidation(
        trackingGeometry,
        decorators,
        field,
        outputDir=str(tmp_path),
        outputName="propagation-volume-material",
        s=s,
    )

    s.run()

    assert val_file.exists()
    assert_root_hash(val_file.name, val_file)


@pytest.mark.parametrize(
    "geoFactory,nobj",
    [
        (GenericDetector.create, 450),
        pytest.param(
            functools.partial(getOpenDataDetector, getOpenDataDetectorDirectory()),
            540,
            marks=[
                pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up"),
                pytest.mark.slow,
                pytest.mark.odd,
            ],
        ),
        (functools.partial(AlignedDetector.create, iovSize=1), 450),
    ],
)
@pytest.mark.slow
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


DIGI_SHARE_DIR = (
    Path(__file__).parent.parent.parent.parent
    / "Examples/Algorithms/Digitization/share"
)


@pytest.mark.parametrize(
    "digi_config_file",
    [
        DIGI_SHARE_DIR / "default-smearing-config-generic.json",
        DIGI_SHARE_DIR / "default-geometric-config-generic.json",
    ],
    ids=["smeared", "geometric"],
)
def test_digitization_example(trk_geo, tmp_path, assert_root_hash, digi_config_file):
    from digitization import runDigitization

    s = Sequencer(events=10, numThreads=-1)

    csv_dir = tmp_path / "csv"
    root_file = tmp_path / "measurements.root"

    assert not root_file.exists()
    assert not csv_dir.exists()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    runDigitization(
        trk_geo, field, outputDir=tmp_path, digiConfigFile=digi_config_file, s=s
    )

    s.run()

    assert root_file.exists()
    assert csv_dir.exists()

    assert len(list(csv_dir.iterdir())) == 3 * s.config.events
    assert all(f.stat().st_size > 50 for f in csv_dir.iterdir())

    assert_entries(root_file, "vol9", 0)
    assert_entries(root_file, "vol14", 0)

    if "smearing" in digi_config_file.name:
        filled_entries = [f"vol{tn}" for tn in (8, 12, 13, 16, 17, 18)]
    else:
        # fmt: off
        filled_entries = [
            'vol8', 'vol8_lay2', 'vol12_lay8_mod117', 'vol12_lay10', 'vol12_lay10_mod154',
            'vol12_lay10_mod163', 'vol12_lay12', 'vol12_lay12_mod150', 'vol13',
            'vol13_lay2', 'vol16_lay2_mod53', 'vol16_lay4', 'vol16_lay6', 'vol16_lay8',
            'vol16_lay10', 'vol16_lay12', 'vol17', 'vol17_lay2', 'vol18_lay2',
            'vol18_lay2_mod1', 'vol18_lay2_mod49', 'vol18_lay2_mod86', 'vol18_lay4',
        ]
        # fmt: on

    for entry in filled_entries:
        assert_has_entries(root_file, entry)

    assert_root_hash(root_file.name, root_file)


@pytest.mark.parametrize(
    "digi_config_file",
    [
        DIGI_SHARE_DIR / "default-smearing-config-generic.json",
        DIGI_SHARE_DIR / "default-geometric-config-generic.json",
    ],
    ids=["smeared", "geometric"],
)
def test_digitization_example_input(
    trk_geo, tmp_path, assert_root_hash, digi_config_file
):
    from particle_gun import runParticleGun
    from digitization import runDigitization

    ptcl_dir = tmp_path / "ptcl"
    ptcl_dir.mkdir()
    pgs = Sequencer(events=20, numThreads=-1)
    runParticleGun(str(ptcl_dir), s=pgs)
    pgs.run()

    s = Sequencer(numThreads=-1)

    csv_dir = tmp_path / "csv"
    root_file = tmp_path / "measurements.root"

    assert not root_file.exists()
    assert not csv_dir.exists()

    assert_root_hash(
        "particles.root",
        ptcl_dir / "particles.root",
    )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    runDigitization(
        trk_geo,
        field,
        outputDir=tmp_path,
        digiConfigFile=digi_config_file,
        particlesInput=ptcl_dir / "particles.root",
        s=s,
        doMerge=True,
    )

    s.run()

    assert root_file.exists()
    assert csv_dir.exists()

    assert len(list(csv_dir.iterdir())) == 3 * pgs.config.events
    assert all(f.stat().st_size > 50 for f in csv_dir.iterdir())

    assert_entries(root_file, "vol7", 0)
    assert_entries(root_file, "vol9", 0)

    if "smearing" in digi_config_file.name:
        filled_entries = [f"vol{tn}" for tn in (8, 12, 13, 16, 17, 18)]
    else:
        # fmt: off
        filled_entries = [
            "vol8", "vol8_lay2", "vol12_lay8_mod120", "vol12_lay10_mod120",
            "vol12_lay10_mod144", "vol12_lay12", "vol12_lay12_mod111",
            "vol12_lay12_mod137", "vol12_lay12_mod170", "vol13", "vol13_lay2",
            "vol14_lay2_mod93", "vol14_lay2_mod102", "vol14_lay2_mod112",
            "vol14_lay2_mod118", "vol14_lay4_mod112", "vol14_lay4_mod118",
            "vol14_lay4_mod152", "vol14_lay4_mod161", "vol16_lay4", "vol16_lay6",
            "vol16_lay8", "vol16_lay10", "vol16_lay12", "vol17", "vol17_lay2",
            "vol18_lay2", "vol18_lay2_mod71", "vol18_lay4", "vol18_lay6",
            "vol18_lay8", "vol18_lay10"
        ]
        # fmt: on

    for entry in filled_entries:
        assert_has_entries(root_file, entry)

    assert_root_hash(root_file.name, root_file)


def test_digitization_config_example(trk_geo, tmp_path):
    from digitization_config import runDigitizationConfig

    out_file = tmp_path / "output.json"
    assert not out_file.exists()

    input = (
        Path(__file__).parent
        / "../../../Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
    )
    assert input.exists(), input.resolve()

    runDigitizationConfig(trk_geo, input=input, output=out_file)

    assert out_file.exists()

    with out_file.open() as fh:
        data = json.load(fh)
    assert len(data.keys()) == 2
    assert data["acts-geometry-hierarchy-map"]["format-version"] == 0
    assert (
        data["acts-geometry-hierarchy-map"]["value-identifier"]
        == "digitization-configuration"
    )
    assert len(data["entries"]) == 27


@pytest.mark.parametrize(
    "truthSmeared,truthEstimated",
    [
        [False, False],
        [False, True],
        [True, False],
    ],
    ids=["full_seeding", "truth_estimated", "truth_smeared"],
)
@pytest.mark.slow
def test_ckf_tracks_example(
    tmp_path, assert_root_hash, truthSmeared, truthEstimated, detector_config
):
    csv = tmp_path / "csv"

    assert not csv.exists()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    events = 100
    s = Sequencer(events=events, numThreads=1)  # Digitization is not thread-safe

    root_files = [
        (
            "performance_ckf.root",
            None,
        ),
        (
            "trackstates_ckf.root",
            "trackstates",
        ),
        (
            "tracksummary_ckf.root",
            "tracksummary",
        ),
    ]

    if not truthSmeared:
        root_files += [
            (
                "performance_seeding.root",
                None,
            ),
        ]

    for rf, _ in root_files:
        assert not (tmp_path / rf).exists()

    from ckf_tracks import runCKFTracks

    runCKFTracks(
        detector_config.trackingGeometry,
        detector_config.decorators,
        field=field,
        outputCsv=True,
        outputDir=tmp_path,
        geometrySelection=detector_config.geometrySelection,
        digiConfigFile=detector_config.digiConfigFile,
        truthSmearedSeeded=truthSmeared,
        truthEstimatedSeeded=truthEstimated,
        s=s,
    )

    s.run()

    del s  # files are closed in destructors, not great

    assert csv.exists()
    for rf, tn in root_files:
        rp = tmp_path / rf
        assert rp.exists()
        if tn is not None:
            assert_root_hash(rf, rp)

    assert (
        len([f for f in csv.iterdir() if f.name.endswith("tracks_ckf.csv")]) == events
    )
    assert all([f.stat().st_size > 300 for f in csv.iterdir()])


@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
@pytest.mark.odd
@pytest.mark.slow
def test_full_chain_odd_example(tmp_path):
    # This test literally only ensures that the full chain example can run without erroring out
    getOpenDataDetector(
        getOpenDataDetectorDirectory()
    )  # just to make sure it can build

    script = (
        Path(__file__).parent.parent.parent.parent
        / "Examples"
        / "Scripts"
        / "Python"
        / "full_chain_odd.py"
    )
    assert script.exists()
    env = os.environ.copy()
    env["ACTS_LOG_FAILURE_THRESHOLD"] = "ERROR"
    try:
        subprocess.check_call(
            [sys.executable, str(script), "-n1"],
            cwd=tmp_path,
            env=env,
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as e:
        print(e.output.decode("utf-8"))
        raise


@pytest.mark.skipif(
    not dd4hepEnabled or not geant4Enabled, reason="DD4hep and/or Geant4 not set up"
)
@pytest.mark.slow
def test_full_chain_odd_example_pythia_geant4(tmp_path):
    # This test literally only ensures that the full chain example can run without erroring out
    getOpenDataDetector(
        getOpenDataDetectorDirectory()
    )  # just to make sure it can build

    script = (
        Path(__file__).parent.parent.parent.parent
        / "Examples"
        / "Scripts"
        / "Python"
        / "full_chain_odd.py"
    )
    assert script.exists()
    env = os.environ.copy()
    env["ACTS_LOG_FAILURE_THRESHOLD"] = "ERROR"
    try:
        stdout = subprocess.check_output(
            [sys.executable, str(script), "-n1", "--geant4", "--ttbar"],
            cwd=tmp_path,
            env=env,
            stderr=subprocess.STDOUT,
        )
        stdout = stdout.decode("utf-8")
    except subprocess.CalledProcessError as e:
        print(e.output.decode("utf-8"))
        raise

    # collect and compare known errors
    errors = []
    error_regex = re.compile(r"^\d\d:\d\d:\d\d\s+(\w+)\s+ERROR\s+", re.MULTILINE)
    for match in error_regex.finditer(stdout):
        (algo,) = match.groups()
        errors.append(algo)
    errors = collections.Counter(errors)
    assert dict(errors) == {}, stdout


@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
@pytest.mark.skipif(not onnxEnabled, reason="ONNX plugin not enabled")
@pytest.mark.slow
def test_ML_Ambiguity_Solver(tmp_path, assert_root_hash):
    root_file = "performance_ambiML.root"
    output_dir = "odd_output"
    assert not (tmp_path / root_file).exists()
    # This test literally only ensures that the full chain example can run without erroring out
    getOpenDataDetector(
        getOpenDataDetectorDirectory()
    )  # just to make sure it can build

    script = (
        Path(__file__).parent.parent.parent.parent
        / "Examples"
        / "Scripts"
        / "Python"
        / "full_chain_odd.py"
    )
    assert script.exists()
    env = os.environ.copy()
    env["ACTS_LOG_FAILURE_THRESHOLD"] = "ERROR"
    try:
        subprocess.check_call(
            [sys.executable, str(script), "-n5", "--MLSolver"],
            cwd=tmp_path,
            env=env,
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as e:
        print(e.output.decode("utf-8"))
        raise

    rfp = tmp_path / output_dir / root_file
    assert rfp.exists()

    assert_root_hash(root_file, rfp)


def test_bfield_writing(tmp_path, seq, assert_root_hash):
    from bfield_writing import runBFieldWriting

    root_files = [
        ("solenoid.root", "solenoid", 100),
        ("solenoid2.root", "solenoid", 100),
    ]

    for fn, _, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    runBFieldWriting(outputDir=tmp_path, rewrites=1)

    for fn, tn, ee in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 2**10 * 2
        assert_entries(fp, tn, ee)
        assert_root_hash(fn, fp)


@pytest.mark.parametrize("backend", ["onnx", "torch"])
@pytest.mark.parametrize("hardware", ["cpu", "gpu"])
@pytest.mark.skipif(not exatrkxEnabled, reason="ExaTrkX environment not set up")
def test_exatrkx(tmp_path, trk_geo, field, assert_root_hash, backend, hardware):
    if backend == "onnx" and hardware == "cpu":
        pytest.skip("Combination of ONNX and CPU not yet supported")

    root_file = "performance_track_finding.root"
    assert not (tmp_path / root_file).exists()

    if backend == "onnx":
        url = "https://acts.web.cern.ch/ci/exatrkx/onnx_models_v01.tar"
    else:
        url = "https://acts.web.cern.ch/ci/exatrkx/torchscript_models_v01.tar"

    tarfile_name = tmp_path / "models.tar"
    urllib.request.urlretrieve(url, tarfile_name)
    tarfile.open(tarfile_name).extractall(tmp_path)
    script = (
        Path(__file__).parent.parent.parent.parent
        / "Examples"
        / "Scripts"
        / "Python"
        / "exatrkx.py"
    )
    assert script.exists()
    env = os.environ.copy()
    env["ACTS_LOG_FAILURE_THRESHOLD"] = "WARNING"

    if hardware == "cpu":
        env["CUDA_VISIBLE_DEVICES"] = ""

    try:
        subprocess.check_call(
            [sys.executable, str(script), backend],
            cwd=tmp_path,
            env=env,
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as e:
        print(e.output.decode("utf-8"))
        raise

    rfp = tmp_path / root_file
    assert rfp.exists()

    assert_root_hash(root_file, rfp)
