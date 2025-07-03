import os
import inspect
from pathlib import Path
import shutil
import math
import sys
import tempfile

import pytest

from helpers import (
    dd4hepEnabled,
    hepmc3Enabled,
    geant4Enabled,
    AssertCollectionExistsAlg,
)

import acts
from acts import UnitConstants as u
from acts.examples import (
    ObjPropagationStepsWriter,
    TrackFinderNTupleWriter,
    RootPropagationStepsWriter,
    RootParticleWriter,
    RootTrackParameterWriter,
    RootMaterialTrackWriter,
    RootMaterialWriter,
    RootSimHitWriter,
    RootTrackStatesWriter,
    RootTrackSummaryWriter,
    VertexNTupleWriter,
    RootMeasurementWriter,
    CsvParticleWriter,
    CsvSimHitWriter,
    CsvTrackWriter,
    CsvTrackingGeometryWriter,
    CsvMeasurementWriter,
    JsonMaterialWriter,
    JsonFormat,
    Sequencer,
    GenericDetector,
)
from acts.examples.odd import getOpenDataDetectorDirectory


@pytest.mark.obj
def test_obj_propagation_step_writer(tmp_path, trk_geo, conf_const, basic_prop_seq):
    with pytest.raises(TypeError):
        ObjPropagationStepsWriter()

    obj = tmp_path / "obj"
    obj.mkdir()

    s, alg = basic_prop_seq(trk_geo)
    w = conf_const(
        ObjPropagationStepsWriter,
        acts.logging.INFO,
        collection=alg.config.outputSummaryCollection,
        outputDir=str(obj),
    )

    s.addWriter(w)

    s.run()

    assert len([f for f in obj.iterdir() if f.is_file()]) == s.config.events
    for f in obj.iterdir():
        assert f.stat().st_size > 1024


@pytest.mark.csv
def test_csv_particle_writer(tmp_path, conf_const, ptcl_gun):
    s = Sequencer(numThreads=1, events=10)
    _, h3conv = ptcl_gun(s)

    out = tmp_path / "csv"

    out.mkdir()

    s.addWriter(
        conf_const(
            CsvParticleWriter,
            acts.logging.INFO,
            inputParticles=h3conv.config.outputParticles,
            outputStem="particle",
            outputDir=str(out),
        )
    )

    s.run()

    assert len([f for f in out.iterdir() if f.is_file()]) == s.config.events
    assert all(f.stat().st_size > 200 for f in out.iterdir())


@pytest.mark.root
def test_root_prop_step_writer(
    tmp_path, trk_geo, conf_const, basic_prop_seq, assert_root_hash
):
    with pytest.raises(TypeError):
        RootPropagationStepsWriter()

    file = tmp_path / "prop_steps.root"
    assert not file.exists()

    s, alg = basic_prop_seq(trk_geo)
    w = conf_const(
        RootPropagationStepsWriter,
        acts.logging.INFO,
        collection=alg.config.outputSummaryCollection,
        filePath=str(file),
    )

    s.addWriter(w)

    s.run()

    assert file.exists()
    assert file.stat().st_size > 2**10 * 50
    assert_root_hash(file.name, file)


@pytest.mark.root
def test_root_particle_writer(tmp_path, conf_const, ptcl_gun, assert_root_hash):
    s = Sequencer(numThreads=1, events=10)
    _, h3conv = ptcl_gun(s)

    file = tmp_path / "particles.root"

    assert not file.exists()

    s.addWriter(
        conf_const(
            RootParticleWriter,
            acts.logging.INFO,
            inputParticles=h3conv.config.outputParticles,
            filePath=str(file),
        )
    )

    s.run()

    assert file.exists()
    assert file.stat().st_size > 1024 * 10
    assert_root_hash(file.name, file)


@pytest.mark.root
def test_root_meas_writer(tmp_path, fatras, trk_geo, assert_root_hash):
    s = Sequencer(numThreads=1, events=10)
    evGen, simAlg, digiAlg = fatras(s)

    out = tmp_path / "meas.root"

    assert not out.exists()

    config = RootMeasurementWriter.Config(
        inputMeasurements=digiAlg.config.outputMeasurements,
        inputClusters=digiAlg.config.outputClusters,
        inputSimHits=simAlg.config.outputSimHits,
        inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
        filePath=str(out),
        surfaceByIdentifier=trk_geo.geoIdSurfaceMap(),
    )
    s.addWriter(RootMeasurementWriter(level=acts.logging.INFO, config=config))
    s.run()

    assert out.exists()
    assert out.stat().st_size > 40000
    assert_root_hash(out.name, out)


@pytest.mark.root
def test_root_simhits_writer(tmp_path, fatras, conf_const, assert_root_hash):
    s = Sequencer(numThreads=1, events=10)
    evGen, simAlg, digiAlg = fatras(s)

    out = tmp_path / "meas.root"

    assert not out.exists()

    s.addWriter(
        conf_const(
            RootSimHitWriter,
            level=acts.logging.INFO,
            inputSimHits=simAlg.config.outputSimHits,
            filePath=str(out),
        )
    )

    s.run()
    assert out.exists()
    assert out.stat().st_size > 2e4
    assert_root_hash(out.name, out)


@pytest.mark.root
def test_root_tracksummary_writer(tmp_path, fatras, conf_const):
    detector = GenericDetector()
    trackingGeometry = detector.trackingGeometry()
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    s = Sequencer(numThreads=1, events=10)

    from truth_tracking_kalman import runTruthTrackingKalman

    # This also runs the RootTrackSummaryWriter with truth information
    runTruthTrackingKalman(
        trackingGeometry,
        field,
        digiConfigFile=Path(
            str(
                Path(__file__).parent.parent.parent.parent
                / "Examples/Configs/generic-digi-smearing-config.json"
            )
        ),
        outputDir=tmp_path,
        s=s,
    )

    # Run the RootTrackSummaryWriter without the truth information
    s.addWriter(
        conf_const(
            RootTrackSummaryWriter,
            level=acts.logging.INFO,
            inputTracks="tracks",
            filePath=str(tmp_path / "track_summary_kf_no_truth.root"),
        )
    )

    s.run()
    assert (tmp_path / "tracksummary_kf.root").exists()
    assert (tmp_path / "track_summary_kf_no_truth.root").exists()


@pytest.mark.csv
def test_csv_meas_writer(tmp_path, fatras, trk_geo, conf_const):
    s = Sequencer(numThreads=1, events=10)
    evGen, simAlg, digiAlg = fatras(s)

    out = tmp_path / "csv"
    out.mkdir()

    s.addWriter(
        conf_const(
            CsvMeasurementWriter,
            level=acts.logging.INFO,
            inputMeasurements=digiAlg.config.outputMeasurements,
            inputClusters=digiAlg.config.outputClusters,
            inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
            outputDir=str(out),
        )
    )
    s.run()

    assert len([f for f in out.iterdir() if f.is_file()]) == s.config.events * 3
    assert all(f.stat().st_size > 10 for f in out.iterdir())


@pytest.mark.csv
def test_csv_simhits_writer(tmp_path, fatras, conf_const):
    s = Sequencer(numThreads=1, events=10)
    evGen, simAlg, digiAlg = fatras(s)

    out = tmp_path / "csv"
    out.mkdir()

    s.addWriter(
        conf_const(
            CsvSimHitWriter,
            level=acts.logging.INFO,
            inputSimHits=simAlg.config.outputSimHits,
            outputDir=str(out),
            outputStem="hits",
        )
    )

    s.run()
    assert len([f for f in out.iterdir() if f.is_file()]) == s.config.events
    assert all(f.stat().st_size > 200 for f in out.iterdir())


@pytest.mark.parametrize(
    "writer",
    [
        RootPropagationStepsWriter,
        RootParticleWriter,
        TrackFinderNTupleWriter,
        RootTrackParameterWriter,
        RootMaterialTrackWriter,
        RootMeasurementWriter,
        RootMaterialWriter,
        RootSimHitWriter,
        RootTrackStatesWriter,
        RootTrackSummaryWriter,
        VertexNTupleWriter,
    ],
)
@pytest.mark.root
def test_root_writer_interface(writer, conf_const, tmp_path, trk_geo):
    assert hasattr(writer, "Config")

    config = writer.Config

    assert hasattr(config, "filePath")
    assert hasattr(config, "fileMode")

    f = tmp_path / "target.root"
    assert not f.exists()

    kw = {"level": acts.logging.INFO, "filePath": str(f)}

    for k, _ in inspect.getmembers(config):
        if k.startswith("input"):
            kw[k] = "collection"
        if k == "surfaceByIdentifier":
            kw[k] = trk_geo.geoIdSurfaceMap()

    assert conf_const(writer, **kw)

    assert f.exists()


@pytest.mark.parametrize(
    "writer",
    [
        CsvParticleWriter,
        CsvMeasurementWriter,
        CsvSimHitWriter,
        CsvTrackWriter,
        CsvTrackingGeometryWriter,
    ],
)
@pytest.mark.csv
def test_csv_writer_interface(writer, conf_const, tmp_path, trk_geo):
    assert hasattr(writer, "Config")

    config = writer.Config

    assert hasattr(config, "outputDir")

    kw = {"level": acts.logging.INFO, "outputDir": str(tmp_path)}

    for k, _ in inspect.getmembers(config):
        if k.startswith("input"):
            kw[k] = "collection"
        if k == "trackingGeometry":
            kw[k] = trk_geo
        if k == "outputStem":
            kw[k] = "stem"

    assert conf_const(writer, **kw)


@pytest.mark.root
@pytest.mark.odd
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_root_material_writer(tmp_path, assert_root_hash):
    from acts.examples.dd4hep import DD4hepDetector

    detector = DD4hepDetector(
        xmlFileNames=[str(getOpenDataDetectorDirectory() / "xml/OpenDataDetector.xml")]
    )
    trackingGeometry = detector.trackingGeometry()

    out = tmp_path / "material.root"

    assert not out.exists()

    rmw = RootMaterialWriter(level=acts.logging.WARNING, filePath=str(out))
    assert out.exists()
    assert out.stat().st_size > 0 and out.stat().st_size < 500
    rmw.write(trackingGeometry)

    assert out.stat().st_size > 1000
    assert_root_hash(out.name, out)


@pytest.mark.json
@pytest.mark.odd
@pytest.mark.parametrize("fmt", [JsonFormat.Json, JsonFormat.Cbor])
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_json_material_writer(tmp_path, fmt):
    from acts.examples.dd4hep import DD4hepDetector

    detector = DD4hepDetector(
        xmlFileNames=[str(getOpenDataDetectorDirectory() / "xml/OpenDataDetector.xml")]
    )
    trackingGeometry = detector.trackingGeometry()

    out = (tmp_path / "material").with_suffix("." + fmt.name.lower())

    assert not out.exists()

    jmw = JsonMaterialWriter(
        level=acts.logging.WARNING, fileName=str(out.with_suffix("")), writeFormat=fmt
    )
    assert not out.exists()
    jmw.write(trackingGeometry)

    assert out.stat().st_size > 1000


@pytest.mark.csv
def test_csv_multitrajectory_writer(tmp_path):
    detector = GenericDetector()
    trackingGeometry = detector.trackingGeometry()
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    from truth_tracking_kalman import runTruthTrackingKalman

    s = Sequencer(numThreads=1, events=10)
    runTruthTrackingKalman(
        trackingGeometry,
        field,
        digiConfigFile=Path(
            str(
                Path(__file__).parent.parent.parent.parent
                / "Examples/Configs/generic-digi-smearing-config.json"
            )
        ),
        outputDir=tmp_path,
        s=s,
    )

    csv_dir = tmp_path / "csv"
    csv_dir.mkdir()
    s.addWriter(
        CsvTrackWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputDir=str(csv_dir),
        )
    )
    s.run()
    assert len([f for f in csv_dir.iterdir() if f.is_file()]) == 10
    assert all(f.stat().st_size > 20 for f in csv_dir.iterdir())
