from typing import Type
import inspect

import pytest

from helpers import dd4hepEnabled


import acts

from acts import PlanarModuleStepper


from acts.examples import (
    ObjPropagationStepsWriter,
    TrackFinderPerformanceWriter,
    SeedingPerformanceWriter,
    RootPropagationStepsWriter,
    RootParticleWriter,
    RootTrackParameterWriter,
    RootMaterialTrackWriter,
    RootMaterialWriter,
    RootPlanarClusterWriter,
    RootSimHitWriter,
    RootTrajectoryStatesWriter,
    RootTrajectorySummaryWriter,
    RootVertexPerformanceWriter,
    RootMeasurementWriter,
    CsvParticleWriter,
    CsvPlanarClusterWriter,
    CsvSimHitWriter,
    CsvMultiTrajectoryWriter,
    CsvTrackingGeometryWriter,
    CsvMeasurementWriter,
    TrackParamsEstimationAlgorithm,
    PlanarSteppingAlgorithm,
    JsonMaterialWriter,
    JsonFormat,
    Sequencer,
)


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
        collection=alg.config.propagationStepCollection,
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
    evGen = ptcl_gun(s)

    out = tmp_path / "csv"

    out.mkdir()

    s.addWriter(
        conf_const(
            CsvParticleWriter,
            acts.logging.INFO,
            inputParticles=evGen.config.outputParticles,
            outputStem="particle",
            outputDir=str(out),
        )
    )

    s.run()

    assert len([f for f in out.iterdir() if f.is_file()]) == s.config.events
    assert all(f.stat().st_size > 200 for f in out.iterdir())


@pytest.mark.root
def test_root_prop_step_writer(tmp_path, trk_geo, conf_const, basic_prop_seq):
    with pytest.raises(TypeError):
        RootPropagationStepsWriter()

    file = tmp_path / "prop_steps.root"
    assert not file.exists()

    s, alg = basic_prop_seq(trk_geo)
    w = conf_const(
        RootPropagationStepsWriter,
        acts.logging.INFO,
        collection=alg.config.propagationStepCollection,
        filePath=str(file),
    )

    s.addWriter(w)

    s.run()

    assert file.exists()
    assert file.stat().st_size > 2 ** 10 * 50


@pytest.mark.root
def test_root_particle_writer(tmp_path, conf_const, ptcl_gun):
    s = Sequencer(numThreads=1, events=10)
    evGen = ptcl_gun(s)

    file = tmp_path / "particles.root"

    assert not file.exists()

    s.addWriter(
        conf_const(
            RootParticleWriter,
            acts.logging.INFO,
            inputParticles=evGen.config.outputParticles,
            filePath=str(file),
        )
    )

    s.run()

    assert file.exists()
    assert file.stat().st_size > 1024 * 10


@pytest.mark.root
def test_root_meas_writer(tmp_path, fatras, trk_geo):
    s = Sequencer(numThreads=1, events=10)
    evGen, simAlg, digiAlg = fatras(s)

    out = tmp_path / "meas.root"

    assert not out.exists()

    config = RootMeasurementWriter.Config(
        inputMeasurements=digiAlg.config.outputMeasurements,
        inputClusters=""
        if digiAlg.config.isSimpleSmearer
        else digiAlg.config.outputClusters,
        inputSimHits=simAlg.config.outputSimHits,
        inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
        filePath=str(out),
        trackingGeometry=trk_geo,
    )
    config.addBoundIndicesFromDigiConfig(digiAlg.config)
    s.addWriter(RootMeasurementWriter(level=acts.logging.INFO, config=config))
    s.run()

    assert out.exists()
    assert out.stat().st_size > 2 ** 10 * 50


@pytest.mark.root
def test_root_simhits_writer(tmp_path, fatras, conf_const):
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
    assert out.stat().st_size > 3e4


@pytest.mark.root
def test_root_clusters_writer(tmp_path, fatras, conf_const, trk_geo, rng):
    s = Sequencer(numThreads=1, events=10)  # we're not going to use this one
    evGen, simAlg, _ = fatras(s)
    s = Sequencer(numThreads=1, events=10)
    s.addReader(evGen)
    s.addAlgorithm(simAlg)
    digiAlg = PlanarSteppingAlgorithm(
        level=acts.logging.ERROR,
        inputSimHits=simAlg.config.outputSimHits,
        outputClusters="clusters",
        outputSourceLinks="sourcelinks",
        outputMeasurements="measurements",
        outputMeasurementParticlesMap="meas_ptcl_map",
        outputMeasurementSimHitsMap="meas_sh_map",
        trackingGeometry=trk_geo,
        randomNumbers=rng,
        planarModuleStepper=PlanarModuleStepper(),
    )
    s.addAlgorithm(digiAlg)

    out = tmp_path / "clusters.root"

    assert not out.exists()

    s.addWriter(
        conf_const(
            RootPlanarClusterWriter,
            level=acts.logging.ERROR,
            filePath=str(out),
            inputSimHits=simAlg.config.outputSimHits,
            inputClusters=digiAlg.config.outputClusters,
            trackingGeometry=trk_geo,
        )
    )

    s.run()
    assert out.exists()
    assert out.stat().st_size > 2 ** 10 * 50


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
            inputClusters=""
            if digiAlg.config.isSimpleSmearer
            else digiAlg.config.outputClusters,
            inputSimHits=simAlg.config.outputSimHits,
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


@pytest.mark.csv
def test_csv_clusters_writer(tmp_path, fatras, conf_const, trk_geo, rng):
    s = Sequencer(numThreads=1, events=10)  # we're not going to use this one
    evGen, simAlg, _ = fatras(s)
    s = Sequencer(numThreads=1, events=10)
    s.addReader(evGen)
    s.addAlgorithm(simAlg)
    digiAlg = PlanarSteppingAlgorithm(
        level=acts.logging.ERROR,
        inputSimHits=simAlg.config.outputSimHits,
        outputClusters="clusters",
        outputSourceLinks="sourcelinks",
        outputMeasurements="measurements",
        outputMeasurementParticlesMap="meas_ptcl_map",
        outputMeasurementSimHitsMap="meas_sh_map",
        trackingGeometry=trk_geo,
        randomNumbers=rng,
        planarModuleStepper=PlanarModuleStepper(),
    )
    s.addAlgorithm(digiAlg)

    out = tmp_path / "csv"
    out.mkdir()

    s.addWriter(
        conf_const(
            CsvPlanarClusterWriter,
            level=acts.logging.ERROR,
            outputDir=str(out),
            inputSimHits=simAlg.config.outputSimHits,
            inputClusters=digiAlg.config.outputClusters,
            trackingGeometry=trk_geo,
        )
    )

    s.run()
    assert len([f for f in out.iterdir() if f.is_file()]) == s.config.events * 3
    assert all(f.stat().st_size > 1024 for f in out.iterdir())


@pytest.mark.parametrize(
    "writer",
    [
        RootPropagationStepsWriter,
        RootParticleWriter,
        TrackFinderPerformanceWriter,
        SeedingPerformanceWriter,
        RootTrackParameterWriter,
        RootMaterialTrackWriter,
        RootMeasurementWriter,
        RootMaterialWriter,
        RootPlanarClusterWriter,
        RootSimHitWriter,
        RootTrajectoryStatesWriter,
        RootTrajectorySummaryWriter,
        RootVertexPerformanceWriter,
        SeedingPerformanceWriter,
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
        if k == "trackingGeometry":
            kw[k] = trk_geo

    assert conf_const(writer, **kw)

    assert f.exists()


@pytest.mark.parametrize(
    "writer",
    [
        CsvParticleWriter,
        CsvMeasurementWriter,
        CsvPlanarClusterWriter,
        CsvSimHitWriter,
        CsvMultiTrajectoryWriter,
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
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_root_material_writer(tmp_path):
    from acts.examples.dd4hep import DD4hepDetector

    detector, trackingGeometry, _ = DD4hepDetector.create(
        xmlFileNames=["thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
    )

    out = tmp_path / "material.root"

    assert not out.exists()

    rmw = RootMaterialWriter(level=acts.logging.ERROR, filePath=str(out))
    assert out.exists()
    assert out.stat().st_size > 0 and out.stat().st_size < 500
    rmw.write(trackingGeometry)

    assert out.stat().st_size > 1000


@pytest.mark.json
@pytest.mark.parametrize("fmt", [JsonFormat.Json, JsonFormat.Cbor])
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_json_material_writer(tmp_path, fmt):
    from acts.examples.dd4hep import DD4hepDetector

    detector, trackingGeometry, _ = DD4hepDetector.create(
        xmlFileNames=["thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
    )

    out = (tmp_path / "material").with_suffix("." + fmt.name.lower())

    assert not out.exists()

    jmw = JsonMaterialWriter(
        level=acts.logging.ERROR, fileName=str(out.with_suffix("")), writeFormat=fmt
    )
    assert not out.exists()
    jmw.write(trackingGeometry)

    assert out.stat().st_size > 1000
