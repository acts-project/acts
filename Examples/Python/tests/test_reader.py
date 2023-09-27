import pytest
import os
import multiprocessing

from helpers import (
    geant4Enabled,
    edm4hepEnabled,
    AssertCollectionExistsAlg,
)
from common import getOpenDataDetectorDirectory

from acts.examples.odd import getOpenDataDetector

import acts
from acts import PlanarModuleStepper
from acts.examples import (
    RootParticleWriter,
    RootParticleReader,
    RootMaterialTrackReader,
    RootTrajectorySummaryReader,
    CsvParticleWriter,
    CsvParticleReader,
    CsvMeasurementWriter,
    CsvMeasurementReader,
    CsvSimHitWriter,
    CsvSimHitReader,
    CsvPlanarClusterWriter,
    CsvPlanarClusterReader,
    PlanarSteppingAlgorithm,
    Sequencer,
)


@pytest.mark.root
def test_root_particle_reader(tmp_path, conf_const, ptcl_gun):
    # need to write out some particles first
    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.WARNING)
    evGen = ptcl_gun(s)

    file = tmp_path / "particles.root"
    s.addWriter(
        conf_const(
            RootParticleWriter,
            acts.logging.WARNING,
            inputParticles=evGen.config.outputParticles,
            filePath=str(file),
        )
    )

    s.run()

    del s  # to properly close the root file

    # reset sequencer for reading

    s2 = Sequencer(numThreads=1, logLevel=acts.logging.WARNING)

    s2.addReader(
        conf_const(
            RootParticleReader,
            acts.logging.WARNING,
            particleCollection="input_particles",
            filePath=str(file),
        )
    )

    alg = AssertCollectionExistsAlg(
        "input_particles", "check_alg", acts.logging.WARNING
    )
    s2.addAlgorithm(alg)

    s2.run()

    assert alg.events_seen == 10


@pytest.mark.csv
def test_csv_particle_reader(tmp_path, conf_const, ptcl_gun):
    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.WARNING)
    evGen = ptcl_gun(s)

    out = tmp_path / "csv"

    out.mkdir()

    s.addWriter(
        conf_const(
            CsvParticleWriter,
            acts.logging.WARNING,
            inputParticles=evGen.config.outputParticles,
            outputStem="particle",
            outputDir=str(out),
        )
    )

    s.run()

    # reset the seeder
    s = Sequencer(numThreads=1, logLevel=acts.logging.WARNING)

    s.addReader(
        conf_const(
            CsvParticleReader,
            acts.logging.WARNING,
            inputDir=str(out),
            inputStem="particle",
            outputParticles="input_particles",
        )
    )

    alg = AssertCollectionExistsAlg(
        "input_particles", "check_alg", acts.logging.WARNING
    )

    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 10


@pytest.mark.parametrize(
    "reader",
    [RootParticleReader, RootTrajectorySummaryReader],
)
@pytest.mark.root
def test_root_reader_interface(reader, conf_const, tmp_path):
    assert hasattr(reader, "Config")

    config = reader.Config

    assert hasattr(config, "filePath")

    kw = {"level": acts.logging.INFO, "filePath": str(tmp_path / "file.root")}

    assert conf_const(reader, **kw)


@pytest.mark.slow
@pytest.mark.root
@pytest.mark.odd
@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
def test_root_material_track_reader(material_recording):

    input_tracks = material_recording / "geant4_material_tracks.root"
    assert input_tracks.exists()

    s = Sequencer(numThreads=1)

    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            fileList=[str(input_tracks)],
        )
    )

    alg = AssertCollectionExistsAlg(
        "material-tracks", "check_alg", acts.logging.WARNING
    )
    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 2


@pytest.mark.csv
def test_csv_meas_reader(tmp_path, fatras, trk_geo, conf_const):
    s = Sequencer(numThreads=1, events=10)
    evGen, simAlg, digiAlg = fatras(s)

    out = tmp_path / "csv"
    out.mkdir()

    config = CsvMeasurementWriter.Config(
        inputMeasurements=digiAlg.config.outputMeasurements,
        inputClusters=digiAlg.config.outputClusters,
        inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
        outputDir=str(out),
    )
    s.addWriter(CsvMeasurementWriter(level=acts.logging.INFO, config=config))
    s.run()

    # read back in
    s = Sequencer(numThreads=1)

    s.addReader(
        conf_const(
            CsvMeasurementReader,
            level=acts.logging.WARNING,
            outputMeasurements="measurements",
            outputMeasurementSimHitsMap="simhitsmap",
            outputSourceLinks="sourcelinks",
            inputDir=str(out),
        )
    )

    algs = [
        AssertCollectionExistsAlg(k, f"check_alg_{k}", acts.logging.WARNING)
        for k in ("measurements", "simhitsmap", "sourcelinks")
    ]
    for alg in algs:
        s.addAlgorithm(alg)

    s.run()

    for alg in algs:
        assert alg.events_seen == 10


@pytest.mark.csv
def test_csv_simhits_reader(tmp_path, fatras, conf_const):
    s = Sequencer(numThreads=1, events=10)
    evGen, simAlg, digiAlg = fatras(s)

    out = tmp_path / "csv"
    out.mkdir()

    s.addWriter(
        CsvSimHitWriter(
            level=acts.logging.INFO,
            inputSimHits=simAlg.config.outputSimHits,
            outputDir=str(out),
            outputStem="hits",
        )
    )

    s.run()

    s = Sequencer(numThreads=1)

    s.addReader(
        conf_const(
            CsvSimHitReader,
            level=acts.logging.INFO,
            inputDir=str(out),
            inputStem="hits",
            outputSimHits="simhits",
        )
    )

    alg = AssertCollectionExistsAlg("simhits", "check_alg", acts.logging.WARNING)
    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 10


@pytest.mark.csv
def test_csv_clusters_reader(tmp_path, fatras, conf_const, trk_geo, rng):
    s = Sequencer(numThreads=1, events=10)  # we're not going to use this one
    evGen, simAlg, _ = fatras(s)
    s = Sequencer(numThreads=1, events=10)
    s.addReader(evGen)
    s.addAlgorithm(simAlg)
    digiAlg = PlanarSteppingAlgorithm(
        level=acts.logging.WARNING,
        inputSimHits=simAlg.config.outputSimHits,
        outputClusters="clusters",
        outputSourceLinks="sourcelinks",
        outputDigiSourceLinks="digiSourceLink",
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
        CsvPlanarClusterWriter(
            level=acts.logging.WARNING,
            outputDir=str(out),
            inputSimHits=simAlg.config.outputSimHits,
            inputClusters=digiAlg.config.outputClusters,
            trackingGeometry=trk_geo,
        )
    )

    s.run()

    s = Sequencer(numThreads=1)

    s.addReader(
        conf_const(
            CsvPlanarClusterReader,
            level=acts.logging.WARNING,
            outputClusters="clusters",
            inputDir=str(out),
            outputHitIds="hits",
            outputMeasurementParticlesMap="meas_ptcl_map",
            outputSimHits="simhits",
            trackingGeometry=trk_geo,
        )
    )

    algs = [
        AssertCollectionExistsAlg(k, f"check_alg_{k}", acts.logging.WARNING)
        for k in ("clusters", "simhits", "meas_ptcl_map")
    ]
    for alg in algs:
        s.addAlgorithm(alg)

    s.run()

    for alg in algs:
        assert alg.events_seen == 10


def generate_input_test_edm4hep_simhit_reader(input, output):
    from DDSim.DD4hepSimulation import DD4hepSimulation

    ddsim = DD4hepSimulation()
    if isinstance(ddsim.compactFile, list):
        ddsim.compactFile = [input]
    else:
        ddsim.compactFile = input
    ddsim.enableGun = True
    ddsim.gun.direction = (1, 0, 0)
    ddsim.gun.distribution = "eta"
    ddsim.numberOfEvents = 10
    ddsim.outputFile = output
    ddsim.run()


@pytest.mark.slow
@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_simhit_reader(tmp_path):
    from acts.examples.edm4hep import EDM4hepSimHitReader

    tmp_file = str(tmp_path / "output_edm4hep.root")
    odd_xml_file = str(getOpenDataDetectorDirectory() / "xml" / "OpenDataDetector.xml")

    with multiprocessing.get_context("spawn").Pool() as pool:
        pool.apply(generate_input_test_edm4hep_simhit_reader, (odd_xml_file, tmp_file))

    assert os.path.exists(tmp_file)

    detector, trackingGeometry, decorators = getOpenDataDetector(
        getOpenDataDetectorDirectory()
    )

    s = Sequencer(numThreads=1)

    s.addReader(
        EDM4hepSimHitReader(
            level=acts.logging.INFO,
            inputPath=tmp_file,
            outputSimHits="simhits",
            dd4hepDetector=detector,
        )
    )

    alg = AssertCollectionExistsAlg("simhits", "check_alg", acts.logging.WARNING)
    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 10


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_measurement_reader(tmp_path, fatras, conf_const):
    from acts.examples.edm4hep import (
        EDM4hepMeasurementWriter,
        EDM4hepMeasurementReader,
    )

    s = Sequencer(numThreads=1, events=10)
    _, simAlg, digiAlg = fatras(s)

    out = tmp_path / "measurements_edm4hep.root"

    config = EDM4hepMeasurementWriter.Config(
        inputMeasurements=digiAlg.config.outputMeasurements,
        inputClusters=digiAlg.config.outputClusters,
        outputPath=str(out),
    )
    s.addWriter(EDM4hepMeasurementWriter(level=acts.logging.INFO, config=config))
    s.run()

    # read back in
    s = Sequencer(numThreads=1)

    s.addReader(
        conf_const(
            EDM4hepMeasurementReader,
            level=acts.logging.WARNING,
            outputMeasurements="measurements",
            outputMeasurementSimHitsMap="simhitsmap",
            outputSourceLinks="sourcelinks",
            inputPath=str(out),
        )
    )

    algs = [
        AssertCollectionExistsAlg(k, f"check_alg_{k}", acts.logging.WARNING)
        for k in ("measurements", "simhitsmap", "sourcelinks")
    ]
    for alg in algs:
        s.addAlgorithm(alg)

    s.run()

    for alg in algs:
        assert alg.events_seen == 10


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_particle_reader(tmp_path, conf_const, ptcl_gun):
    from acts.examples.edm4hep import (
        EDM4hepParticleWriter,
        EDM4hepParticleReader,
    )

    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.WARNING)
    evGen = ptcl_gun(s)

    out = tmp_path / "particles_edm4hep.root"

    out.mkdir()

    s.addWriter(
        conf_const(
            EDM4hepParticleWriter,
            acts.logging.WARNING,
            inputParticles=evGen.config.outputParticles,
            outputPath=str(out),
        )
    )

    s.run()

    # reset the seeder
    s = Sequencer(numThreads=1, logLevel=acts.logging.WARNING)

    s.addReader(
        conf_const(
            EDM4hepParticleReader,
            acts.logging.WARNING,
            inputPath=str(out),
            outputParticles="input_particles",
        )
    )

    alg = AssertCollectionExistsAlg(
        "input_particles", "check_alg", acts.logging.WARNING
    )

    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 10
