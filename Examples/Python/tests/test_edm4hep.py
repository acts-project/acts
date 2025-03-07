import os
import multiprocessing
from pathlib import Path

import pytest

from helpers import (
    edm4hepEnabled,
    podioEnabled,
    AssertCollectionExistsAlg,
)

import acts
from acts import UnitConstants as u

from acts.examples import (
    Sequencer,
    GenericDetector,
)

from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_measurement_writer(tmp_path, fatras):
    from acts.examples.edm4hep import EDM4hepMeasurementWriter

    s = Sequencer(numThreads=1, events=10)
    _, simAlg, digiAlg = fatras(s)

    out = tmp_path / "measurements_edm4hep.root"

    s.addWriter(
        EDM4hepMeasurementWriter(
            level=acts.logging.VERBOSE,
            inputMeasurements=digiAlg.config.outputMeasurements,
            inputClusters=digiAlg.config.outputClusters,
            outputPath=str(out),
        )
    )

    s.run()

    assert os.path.isfile(out)
    assert os.stat(out).st_size > 10


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_simhit_writer(tmp_path, fatras, conf_const):
    from acts.examples.edm4hep import EDM4hepSimHitWriter

    s = Sequencer(numThreads=1, events=10)
    _, simAlg, _ = fatras(s)

    out = tmp_path / "simhits_edm4hep.root"

    s.addWriter(
        conf_const(
            EDM4hepSimHitWriter,
            level=acts.logging.INFO,
            inputSimHits=simAlg.config.outputSimHits,
            outputPath=str(out),
        )
    )

    s.run()

    assert os.path.isfile(out)
    assert os.stat(out).st_size > 200


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_particle_writer(tmp_path, conf_const, ptcl_gun):
    from acts.examples.edm4hep import EDM4hepParticleWriter

    s = Sequencer(numThreads=1, events=10)
    evGen = ptcl_gun(s)

    out = tmp_path / "particles_edm4hep.root"

    out.mkdir()

    s.addWriter(
        conf_const(
            EDM4hepParticleWriter,
            acts.logging.INFO,
            inputParticles=evGen.config.outputParticles,
            outputPath=str(out),
        )
    )

    s.run()

    assert os.path.isfile(out)
    assert os.stat(out).st_size > 200


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_multitrajectory_writer(tmp_path):
    from acts.examples.edm4hep import EDM4hepMultiTrajectoryWriter

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
                / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
            )
        ),
        outputDir=tmp_path,
        s=s,
    )

    s.addAlgorithm(
        acts.examples.TracksToTrajectories(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTrajectories="trajectories",
        )
    )

    out = tmp_path / "trajectories_edm4hep.root"

    s.addWriter(
        EDM4hepMultiTrajectoryWriter(
            level=acts.logging.VERBOSE,
            inputTrajectories="trajectories",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputPath=str(out),
        )
    )

    s.run()

    assert os.path.isfile(out)
    assert os.stat(out).st_size > 200


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_tracks_writer(tmp_path):
    from acts.examples.edm4hep import EDM4hepTrackWriter

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
                / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
            )
        ),
        outputDir=tmp_path,
        s=s,
    )

    out = tmp_path / "tracks_edm4hep.root"

    s.addWriter(
        EDM4hepTrackWriter(
            level=acts.logging.VERBOSE,
            inputTracks="kf_tracks",
            outputPath=str(out),
            Bz=2 * u.T,
        )
    )

    s.run()

    assert os.path.isfile(out)
    assert os.stat(out).st_size > 200

    if not podioEnabled:
        import warnings

        warnings.warn(
            "edm4hep output checks were skipped, because podio was not on the python path"
        )
        return

    from podio.root_io import Reader
    import cppyy

    reader = Reader(str(out))

    actual = []

    for frame in reader.get("events"):
        tracks = frame.get("ActsTracks")
        for track in tracks:
            actual.append(
                (track.getChi2(), track.getNdf(), len(track.getTrackStates()))
            )

            locs = []

            perigee = None
            for ts in track.getTrackStates():
                if ts.location == cppyy.gbl.edm4hep.TrackState.AtIP:
                    perigee = ts
                    continue
                locs.append(ts.location)

                rp = ts.referencePoint
                r = math.sqrt(rp.x**2 + rp.y**2)
                assert r > 25

            assert locs[0] == cppyy.gbl.edm4hep.TrackState.AtLastHit
            assert locs[-1] == cppyy.gbl.edm4hep.TrackState.AtFirstHit

            assert perigee is not None
            rp = perigee.referencePoint
            assert rp.x == 0.0
            assert rp.y == 0.0
            assert rp.z == 0.0
            assert abs(perigee.D0) < 1e0
            assert abs(perigee.Z0) < 1e1


def generate_input_test_edm4hep_simhit_reader(input, output):
    from DDSim.DD4hepSimulation import DD4hepSimulation

    ddsim = DD4hepSimulation()
    if isinstance(ddsim.compactFile, list):
        ddsim.compactFile = [input]
    else:
        ddsim.compactFile = input
    ddsim.enableGun = True
    ddsim.gun.direction = (1, 0, 0)
    ddsim.gun.particle = "pi-"
    ddsim.gun.distribution = "eta"
    ddsim.numberOfEvents = 10
    ddsim.outputFile = output
    ddsim.outputConfig.forceEDM4HEP = True
    ddsim.run()


@pytest.mark.slow
@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_simhit_particle_reader(tmp_path):
    from acts.examples.edm4hep import EDM4hepSimReader

    tmp_file = str(tmp_path / "output_edm4hep.root")
    odd_xml_file = str(getOpenDataDetectorDirectory() / "xml" / "OpenDataDetector.xml")

    # explicitly ask for "spawn" as CI failures were observed with "fork"
    spawn_context = multiprocessing.get_context("spawn")
    p = spawn_context.Process(
        target=generate_input_test_edm4hep_simhit_reader, args=(odd_xml_file, tmp_file)
    )
    p.start()
    p.join()

    assert os.path.exists(tmp_file)

    s = Sequencer(numThreads=1)

    with getOpenDataDetector() as detector:
        trackingGeometry = detector.trackingGeometry()

        s.addReader(
            EDM4hepSimReader(
                level=acts.logging.INFO,
                inputPath=tmp_file,
                inputSimHits=[
                    "PixelBarrelReadout",
                    "PixelEndcapReadout",
                    "ShortStripBarrelReadout",
                    "ShortStripEndcapReadout",
                    "LongStripBarrelReadout",
                    "LongStripEndcapReadout",
                ],
                outputParticlesGenerator="particles_generated",
                outputParticlesSimulation="particles_simulated",
                outputSimHits="simhits",
                dd4hepDetector=detector,
                trackingGeometry=trackingGeometry,
            )
        )

        alg = AssertCollectionExistsAlg("simhits", "check_alg", acts.logging.WARNING)
        s.addAlgorithm(alg)

        alg = AssertCollectionExistsAlg(
            "particles_generated", "check_alg", acts.logging.WARNING
        )
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
            inputPath=str(out),
        )
    )

    algs = [
        AssertCollectionExistsAlg(k, f"check_alg_{k}", acts.logging.WARNING)
        for k in ("measurements", "simhitsmap")
    ]
    for alg in algs:
        s.addAlgorithm(alg)

    s.run()

    for alg in algs:
        assert alg.events_seen == 10


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_tracks_reader(tmp_path):
    from acts.examples.edm4hep import EDM4hepTrackWriter, EDM4hepTrackReader

    detector = acts.examples.GenericDetector()
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
                / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
            )
        ),
        outputDir=tmp_path,
        s=s,
    )

    out = tmp_path / "tracks_edm4hep.root"

    s.addWriter(
        EDM4hepTrackWriter(
            level=acts.logging.VERBOSE,
            inputTracks="kf_tracks",
            outputPath=str(out),
            Bz=2 * u.T,
        )
    )

    s.run()

    s = Sequencer(numThreads=1)
    s.addReader(
        EDM4hepTrackReader(
            level=acts.logging.VERBOSE,
            outputTracks="kf_tracks",
            inputPath=str(out),
            Bz=2 * u.T,
        )
    )

    s.run()
