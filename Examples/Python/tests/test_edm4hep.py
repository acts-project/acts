import os
import multiprocessing
from pathlib import Path
import math
import tempfile
import shutil
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


def assert_podio(
    target_file: Path,
    category: str,
    collections: set[str] | None = None,
    nevents: int | None = None,
):
    assert target_file.exists(), f"File {target_file} does not exist"
    from podio.root_io import Reader
    import cppyy

    reader = Reader(str(target_file))
    assert (
        category in reader.categories
    ), f"Category {category} not found in {target_file} ({reader.categories})"

    if nevents is not None:
        assert (
            len(reader.get(category)) == nevents
        ), f"Expected {nevents} events in {target_file} ({category}) but got {len(reader.get(category))}"

    if collections is not None:
        for frame in reader.get(category):
            assert (
                set(frame.getAvailableCollections()) == collections
            ), f"Expected collections {collections} in {target_file} ({category}) but got {frame.getAvailableCollections()}"
            break


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_measurement_writer(tmp_path, fatras):
    from acts.examples.edm4hep import EDM4hepMeasurementOutputConverter
    from acts.examples.podio import PodioWriter

    s = Sequencer(numThreads=1, events=10)
    _, simAlg, digiAlg = fatras(s)

    out = tmp_path / "measurements_edm4hep.root"

    converter = EDM4hepMeasurementOutputConverter(
        level=acts.logging.VERBOSE,
        inputMeasurements=digiAlg.config.outputMeasurements,
        inputClusters=digiAlg.config.outputClusters,
        outputTrackerHitsPlane="tracker_hits_plane",
        outputTrackerHitsRaw="tracker_hits_raw",
    )
    s.addAlgorithm(converter)

    s.addWriter(
        PodioWriter(
            level=acts.logging.VERBOSE,
            outputPath=str(out),
            category="events",
            collections=converter.collections,
        )
    )

    s.run()

    assert os.path.isfile(out)
    assert os.stat(out).st_size > 10

    assert_podio(
        out,
        "events",
        collections=set(["tracker_hits_plane", "tracker_hits_raw"]),
        nevents=10,
    )


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_simhit_writer(tmp_path, fatras, conf_const):
    from acts.examples.edm4hep import EDM4hepSimHitOutputConverter
    from acts.examples.podio import PodioWriter

    s = Sequencer(numThreads=1, events=10)
    _, simAlg, _ = fatras(s)

    out = tmp_path / "simhits_edm4hep.root"

    converter = EDM4hepSimHitOutputConverter(
        level=acts.logging.INFO,
        inputSimHits=simAlg.config.outputSimHits,
        outputSimTrackerHits="sim_tracker_hits",
    )
    s.addAlgorithm(converter)

    s.addWriter(
        PodioWriter(
            level=acts.logging.INFO,
            outputPath=str(out),
            category="events",
            collections=converter.collections,
        )
    )
    s.run()

    assert os.path.isfile(out)
    assert os.stat(out).st_size > 200


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_particle_writer(tmp_path, ptcl_gun):
    from acts.examples.edm4hep import EDM4hepParticleOutputConverter
    from acts.examples.podio import PodioWriter

    s = Sequencer(numThreads=1, events=10)
    _, h3conv = ptcl_gun(s)

    out = tmp_path / "particles_edm4hep.root"

    out.mkdir()

    converter = EDM4hepParticleOutputConverter(
        acts.logging.INFO,
        inputParticles=h3conv.config.outputParticles,
        outputParticles="MCParticles",
    )
    s.addAlgorithm(converter)

    s.addWriter(
        PodioWriter(
            level=acts.logging.INFO,
            outputPath=str(out),
            category="events",
            collections=converter.collections,
        )
    )

    s.run()

    assert os.path.isfile(out)
    assert os.stat(out).st_size > 200

    assert_podio(out, "events", collections=set(["MCParticles"]), nevents=10)


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_multitrajectory_writer(tmp_path):
    from acts.examples.edm4hep import EDM4hepMultiTrajectoryOutputConverter
    from acts.examples.podio import PodioWriter

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

    s.addAlgorithm(
        acts.examples.TracksToTrajectories(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTrajectories="trajectories",
        )
    )

    out = tmp_path / "trajectories_edm4hep.root"

    converter = EDM4hepMultiTrajectoryOutputConverter(
        level=acts.logging.VERBOSE,
        inputTrajectories="trajectories",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTracks="ActsTrajectories",
    )
    s.addAlgorithm(converter)

    s.addWriter(
        PodioWriter(
            level=acts.logging.VERBOSE,
            outputPath=str(out),
            category="events",
            collections=converter.collections,
        )
    )
    s.run()

    assert os.path.isfile(out)
    assert os.stat(out).st_size > 200


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_tracks_writer(tmp_path):
    from acts.examples.edm4hep import EDM4hepTrackOutputConverter
    from acts.examples.podio import PodioWriter

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

    out = tmp_path / "tracks_edm4hep.root"

    converter = EDM4hepTrackOutputConverter(
        level=acts.logging.VERBOSE,
        inputTracks="kf_tracks",
        outputTracks="ActsTracks",
        Bz=2 * u.T,
    )
    s.addAlgorithm(converter)

    s.addWriter(
        PodioWriter(
            level=acts.logging.VERBOSE,
            outputPath=str(out),
            category="events",
            collections=converter.collections,
        )
    )
    s.run()

    assert os.path.isfile(out), f"File {out} does not exist"
    assert os.stat(out).st_size > 200, f"File {out} is too small"

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

    ddsim.random.seed = 37

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


# Session scoped fixture that uses a temp folder
@pytest.fixture(scope="session")
def ddsim_input_session():
    with tempfile.TemporaryDirectory() as tmp_dir:
        odd_xml_file = str(
            getOpenDataDetectorDirectory() / "xml" / "OpenDataDetector.xml"
        )

        output_file = str(Path(tmp_dir) / "output_edm4hep.root")

        if not os.path.exists(output_file):
            # explicitly ask for "spawn" as CI failures were observed with "fork"
            spawn_context = multiprocessing.get_context("spawn")
            p = spawn_context.Process(
                target=generate_input_test_edm4hep_simhit_reader,
                args=(odd_xml_file, output_file),
            )
            p.start()
            p.join()
            if p.exitcode != 0:
                raise RuntimeError("ddsim process failed")

            assert os.path.exists(output_file)

        yield output_file


# Function scoped fixture that uses a temp folder
@pytest.fixture(scope="function")
def ddsim_input(ddsim_input_session, tmp_path):
    tmp_file = str(tmp_path / "output_edm4hep.root")
    shutil.copy(ddsim_input_session, tmp_file)
    return tmp_file


@pytest.mark.slow
@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_simhit_particle_reader(tmp_path, ddsim_input):
    from acts.examples.edm4hep import EDM4hepSimInputConverter
    from acts.examples.podio import PodioReader

    s = Sequencer(numThreads=1)

    with getOpenDataDetector() as detector:
        trackingGeometry = detector.trackingGeometry()

        s.addReader(
            PodioReader(
                level=acts.logging.VERBOSE,
                inputPath=ddsim_input,
                outputFrame="events",
                category="events",
            )
        )

        s.addAlgorithm(
            EDM4hepSimInputConverter(
                level=acts.logging.INFO,
                inputFrame="events",
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
                outputSimVertices="simvertices",
                dd4hepDetector=detector,
                trackingGeometry=trackingGeometry,
            )
        )

        alg = AssertCollectionExistsAlg(
            ["simhits", "simvertices", "particles_generated", "particles_simulated"],
            "check_alg",
            acts.logging.WARNING,
        )
        s.addAlgorithm(alg)

        alg = AssertCollectionExistsAlg(
            "particles_generated", "check_alg", acts.logging.WARNING
        )
        s.addAlgorithm(alg)

        s.run()

    assert alg.events_seen == 10


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_measurement_reader(tmp_path, fatras):
    from acts.examples.edm4hep import (
        EDM4hepMeasurementOutputConverter,
        EDM4hepMeasurementInputConverter,
    )
    from acts.examples.podio import PodioWriter, PodioReader

    s = Sequencer(numThreads=1, events=10)
    _, simAlg, digiAlg = fatras(s)

    out = tmp_path / "measurements_edm4hep.root"

    converter = EDM4hepMeasurementOutputConverter(
        level=acts.logging.INFO,
        inputMeasurements=digiAlg.config.outputMeasurements,
        inputClusters=digiAlg.config.outputClusters,
    )
    s.addAlgorithm(converter)
    s.addWriter(
        PodioWriter(
            level=acts.logging.INFO,
            outputPath=str(out),
            category="events",
            collections=converter.collections,
        )
    )
    s.run()

    # read back in
    s = Sequencer(numThreads=1)

    s.addReader(
        PodioReader(
            level=acts.logging.WARNING,
            inputPath=str(out),
            outputFrame="events",
            category="events",
        )
    )
    s.addAlgorithm(
        EDM4hepMeasurementInputConverter(
            level=acts.logging.WARNING,
            inputFrame="events",
            outputMeasurements="measurements",
            outputMeasurementSimHitsMap="simhitsmap",
        )
    )

    alg = AssertCollectionExistsAlg(
        ["measurements", "simhitsmap"], "check_alg", acts.logging.WARNING
    )
    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 10


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_tracks_reader(tmp_path):
    from acts.examples.edm4hep import (
        EDM4hepTrackOutputConverter,
        EDM4hepTrackInputConverter,
    )
    from acts.examples.podio import PodioWriter, PodioReader

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
                / "Examples/Configs/generic-digi-smearing-config.json"
            )
        ),
        outputDir=tmp_path,
        s=s,
    )

    out = tmp_path / "tracks_edm4hep.root"

    converter = EDM4hepTrackOutputConverter(
        level=acts.logging.VERBOSE,
        inputTracks="kf_tracks",
        Bz=2 * u.T,
    )
    s.addAlgorithm(converter)

    s.addWriter(
        PodioWriter(
            level=acts.logging.VERBOSE,
            outputPath=str(out),
            category="events",
            collections=converter.collections,
        )
    )
    s.run()

    s = Sequencer(numThreads=1)

    s.addReader(
        PodioReader(
            level=acts.logging.VERBOSE,
            inputPath=str(out),
            outputFrame="events",
            category="events",
        )
    )

    s.addAlgorithm(
        EDM4hepTrackInputConverter(
            level=acts.logging.VERBOSE,
            inputFrame="events",
            inputTracks="kf_tracks",
            outputTracks="kf_tracks",
            Bz=2 * u.T,
        )
    )

    s.run()


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_reader(ddsim_input):
    from acts.examples.podio import PodioReader

    s = Sequencer(numThreads=1)
    s.addReader(
        PodioReader(
            level=acts.logging.VERBOSE,
            inputPath=ddsim_input,
            outputFrame="frame",
        )
    )

    alg = AssertCollectionExistsAlg("frame", "check_alg", acts.logging.WARNING)
    s.addAlgorithm(alg)

    s.run()


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_writer_copy(ddsim_input, tmp_path):
    from acts.examples.podio import PodioWriter, PodioReader

    target_file = tmp_path / "rewritten_edm4hep.root"

    s = Sequencer(numThreads=1)
    s.addReader(
        PodioReader(
            level=acts.logging.VERBOSE,
            inputPath=ddsim_input,
            outputFrame="myframe",
            category="events",
        )
    )
    s.addWriter(
        PodioWriter(
            level=acts.logging.VERBOSE,
            inputFrame="myframe",
            outputPath=str(target_file),
            category="myevents",
        )
    )

    s.run()

    assert_podio(
        target_file,
        "myevents",
        nevents=10,
        collections=set(
            [
                "ShortStripEndcapReadout",
                "MCParticles",
                "LongStripBarrelReadout",
                "PixelEndcapReadout",
                "HCalEndcapCollectionContributions",
                "HCalEndcapCollection",
                "HCalBarrelCollectionContributions",
                "ECalBarrelCollectionContributions",
                "LongStripEndcapReadout",
                "EventHeader",
                "ECalEndcapCollectionContributions",
                "ShortStripBarrelReadout",
                "PixelBarrelReadout",
                "ECalEndcapCollection",
                "HCalBarrelCollection",
                "ECalBarrelCollection",
            ]
        ),
    )
