from pathlib import Path

import pytest

import acts
import acts.examples
from acts import UnitConstants as u

from helpers import arrowEnabled, dd4hepEnabled

pytestmark = [
    pytest.mark.skipif(not arrowEnabled, reason="Arrow/Parquet bindings not built"),
    pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up"),
]

# Deterministic per-event perigee parameters written by
# _write_colliderml_tracks_dataset() below, checked for round-trip fidelity
# by test_colliderml_track_reading().
_TRACK_D0_BASE = 0.1
_TRACK_Z0_BASE = 0.2
_TRACK_PHI = 0.3
_TRACK_THETA = 1.2
_TRACK_QOP = 0.5


def _read_event_hit_counts(hits_dir: Path, n_events: int) -> dict:
    """Number of tracker hits per event in a ColliderML-format hits Parquet
    dataset, read back through the same ParquetReader/ArrowTable machinery
    ColliderMLRelease1InputConverter itself uses.

    Deliberately not read with `pyarrow.parquet` directly: `ParquetWriter`
    always writes ZSTD-compressed row groups (`ArrowUtil.cpp`), and this
    environment's pyarrow wheel lacks the zstd codec, so a direct
    `pyarrow.parquet` read of the on-disk file fails even though the
    dataset is perfectly valid. Going through `ParquetReader` decompresses
    with ActsPluginArrow's own (zstd-capable) Arrow build; only the
    resulting in-memory table is then handed to pyarrow via the C-Data
    interface, which involves no codec at all.
    """
    from acts.arrow import ArrowTable
    from acts.examples.arrow import ColliderMLRelease1InputConverter, ParquetReader

    counts = {}

    class HitCounter(acts.examples.IAlgorithm):
        def __init__(self):
            super().__init__("HitCounter", acts.logging.WARNING)
            self.hits = acts.examples.ReadDataHandle(self, ArrowTable, "InputHits")
            self.hits.initialize("cml_hits")

        def execute(self, context):
            table = self.hits(context.eventStore).as_table()
            counts[context.eventNumber] = len(table.column("x")[0].as_py())
            return acts.examples.ProcessCode.SUCCESS

    s = acts.examples.Sequencer(
        events=n_events, numThreads=1, logLevel=acts.logging.WARNING
    )
    s.addReader(
        ParquetReader(
            level=acts.logging.WARNING,
            collections={"cml_hits": str(hits_dir)},
            expectedSchemas={"cml_hits": ColliderMLRelease1InputConverter.hitSchema()},
        )
    )
    s.addAlgorithm(HitCounter())
    s.run()
    return counts


def _write_colliderml_tracks_dataset(
    tracks_dir: Path, n_hits_by_event: dict, n_events: int
):
    """Hand-write a ColliderML-Release-1-format tracks Parquet dataset: one
    track per event, on a schema matching
    ColliderMLRelease1InputConverter.tracksSchema() (own, ColliderML-owned
    contract — deliberately not ActsPlugins::ArrowUtil::trackSchema(), which
    is a different, ACTS-internal round-trip format).

    Each track's `hit_ids` references the event's first tracker hit (row 0
    in the corresponding hits table) when that event has any hits (per
    @p n_hits_by_event, from `_read_event_hit_counts`), so the reader's
    hit_ids -> measurement-index resolution is exercised against real,
    already-digitizable ColliderML hits; otherwise it's left empty.
    """
    pa = pytest.importorskip("pyarrow")
    pq = pytest.importorskip("pyarrow.parquet")

    tracks_dir.mkdir(parents=True, exist_ok=True)

    schema = pa.schema(
        [
            pa.field("event_id", pa.uint32(), nullable=False),
            pa.field("d0", pa.list_(pa.float32()), nullable=False),
            pa.field("z0", pa.list_(pa.float32()), nullable=False),
            pa.field("phi", pa.list_(pa.float32()), nullable=False),
            pa.field("theta", pa.list_(pa.float32()), nullable=False),
            pa.field("qop", pa.list_(pa.float32()), nullable=False),
            pa.field("majority_particle_id", pa.list_(pa.uint64()), nullable=False),
            pa.field("hit_ids", pa.list_(pa.list_(pa.uint32())), nullable=False),
            pa.field("track_id", pa.list_(pa.uint16()), nullable=False),
        ]
    )

    def make_event_table(event_id: int) -> "pa.Table":
        hit_ids = [0] if n_hits_by_event.get(event_id, 0) > 0 else []
        return pa.table(
            {
                "event_id": pa.array([event_id], type=pa.uint32()),
                "d0": pa.array(
                    [[_TRACK_D0_BASE + event_id]], type=schema.field("d0").type
                ),
                "z0": pa.array(
                    [[_TRACK_Z0_BASE + event_id]], type=schema.field("z0").type
                ),
                "phi": pa.array([[_TRACK_PHI]], type=schema.field("phi").type),
                "theta": pa.array([[_TRACK_THETA]], type=schema.field("theta").type),
                "qop": pa.array([[_TRACK_QOP]], type=schema.field("qop").type),
                "majority_particle_id": pa.array(
                    [[0]], type=schema.field("majority_particle_id").type
                ),
                "hit_ids": pa.array([[hit_ids]], type=schema.field("hit_ids").type),
                "track_id": pa.array([[0]], type=schema.field("track_id").type),
            },
            schema=schema,
        )

    shard_path = tracks_dir / "tracks_000000.parquet"
    with pq.ParquetWriter(str(shard_path), schema) as writer:
        for event_id in range(n_events):
            writer.write_table(make_event_table(event_id))

    return n_hits_by_event


_srcdir = Path(__file__).resolve().parent.parent.parent.parent

# Number of events generated by the FATRAS fixture and processed by the test.
_N_EVENTS = 3


@pytest.fixture(scope="session")
def colliderml_fatras_sample(tmp_path_factory):
    """Generate a small ColliderML-compatible parquet dataset using FATRAS + ODD.

    Runs a minimal sequencer: ParametricParticleGenerator → FATRAS → Digitization
    → ArrowParticleOutputConverter + ArrowSimHitOutputConverter → ParquetWriter.

    Returns (base_dir, sample_name).
    """
    from acts.arrow import particleSchema, simHitSchema
    from acts.examples.arrow import (
        ArrowParticleOutputConverter,
        ArrowSimHitOutputConverter,
        ParquetWriter,
    )
    from acts.examples.json import readDigiConfigFromJson
    from acts.examples.odd import getOpenDataDetector

    sample = "fatras_test"
    tmp = tmp_path_factory.mktemp("colliderml_fatras")

    particles_dir = tmp / f"{sample}_particles" / "data" / f"{sample}_particles"
    hits_dir = tmp / f"{sample}_tracker_hits" / "data" / f"{sample}_tracker_hits"
    particles_dir.mkdir(parents=True)
    hits_dir.mkdir(parents=True)

    rng = acts.examples.RandomNumbers(seed=42)
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    with getOpenDataDetector() as detector:
        tgeo = detector.trackingGeometry()

        s = acts.examples.Sequencer(
            events=_N_EVENTS,
            numThreads=1,
            logLevel=acts.logging.WARNING,
            outputDir=str(tmp),
            failOnUnmaskedFpe=False,
        )

        # --- Particle gun: 10 muons per event, pT=2 GeV, full eta/phi ---
        evGen = acts.examples.EventGenerator(
            level=acts.logging.WARNING,
            generators=[
                acts.examples.EventGenerator.Generator(
                    multiplicity=acts.examples.FixedMultiplicityGenerator(n=1),
                    vertex=acts.examples.GaussianVertexGenerator(
                        stddev=acts.Vector4(0, 0, 0, 0),
                        mean=acts.Vector4(0, 0, 0, 0),
                    ),
                    particles=acts.examples.ParametricParticleGenerator(
                        p=(2 * u.GeV, 2 * u.GeV),
                        eta=(-2, 2),
                        phi=(0, 360 * u.degree),
                        pdg=acts.PdgParticle.eMuon,
                        randomizeCharge=True,
                        numParticles=10,
                    ),
                )
            ],
            outputEvent="gun_event",
            randomNumbers=rng,
        )
        s.addReader(evGen)

        hepmc3Conv = acts.examples.hepmc3.HepMC3InputConverter(
            level=acts.logging.WARNING,
            inputEvent=evGen.config.outputEvent,
            outputParticles="particles_gen",
            outputVertices="vertices_gen",
        )
        s.addAlgorithm(hepmc3Conv)

        # --- FATRAS simulation ---
        fatrasAlg = acts.examples.FatrasSimulation(
            level=acts.logging.WARNING,
            inputParticles="particles_gen",
            outputParticles="particles_sim",
            outputSimHits="simhits",
            randomNumbers=rng,
            trackingGeometry=tgeo,
            magneticField=field,
            generateHitsOnSensitive=True,
            emScattering=False,
            emEnergyLossIonisation=False,
            emEnergyLossRadiation=False,
            emPhotonConversion=False,
        )
        s.addAlgorithm(fatrasAlg)

        # --- Digitization (needed for valid x/y/z cluster positions) ---
        digiCfg = acts.examples.DigitizationAlgorithm.Config(
            digitizationConfigs=readDigiConfigFromJson(
                str(_srcdir / "Examples/Configs/odd-digi-smearing-config-notime.json")
            ),
            surfaceByIdentifier=tgeo.geoIdSurfaceMap(),
            randomNumbers=rng,
            inputSimHits="simhits",
        )
        s.addAlgorithm(
            acts.examples.DigitizationAlgorithm(digiCfg, acts.logging.WARNING)
        )

        # --- Arrow output converters ---
        s.addAlgorithm(
            ArrowParticleOutputConverter(
                level=acts.logging.WARNING,
                inputParticles="particles_sim",
                outputTable="cml_particles_arrow",
            )
        )
        s.addAlgorithm(
            ArrowSimHitOutputConverter(
                level=acts.logging.WARNING,
                inputSimHits="simhits",
                inputParticles="particles_sim",
                inputClusters="clusters",
                inputSimHitMeasurementsMap="simhit_measurements_map",
                outputTable="cml_hits_arrow",
            )
        )

        # --- ParquetWriter: absolute paths so files land in the nested layout ---
        s.addWriter(
            ParquetWriter(
                level=acts.logging.WARNING,
                outputDir=str(tmp),
                collections={
                    "cml_particles_arrow": str(particles_dir),
                    "cml_hits_arrow": str(hits_dir),
                },
                expectedSchemas={
                    "cml_particles_arrow": particleSchema(),
                    "cml_hits_arrow": simHitSchema(),
                },
                eventsPerShard=_N_EVENTS,
            )
        )

        s.run()

    return tmp, sample


@pytest.mark.parametrize("reco_geo", ["gen1", "gen3"])
def test_colliderml_truth_tracking(tmp_path, colliderml_fatras_sample, reco_geo):
    """Read FATRAS-generated ColliderML-format data and run truth-seeded KF tracking.

    Parametrized on the reconstruction geometry:
      gen1 — same Gen1 ODD used for simulation; no geo-id map needed.
      gen3 — Gen3 ODD; geo-id map and digi config generated on the fly.

    Verifies that the ColliderML reader pipeline (ParquetReader +
    ColliderMLRelease1InputConverter + TruthEstimated seeding + KF) runs without
    error.
    """
    import sys

    sys.path.insert(0, str(_srcdir / "Examples/Scripts/Python"))
    from acts.examples.odd import getOpenDataDetector
    from colliderml_truth_tracking import runColliderMLTruthTracking
    from generate_geoid_map import generate_geoid_map

    inputDir, sample = colliderml_fatras_sample
    particlesDir = inputDir / f"{sample}_particles" / "data" / f"{sample}_particles"
    hitsDir = inputDir / f"{sample}_tracker_hits" / "data" / f"{sample}_tracker_hits"
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

    def _run():
        s = runColliderMLTruthTracking(
            trackingGeometry=trackingGeometry,
            field=field,
            outputDir=tmp_path,
            particlesDir=particlesDir,
            hitsDir=hitsDir,
            geoIdMapPath=geoid_map_path if reco_geo == "gen3" else None,
            decorators=decorators,
            events=_N_EVENTS,
            numThreads=1,
            sample=sample,
        )
        s.run()

    if ctx is not None:
        with ctx:
            _run()
    else:
        _run()


@pytest.fixture(scope="session")
def colliderml_tracks_dataset(colliderml_fatras_sample):
    """Hand-written ColliderML-format published-tracks dataset, derived from
    the FATRAS fixture's already-written hits (one track per event,
    referencing the event's first hit if any)."""
    inputDir, sample = colliderml_fatras_sample
    hits_dir = inputDir / f"{sample}_tracker_hits" / "data" / f"{sample}_tracker_hits"
    tracks_dir = inputDir / f"{sample}_tracks" / "data" / f"{sample}_tracks"

    n_hits_by_event = _read_event_hit_counts(hits_dir, _N_EVENTS)
    _write_colliderml_tracks_dataset(tracks_dir, n_hits_by_event, _N_EVENTS)

    return inputDir, sample, tracks_dir, n_hits_by_event


@pytest.mark.pypi
def test_colliderml_track_reading(colliderml_tracks_dataset):
    """Read the published-tracks table via
    ColliderMLRelease1InputConverter(inputTracksTable=..., outputTracks=...)
    and verify the resulting ConstTrackContainer: every track carries a
    (0,0,0)-centered Perigee reference surface, its perigee parameters
    round-trip the values written to the tracks dataset, and
    nMeasurements() is nonzero exactly for tracks whose hit_ids resolved to
    a real measurement.
    """
    from acts.examples.arrow import ColliderMLRelease1InputConverter, ParquetReader
    from acts.examples.odd import getOpenDataDetector

    inputDir, sample, tracks_dir, n_hits_by_event = colliderml_tracks_dataset
    particles_dir = inputDir / f"{sample}_particles" / "data" / f"{sample}_particles"
    hits_dir = inputDir / f"{sample}_tracker_hits" / "data" / f"{sample}_tracker_hits"

    odd_dir = acts.examples.odd.getOpenDataDetectorDirectory()
    matDeco = acts.IMaterialDecorator.fromFile(
        odd_dir / "data/odd-material-maps.root", level=acts.logging.WARNING
    )
    detector = getOpenDataDetector(matDeco)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    results = {}

    class TrackCheck(acts.examples.IAlgorithm):
        def __init__(self):
            super().__init__("TrackCheck", acts.logging.WARNING)
            self.tracks = acts.examples.ReadDataHandle(
                self, acts.examples.ConstTrackContainer, "InputTracks"
            )
            self.tracks.initialize("colliderml_tracks")

        def execute(self, context):
            tracks = self.tracks(context.eventStore)
            entries = []
            for track in tracks:
                center = track.referenceSurface.center(context.geoContext)
                entries.append(
                    (
                        list(center),
                        list(track.parameters),
                        track.nMeasurements,
                    )
                )
            results[context.eventNumber] = entries
            return acts.examples.ProcessCode.SUCCESS

    with detector:
        s = acts.examples.Sequencer(
            events=_N_EVENTS, numThreads=1, logLevel=acts.logging.WARNING
        )
        for d in decorators:
            s.addContextDecorator(d)

        s.addReader(
            ParquetReader(
                level=acts.logging.WARNING,
                collections={
                    "cml_particles": str(particles_dir),
                    "cml_hits": str(hits_dir),
                    "cml_tracks": str(tracks_dir),
                },
                expectedSchemas={
                    "cml_particles": ColliderMLRelease1InputConverter.particleSchema(),
                    "cml_hits": ColliderMLRelease1InputConverter.hitSchema(),
                    "cml_tracks": ColliderMLRelease1InputConverter.tracksSchema(),
                },
            )
        )

        s.addAlgorithm(
            ColliderMLRelease1InputConverter(
                level=acts.logging.WARNING,
                inputParticlesTable="cml_particles",
                inputHitsTable="cml_hits",
                inputTracksTable="cml_tracks",
                outputMeasurements="measurements",
                outputMeasurementSubset="measurement_subset",
                outputMeasSimHitsMap="measurement_simhits_map",
                outputMeasParticlesMap="measurement_particles_map",
                outputParticleMeasurementsMap="particle_measurements_map",
                outputClusters="clusters",
                outputTracks="colliderml_tracks",
                trackingGeometry=trackingGeometry,
            )
        )
        s.addAlgorithm(TrackCheck())
        s.run()

    assert len(results) == _N_EVENTS
    for event, entries in results.items():
        assert len(entries) == 1
        center, params, n_meas = entries[0]
        assert center == pytest.approx([0.0, 0.0, 0.0], abs=1e-6)
        assert params[0] == pytest.approx(_TRACK_D0_BASE + event, abs=1e-5)
        assert params[1] == pytest.approx(_TRACK_Z0_BASE + event, abs=1e-5)
        assert params[2] == pytest.approx(_TRACK_PHI, abs=1e-5)
        assert params[3] == pytest.approx(_TRACK_THETA, abs=1e-5)
        assert params[4] == pytest.approx(_TRACK_QOP, abs=1e-5)
        if n_hits_by_event.get(event, 0) > 0:
            assert n_meas > 0, f"event {event}: expected a resolved hit"
        else:
            assert n_meas == 0


def test_colliderml_reader_track_output_requires_measurements(colliderml_fatras_sample):
    """outputTracks without outputMeasurements must raise at construction —
    the tracks' hit_ids can only be resolved against measurements built from
    the same hits table."""
    from acts.examples.arrow import ColliderMLRelease1InputConverter

    inputDir, sample = colliderml_fatras_sample
    odd_dir = acts.examples.odd.getOpenDataDetectorDirectory()
    matDeco = acts.IMaterialDecorator.fromFile(
        odd_dir / "data/odd-material-maps.root", level=acts.logging.WARNING
    )
    from acts.examples.odd import getOpenDataDetector

    detector = getOpenDataDetector(matDeco)
    with detector:
        trackingGeometry = detector.trackingGeometry()
        with pytest.raises(ValueError):
            ColliderMLRelease1InputConverter(
                level=acts.logging.WARNING,
                inputParticlesTable="cml_particles",
                inputHitsTable="cml_hits",
                inputTracksTable="cml_tracks",
                outputParticles="particles",
                outputTracks="colliderml_tracks",
                trackingGeometry=trackingGeometry,
            )


def test_colliderml_reader_works_without_tracks(colliderml_fatras_sample):
    """The reader must keep working unchanged (particles/hits/measurements
    only) when inputTracksTable/outputTracks are left unset — tracks are an
    optional capability, not a requirement."""
    from acts.examples.arrow import ColliderMLRelease1InputConverter, ParquetReader
    from acts.examples.odd import getOpenDataDetector

    inputDir, sample = colliderml_fatras_sample
    particles_dir = inputDir / f"{sample}_particles" / "data" / f"{sample}_particles"
    hits_dir = inputDir / f"{sample}_tracker_hits" / "data" / f"{sample}_tracker_hits"

    odd_dir = acts.examples.odd.getOpenDataDetectorDirectory()
    matDeco = acts.IMaterialDecorator.fromFile(
        odd_dir / "data/odd-material-maps.root", level=acts.logging.WARNING
    )
    detector = getOpenDataDetector(matDeco)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    with detector:
        s = acts.examples.Sequencer(
            events=_N_EVENTS, numThreads=1, logLevel=acts.logging.WARNING
        )
        for d in decorators:
            s.addContextDecorator(d)

        s.addReader(
            ParquetReader(
                level=acts.logging.WARNING,
                collections={
                    "cml_particles": str(particles_dir),
                    "cml_hits": str(hits_dir),
                },
                expectedSchemas={
                    "cml_particles": ColliderMLRelease1InputConverter.particleSchema(),
                    "cml_hits": ColliderMLRelease1InputConverter.hitSchema(),
                },
            )
        )
        s.addAlgorithm(
            ColliderMLRelease1InputConverter(
                level=acts.logging.WARNING,
                inputParticlesTable="cml_particles",
                inputHitsTable="cml_hits",
                outputParticles="particles",
                outputSimHits="simhits",
                trackingGeometry=trackingGeometry,
            )
        )
        s.run()
