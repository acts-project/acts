import math
import warnings
from pathlib import Path

import pytest

import acts
import acts.examples
from acts import UnitConstants as u
from acts.examples import Sequencer

from helpers import AssertCollectionExistsAlg, arrowEnabled, isCI, pythia8Enabled

pytestmark = pytest.mark.skipif(
    not arrowEnabled, reason="Arrow/Parquet bindings not built"
)


def test_coexist_with_pyarrow():
    """Verify acts.examples.arrow and pyarrow can be loaded into the same
    Python process. Regression guard for the linker-isolation design: if
    libActsPluginArrow leaked any arrow symbols, or if pyarrow's bundled
    libarrow got shadowed by a differently-built libarrow, this import
    sequence would fail with a missing-vtable / missing-symbol ImportError.

    Skipped if pyarrow isn't installed."""
    pytest.importorskip("pyarrow")

    import pyarrow as pa
    from acts.examples.arrow import ParquetWriter  # forces the arrow island

    # Exercise pyarrow end-to-end to prove its libarrow loaded correctly
    # (not just pyarrow's __init__).
    table = pa.table({"x": [1, 2, 3]})
    assert table.num_rows == 3


PARTICLE_FIELDS = {
    "particle_id",
    "pdg_id",
    "mass",
    "energy",
    "charge",
    "vx",
    "vy",
    "vz",
    "time",
    "px",
    "py",
    "pz",
    "perigee_d0",
    "perigee_z0",
    "vertex_primary",
    "parent_id",
    "primary",
}


def with_pyarrow(fn):
    """
    Some assertions use pyarrow to inspect the written Parquet files.
    Locally, if pyarrow is not available (or incompatible with the linked
    libarrow), we skip these checks with a warning. In CI mode, missing
    pyarrow becomes a failure.
    """
    try:
        import pyarrow
        import pyarrow.parquet as pq

        fn(pyarrow, pq)
    except ImportError as e:
        if isCI:
            raise
        warnings.warn(f"pyarrow not available ({e}), skipping parquet checks")


def _assert_particles_parquet(directory: Path, expected_events: int) -> None:
    """Verify a particles dataset directory has the expected nested schema
    and row count summed across shard files."""
    assert directory.exists(), f"{directory} does not exist"
    assert directory.is_dir(), f"{directory} is not a directory"

    fragments = sorted(directory.glob("*.parquet"))
    assert fragments, f"{directory} contains no parquet shards"
    for f in fragments:
        assert f.stat().st_size > 0, f"{f} is empty"

    @with_pyarrow
    def _check(pa, pq):
        total_rows = 0
        union_event_ids: list[int] = []
        first_schema = None
        for f in fragments:
            pf = pq.ParquetFile(str(f))
            total_rows += pf.metadata.num_rows
            if first_schema is None:
                first_schema = pf.schema_arrow
            else:
                assert pf.schema_arrow.equals(
                    first_schema
                ), f"{f.name}: schema differs from {fragments[0].name}"
            t = pf.read(columns=["event_id"])
            union_event_ids.extend(t.column("event_id").to_pylist())

        assert total_rows == expected_events, (
            f"{directory.name}: expected {expected_events} events, " f"got {total_rows}"
        )

        # Each event id should appear exactly once across the dataset
        # (validates the writer's shard routing).
        assert sorted(union_event_ids) == list(range(expected_events)), (
            f"{directory.name}: event ids not contiguous 0..{expected_events-1}: "
            f"{sorted(union_event_ids)}"
        )

        names = {first_schema.field(i).name for i in range(len(first_schema))}
        assert "event_id" in names, f"{directory.name}: event_id column missing"
        missing = PARTICLE_FIELDS - names
        assert not missing, f"{directory.name}: missing fields {missing}"

        # Every particle column should be a list<T> in the nested layout.
        for field_name in PARTICLE_FIELDS:
            ftype = first_schema.field(field_name).type
            assert pa.types.is_list(
                ftype
            ), f"{directory.name}: field '{field_name}' should be list, got {ftype}"

        # At least one event should contain particles.
        table = pq.ParquetDataset(str(directory)).read()
        counts = [len(row) for row in table.column("particle_id").to_pylist()]
        assert any(
            c > 0 for c in counts
        ), f"{directory.name}: all events are empty ({counts})"


def _assert_parquet_reader_config(
    inputDir: Path,
    collections: dict[str, str],
    expectedSchemas: dict[str, "acts.arrow.ArrowSchema"],
    expected_events: int,
) -> None:
    from acts.examples.arrow import ParquetReader

    reader = ParquetReader(
        level=acts.logging.INFO,
        inputDir=str(inputDir),
        collections=collections,
        expectedSchemas=expectedSchemas,
    )

    assert reader.availableEvents() == (0, expected_events)


def _add_arrow_writer(
    s: Sequencer,
    outputDir: Path,
    inputs_to_tables: dict[str, str],
    eventsPerShard: int = 2,
) -> None:
    """Wire one ArrowParticleOutputConverter per (input, table) pair, and one
    ParquetWriter picking up all the resulting tables.

    @param inputs_to_tables: maps whiteboard key of SimParticleContainer to
        the desired whiteboard key / filename basename of the Arrow table.
        These must differ — the same key can't hold both an ACTS container
        and an arrow::Table.
    """
    from acts.arrow import particleSchema
    from acts.examples.arrow import ArrowParticleOutputConverter, ParquetWriter

    for input_key, table_key in inputs_to_tables.items():
        assert (
            input_key != table_key
        ), "Arrow output key must differ from the SimParticleContainer key"
        s.addAlgorithm(
            ArrowParticleOutputConverter(
                level=acts.logging.INFO,
                inputParticles=input_key,
                outputTable=table_key,
            )
        )

    s.addWriter(
        ParquetWriter(
            level=acts.logging.INFO,
            outputDir=str(outputDir),
            collections={
                table_key: table_key for table_key in inputs_to_tables.values()
            },
            expectedSchemas={
                table_key: particleSchema() for table_key in inputs_to_tables.values()
            },
            eventsPerShard=eventsPerShard,
        )
    )


def test_particle_gun_generated(tmp_path, ptcl_gun):
    """Particle gun → generated particles → Parquet."""
    from acts.arrow import particleSchema

    nevents = 5
    s = Sequencer(numThreads=1, events=nevents)
    ptcl_gun(s)
    _add_arrow_writer(s, tmp_path, {"particles_generated": "particles_generated_arrow"})
    s.run()

    _assert_particles_parquet(tmp_path / "particles_generated_arrow", nevents)
    _assert_parquet_reader_config(
        tmp_path,
        {"particles_generated_arrow": "particles_generated_arrow"},
        {"particles_generated_arrow": particleSchema()},
        nevents,
    )


def test_particle_gun_roundtrip(tmp_path, ptcl_gun):
    """Write sharded Parquet, then drive a second Sequencer off ParquetReader
    and check the reader exposes — and processes — the same number of events
    that were written."""
    from acts.arrow import particleSchema
    from acts.examples.arrow import ParquetReader

    # nevents/eventsPerShard chosen so the write phase produces multiple
    # shards with a non-full final shard — exercises shard discovery, the
    # multi-fragment scan, and the partial-shard edge case in one test.
    nevents = 5
    events_per_shard = 2
    expected_shards = (nevents + events_per_shard - 1) // events_per_shard

    s_write = Sequencer(numThreads=1, events=nevents)
    ptcl_gun(s_write)
    _add_arrow_writer(
        s_write,
        tmp_path,
        {"particles_generated": "particles_generated_arrow"},
        eventsPerShard=events_per_shard,
    )
    s_write.run()

    out_dir = tmp_path / "particles_generated_arrow"
    _assert_particles_parquet(out_dir, nevents)

    shards = sorted(out_dir.glob("*.parquet"))
    assert len(shards) == expected_shards, (
        f"expected {expected_shards} shards for {nevents} events at "
        f"{events_per_shard} events/shard, got {len(shards)}: "
        f"{[s.name for s in shards]}"
    )

    reader = ParquetReader(
        level=acts.logging.INFO,
        inputDir=str(tmp_path),
        collections={"particles_generated_arrow": "particles_generated_arrow"},
        expectedSchemas={"particles_generated_arrow": particleSchema()},
    )
    assert reader.availableEvents() == (0, nevents)

    # No `events=` — the sequencer derives the event range from the reader,
    # so a wrong count here would surface as a mismatch with `nevents`.
    s_read = Sequencer(numThreads=1)
    s_read.addReader(reader)
    counter = AssertCollectionExistsAlg(
        collections="particles_generated_arrow",
        name="roundtrip_check",
        level=acts.logging.INFO,
    )
    s_read.addAlgorithm(counter)
    s_read.run()

    assert counter.events_seen == nevents, (
        f"reader-driven sequencer processed {counter.events_seen} events, "
        f"expected {nevents}"
    )


def test_particle_gun_roundtrip_parallel_cache(tmp_path, ptcl_gun):
    """Write a multi-shard Parquet dataset, then read it back with a
    multi-threaded Sequencer and shardCacheCapacity > 1. Exercises the
    ParquetDatasetReader's concurrent shard-cache path (double-checked
    locking + LRU eviction across multiple resident shards) end-to-end,
    which the single-threaded/default-capacity tests above don't reach."""
    from acts.arrow import particleSchema
    from acts.examples.arrow import ParquetReader

    nevents = 12
    events_per_shard = 2  # -> 6 shards, more than numThreads/cacheCapacity
    num_threads = 3
    cache_capacity = 3

    s_write = Sequencer(numThreads=1, events=nevents)
    ptcl_gun(s_write)
    _add_arrow_writer(
        s_write,
        tmp_path,
        {"particles_generated": "particles_generated_arrow"},
        eventsPerShard=events_per_shard,
    )
    s_write.run()

    out_dir = tmp_path / "particles_generated_arrow"
    _assert_particles_parquet(out_dir, nevents)

    reader = ParquetReader(
        level=acts.logging.INFO,
        inputDir=str(tmp_path),
        collections={"particles_generated_arrow": "particles_generated_arrow"},
        expectedSchemas={"particles_generated_arrow": particleSchema()},
        shardCacheCapacity=cache_capacity,
    )
    assert reader.availableEvents() == (0, nevents)

    s_read = Sequencer(numThreads=num_threads)
    s_read.addReader(reader)
    counter = AssertCollectionExistsAlg(
        collections="particles_generated_arrow",
        name="roundtrip_check_parallel",
        level=acts.logging.INFO,
    )
    s_read.addAlgorithm(counter)
    s_read.run()

    assert counter.events_seen == nevents, (
        f"reader-driven parallel sequencer processed {counter.events_seen} events, "
        f"expected {nevents}"
    )


def test_particle_gun_fatras(tmp_path, fatras):
    """Particle gun + Fatras → both generated and simulated particles → Parquet."""
    nevents = 5
    s = Sequencer(numThreads=1, events=nevents)
    fatras(s)
    _add_arrow_writer(
        s,
        tmp_path,
        {
            "particles_generated": "particles_generated_arrow",
            "particles_simulated": "particles_simulated_arrow",
        },
    )
    s.run()

    _assert_particles_parquet(tmp_path / "particles_generated_arrow", nevents)
    _assert_particles_parquet(tmp_path / "particles_simulated_arrow", nevents)


SIMHIT_FIELDS = {
    "x",
    "y",
    "z",
    "true_x",
    "true_y",
    "true_z",
    "time",
    "particle_id",
    "detector",
    "volume_id",
    "layer_id",
    "surface_id",
}


def _add_simhit_arrow_writer(
    s: Sequencer,
    outputDir: Path,
    table_key: str = "simhits_arrow",
    *,
    withClusters: bool,
    eventsPerShard: int = 2,
) -> None:
    """Wire an ArrowSimHitOutputConverter + ParquetWriter for the simhits.

    @param withClusters: if True, feed the cluster container and the
        sim-hit→measurement map so the digitized x,y,z columns are filled from
        the precomputed Cluster::globalPosition; otherwise leave them unwired so
        those columns come out NaN.
    """
    from acts.arrow import simHitSchema
    from acts.examples.arrow import ArrowSimHitOutputConverter, ParquetWriter

    kwargs = {}
    if withClusters:
        kwargs["inputClusters"] = "clusters"
        kwargs["inputSimHitMeasurementsMap"] = "simhit_measurements_map"

    s.addAlgorithm(
        ArrowSimHitOutputConverter(
            level=acts.logging.INFO,
            inputSimHits="simhits",
            inputParticles="particles_simulated",
            outputTable=table_key,
            **kwargs,
        )
    )

    s.addWriter(
        ParquetWriter(
            level=acts.logging.INFO,
            outputDir=str(outputDir),
            collections={table_key: table_key},
            expectedSchemas={table_key: simHitSchema()},
            eventsPerShard=eventsPerShard,
        )
    )


def _assert_simhits_parquet(
    directory: Path, expected_events: int, *, expect_digitized: bool
) -> None:
    """Verify a simhits dataset directory: schema, row count, and the
    digitized-position semantics.

    The truth columns (true_x/true_y/true_z/time) are always finite. The
    digitized columns (x/y/z) are finite only for hits matched to a cluster
    when clusters were wired in; otherwise every value is NaN.
    """
    assert directory.exists(), f"{directory} does not exist"
    assert directory.is_dir(), f"{directory} is not a directory"

    fragments = sorted(directory.glob("*.parquet"))
    assert fragments, f"{directory} contains no parquet shards"
    for f in fragments:
        assert f.stat().st_size > 0, f"{f} is empty"

    @with_pyarrow
    def _check(pa, pq):
        first_schema = pq.ParquetFile(str(fragments[0])).schema_arrow
        names = {first_schema.field(i).name for i in range(len(first_schema))}
        assert "event_id" in names, f"{directory.name}: event_id column missing"
        missing = SIMHIT_FIELDS - names
        assert not missing, f"{directory.name}: missing fields {missing}"
        for field_name in SIMHIT_FIELDS:
            ftype = first_schema.field(field_name).type
            assert pa.types.is_list(
                ftype
            ), f"{directory.name}: field '{field_name}' should be list, got {ftype}"

        table = pq.ParquetDataset(str(directory)).read()
        assert table.num_rows == expected_events, (
            f"{directory.name}: expected {expected_events} events, "
            f"got {table.num_rows}"
        )

        x = table.column("x").to_pylist()
        y = table.column("y").to_pylist()
        z = table.column("z").to_pylist()
        tx = table.column("true_x").to_pylist()
        ty = table.column("true_y").to_pylist()
        tz = table.column("true_z").to_pylist()

        total_hits = sum(len(row) for row in tx)
        assert total_hits > 0, f"{directory.name}: no sim hits across any event"

        # Truth position is always written, so it must be finite everywhere.
        for col in (tx, ty, tz):
            for row in col:
                for v in row:
                    assert math.isfinite(v), f"{directory.name}: non-finite truth {v}"

        n_digitized = 0
        for ex, ey, ez, etx, ety, etz in zip(x, y, z, tx, ty, tz):
            for xv, yv, zv, txv, tyv, tzv in zip(ex, ey, ez, etx, ety, etz):
                if not expect_digitized:
                    assert (
                        math.isnan(xv) and math.isnan(yv) and math.isnan(zv)
                    ), f"{directory.name}: expected NaN digitized pos, got ({xv},{yv},{zv})"
                    continue
                # With clusters wired, matched hits carry a finite digitized
                # position; unmatched ones stay NaN. A finite x must come with a
                # finite y and z.
                if math.isnan(xv):
                    continue
                n_digitized += 1
                assert math.isfinite(yv) and math.isfinite(
                    zv
                ), f"{directory.name}: partial digitized pos ({xv},{yv},{zv})"
                # The position must sit inside the generic-detector envelope
                # (in mm) — a units regression (e.g. m vs mm) or a stale/garbage
                # position would blow well past this.
                r = math.hypot(xv, yv)
                assert r < 1500.0 and abs(zv) < 4000.0, (
                    f"{directory.name}: digitized pos ({xv},{yv},{zv}) outside "
                    f"detector envelope — units or stale-position regression?"
                )
                # It must also stay near the truth hit's module. The measured
                # coordinate(s) pin it to the truth, but the unmeasured strip
                # coordinate defaults to the surface centre, so allow a generous
                # module-scale slack rather than a tight smearing tolerance.
                dist = math.sqrt((xv - txv) ** 2 + (yv - tyv) ** 2 + (zv - tzv) ** 2)
                assert dist < 250.0, (
                    f"{directory.name}: digitized pos {dist:.1f} mm from truth "
                    f"({txv},{tyv},{tzv}) — wrong surface or unit error?"
                )

        if expect_digitized:
            assert n_digitized > 0, (
                f"{directory.name}: clusters were wired but no hit got a "
                f"digitized position"
            )


def test_fatras_simhits_digitized(tmp_path, fatras):
    """Fatras + digitization → ArrowSimHitOutputConverter reads cluster
    positions → Parquet. The matched-hit x,y,z must be the precomputed cluster
    global positions (finite, near truth)."""
    nevents = 3
    s = Sequencer(numThreads=1, events=nevents)
    fatras(s)
    _add_simhit_arrow_writer(s, tmp_path, withClusters=True)
    s.run()

    _assert_simhits_parquet(tmp_path / "simhits_arrow", nevents, expect_digitized=True)


def test_fatras_simhits_no_clusters_are_nan(tmp_path, fatras):
    """Without the cluster container and sim-hit→measurement map wired, the
    digitized x,y,z columns fall back to NaN while truth positions still
    populate."""
    nevents = 3
    s = Sequencer(numThreads=1, events=nevents)
    fatras(s)
    _add_simhit_arrow_writer(s, tmp_path, withClusters=False)
    s.run()

    _assert_simhits_parquet(tmp_path / "simhits_arrow", nevents, expect_digitized=False)


@pytest.mark.skipif(not pythia8Enabled, reason="Pythia8 not built")
def test_pythia8_fatras(tmp_path, rng, trk_geo):
    """Pythia8 ttbar + Fatras → generated AND simulated particles → Parquet."""
    from acts.examples.simulation import (
        addPythia8,
        addGenParticleSelection,
        ParticleSelectorConfig,
    )

    nevents = 3
    s = Sequencer(numThreads=1, events=nevents)

    vtxGen = acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(50 * u.um, 50 * u.um, 150 * u.mm, 20 * u.ns),
        mean=acts.Vector4(0, 0, 0, 0),
    )

    addPythia8(
        s,
        rnd=rng,
        hardProcess=["Top:qqbar2ttbar=on"],
        npileup=2,
        vtxGen=vtxGen,
        outputDirCsv=None,
        outputDirRoot=None,
        logLevel=acts.logging.WARNING,
    )

    # Pythia8 produces many particles per event (including those outside the
    # detector acceptance or below threshold). Prune them before handing to
    # Fatras — this matches the pattern in full_chain_odd.py.
    addGenParticleSelection(
        s,
        ParticleSelectorConfig(
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
            eta=(-3.0, 3.0),
            pt=(150 * u.MeV, None),
        ),
    )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    s.addAlgorithm(
        acts.examples.FatrasSimulation(
            level=acts.logging.WARNING,
            inputParticles="particles_generated_selected",
            outputParticles="particles_simulated",
            outputSimHits="simhits",
            randomNumbers=rng,
            trackingGeometry=trk_geo,
            magneticField=field,
            generateHitsOnSensitive=True,
            emScattering=False,
            emEnergyLossIonisation=False,
            emEnergyLossRadiation=False,
            emPhotonConversion=False,
        )
    )

    _add_arrow_writer(
        s,
        tmp_path,
        {
            "particles_generated": "particles_generated_arrow",
            "particles_simulated": "particles_simulated_arrow",
        },
    )
    s.run()

    _assert_particles_parquet(tmp_path / "particles_generated_arrow", nevents)
    _assert_particles_parquet(tmp_path / "particles_simulated_arrow", nevents)


def test_reader_schema_evolution_added_optional_column(tmp_path):
    """Read shards written without an optional column and verify the reader
    materializes it as null.

    Concretely: hand-write track shards using the production track schema
    *minus* the optional `t` field, then drive a sequencer with
    `ParquetReader` configured with the *full* track schema (which has `t`
    as a nullable column). The dataset scanner should project missing
    columns to null per fragment, so the table on the whiteboard must
    carry `t` and it must be all-null.

    This is the canonical added-optional-column schema-evolution case.
    """
    pa = pytest.importorskip("pyarrow")
    pq = pytest.importorskip("pyarrow.parquet")

    from acts.arrow import ArrowTable, trackSchema
    from acts.examples import ReadDataHandle
    from acts.examples.arrow import ParquetReader

    # Take the production track schema as the consumer's view.
    full_track_schema_pa = pa.schema(trackSchema())
    assert "t" in full_track_schema_pa.names, (
        f"trackSchema() unexpectedly lacks the 't' field; this test relies "
        f"on it being present. Schema:\n{full_track_schema_pa}"
    )

    # "Old" producer schema = production schema MINUS `t`, derived from the
    # full schema rather than rebuilt by hand so the two stay in sync if
    # the underlying definition evolves.
    old_track_schema = full_track_schema_pa.remove(
        full_track_schema_pa.get_field_index("t")
    )

    nevents = 4
    events_per_shard = 2
    collection_dir = tmp_path / "tracks_arrow"
    collection_dir.mkdir()

    # Build one row per event, with a single track per event for simplicity.
    # The on-disk schema additionally needs `event_id` prepended (the reader
    # uses it for filter pushdown).
    event_id_field = pa.field("event_id", pa.uint32(), nullable=False)
    on_disk_schema = pa.schema([event_id_field, *list(old_track_schema)])

    def field_type(name: str) -> "pa.DataType":
        return old_track_schema.field(name).type

    def make_event_table(event_id: int) -> "pa.Table":
        return pa.table(
            {
                "event_id": pa.array([event_id], type=pa.uint32()),
                "d0": pa.array([[0.1]], type=field_type("d0")),
                "z0": pa.array([[0.2]], type=field_type("z0")),
                "phi": pa.array([[0.3]], type=field_type("phi")),
                "theta": pa.array([[0.4]], type=field_type("theta")),
                "qop": pa.array([[0.5]], type=field_type("qop")),
                "majority_particle_id": pa.array(
                    [[1]], type=field_type("majority_particle_id")
                ),
                "hit_ids": pa.array([[[1, 2, 3]]], type=field_type("hit_ids")),
                "track_id": pa.array([[7]], type=field_type("track_id")),
            },
            schema=on_disk_schema,
        )

    # Write two shards, [0,2) and [2,4), each with one row group per event.
    for shard_start in range(0, nevents, events_per_shard):
        shard_end = shard_start + events_per_shard
        shard_path = (
            collection_dir / f"tracks_{shard_start:06d}-{shard_end:06d}.parquet"
        )
        with pq.ParquetWriter(str(shard_path), on_disk_schema) as writer:
            for event_id in range(shard_start, shard_end):
                writer.write_table(make_event_table(event_id))

    # Sanity check: the on-disk shards genuinely lack `t`.
    for shard_path in sorted(collection_dir.glob("*.parquet")):
        on_disk = pq.ParquetFile(str(shard_path)).schema_arrow
        assert "t" not in on_disk.names, (
            f"{shard_path.name}: precondition broken, on-disk shard "
            f"unexpectedly contains 't'. Schema: {on_disk}"
        )

    reader = ParquetReader(
        level=acts.logging.INFO,
        inputDir=str(tmp_path),
        collections={"tracks_arrow": "tracks_arrow"},
        expectedSchemas={"tracks_arrow": trackSchema()},
    )
    assert reader.availableEvents() == (0, nevents)

    # Pure-Python checker: pulls the ArrowTable off the WhiteBoard via the
    # typed registry, crosses into pyarrow zero-copy via __arrow_c_array__,
    # and inspects the result with pyarrow's full API.
    class TrackTableCheck(acts.examples.IAlgorithm):
        events_seen = 0

        def __init__(self, name="TrackTableCheck"):
            super().__init__(name=name, level=acts.logging.INFO)
            self._handle = ReadDataHandle(self, ArrowTable, "tracks_arrow")
            self._handle.initialize("tracks_arrow")

        def execute(self, ctx):
            handle = self._handle(ctx.eventStore)
            t = handle.as_table()
            assert "t" in t.column_names, (
                f"event {ctx.eventNumber}: 't' column missing from "
                f"projected table; schema: {t.schema}"
            )
            t_col = t.column("t")
            assert t_col.null_count == t_col.length(), (
                f"event {ctx.eventNumber}: 't' expected all-null, got "
                f"{t_col.length() - t_col.null_count} non-null of "
                f"{t_col.length()} values"
            )
            for required in ("d0", "z0", "phi", "theta", "qop"):
                assert required in t.column_names, (
                    f"event {ctx.eventNumber}: required column " f"'{required}' missing"
                )
            type(self).events_seen += 1
            return acts.examples.ProcessCode.SUCCESS

    s = Sequencer(numThreads=1)
    s.addReader(reader)
    s.addAlgorithm(TrackTableCheck())
    s.run()

    assert (
        TrackTableCheck.events_seen == nevents
    ), f"checker saw {TrackTableCheck.events_seen} events, expected {nevents}"


def test_python_alg_writes_arrow_table(tmp_path):
    """Smoke test for the write direction.

    A pure-Python algorithm constructs a per-event pyarrow table, wraps it
    via `ArrowTable.from_arrow`, and writes it onto the WhiteBoard through
    a typed `WriteDataHandle`. A second pure-Python algorithm reads it back
    via `ReadDataHandle`, slurps it into pyarrow via `as_table()`, and
    asserts the values survived the round-trip.

    Exercises: C-Data import (pyarrow → ArrowTable), `WhiteBoardRegistry`
    fromPython (ArrowTable → WhiteBoard storage), C-Data export (ArrowTable
    → pyarrow). End-to-end zero-copy across two libarrow instances.
    """
    pa = pytest.importorskip("pyarrow")

    from acts.arrow import ArrowTable, trackSchema
    from acts.examples import ReadDataHandle, WriteDataHandle

    # Use the production track schema as the inter-algorithm contract.
    track_schema_pa = pa.schema(trackSchema())

    def field_type(name):
        return track_schema_pa.field(name).type

    class TrackProducer(acts.examples.IAlgorithm):
        """Builds one row per event in the production track schema and
        writes it onto the whiteboard as an ArrowTable."""

        def __init__(self, key, name="TrackProducer"):
            super().__init__(name=name, level=acts.logging.INFO)
            self._out = WriteDataHandle(self, ArrowTable, key)
            self._out.initialize(key)

        def execute(self, ctx):
            evt = float(ctx.eventNumber)
            pa_table = pa.table(
                {
                    "d0": pa.array([[0.1 + evt]], type=field_type("d0")),
                    "z0": pa.array([[0.2 + evt]], type=field_type("z0")),
                    "phi": pa.array([[0.3]], type=field_type("phi")),
                    "theta": pa.array([[0.4]], type=field_type("theta")),
                    "qop": pa.array([[0.5]], type=field_type("qop")),
                    "majority_particle_id": pa.array(
                        [[1]], type=field_type("majority_particle_id")
                    ),
                    "hit_ids": pa.array([[[1, 2, 3]]], type=field_type("hit_ids")),
                    "track_id": pa.array([[7]], type=field_type("track_id")),
                    "t": pa.array([None], type=field_type("t")),
                },
                schema=track_schema_pa,
            )
            self._out(ctx, ArrowTable.from_arrow(pa_table))
            return acts.examples.ProcessCode.SUCCESS

    class TrackConsumer(acts.examples.IAlgorithm):
        """Reads the table back, exports through C-Data into pyarrow, and
        asserts the per-event values match what TrackProducer wrote."""

        events_seen = 0

        def __init__(self, key, name="TrackConsumer"):
            super().__init__(name=name, level=acts.logging.INFO)
            self._in = ReadDataHandle(self, ArrowTable, key)
            self._in.initialize(key)

        def execute(self, ctx):
            evt = float(ctx.eventNumber)
            t = self._in(ctx.eventStore).as_table()
            assert t.num_rows == 1
            d0 = t.column("d0").to_pylist()[0]
            z0 = t.column("z0").to_pylist()[0]
            assert d0 == [
                pytest.approx(0.1 + evt)
            ], f"event {ctx.eventNumber}: d0 round-trip mismatch: {d0}"
            assert z0 == [
                pytest.approx(0.2 + evt)
            ], f"event {ctx.eventNumber}: z0 round-trip mismatch: {z0}"
            type(self).events_seen += 1
            return acts.examples.ProcessCode.SUCCESS

    nevents = 3
    s = Sequencer(numThreads=1, events=nevents)
    s.addAlgorithm(TrackProducer(key="produced_tracks_arrow"))
    s.addAlgorithm(TrackConsumer(key="produced_tracks_arrow"))
    s.run()

    assert (
        TrackConsumer.events_seen == nevents
    ), f"consumer saw {TrackConsumer.events_seen} events, expected {nevents}"


def test_writer_rejects_missing_schema(tmp_path):
    """ParquetWriter requires an expected schema for every collection.
    Constructing one without one must fail at config time, not at run time.
    """
    from acts.examples.arrow import ParquetWriter

    with pytest.raises(ValueError, match="no expected schema"):
        ParquetWriter(
            level=acts.logging.INFO,
            outputDir=str(tmp_path),
            collections={"some_collection": "some_collection"},
            expectedSchemas={},  # missing entry for "some_collection"
        )


def test_writer_aborts_on_per_event_schema_mismatch(tmp_path):
    """A pure-Python algorithm produces a table whose schema doesn't match
    the writer's declared expectedSchemas. The writer must abort the
    sequencer with a clear message rather than silently writing garbage.
    """
    pa = pytest.importorskip("pyarrow")

    from acts.arrow import ArrowTable, particleSchema
    from acts.examples import WriteDataHandle
    from acts.examples.arrow import ParquetWriter

    # Producer writes a 1-row table with a single int column. Whatever the
    # schema, it isn't particleSchema(), which is what the writer is told
    # to expect below — so the per-event check must fire on event 0.
    class WrongShapeProducer(acts.examples.IAlgorithm):
        def __init__(self, key, name="WrongShapeProducer"):
            super().__init__(name=name, level=acts.logging.INFO)
            self._out = WriteDataHandle(self, ArrowTable, key)
            self._out.initialize(key)

        def execute(self, ctx):
            wrong = pa.table({"unexpected": pa.array([1], type=pa.int32())})
            self._out(ctx, ArrowTable.from_arrow(wrong))
            return acts.examples.ProcessCode.SUCCESS

    s = Sequencer(numThreads=1, events=1)
    s.addAlgorithm(WrongShapeProducer(key="bogus_arrow"))
    s.addWriter(
        ParquetWriter(
            level=acts.logging.INFO,
            outputDir=str(tmp_path),
            collections={"bogus_arrow": "bogus_arrow"},
            expectedSchemas={"bogus_arrow": particleSchema()},
        )
    )

    # The writer ABORTs, which the Sequencer turns into a runtime error.
    with pytest.raises(RuntimeError):
        s.run()
