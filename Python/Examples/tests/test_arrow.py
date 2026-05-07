import warnings
from pathlib import Path

import pytest

import acts
import acts.examples
from acts import UnitConstants as u
from acts.examples import Sequencer

from helpers import arrowEnabled, isCI, pythia8Enabled


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
            f"{directory.name}: expected {expected_events} events, "
            f"got {total_rows}"
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
    inputDir: Path, collections: dict[str, str], expected_events: int
) -> None:
    from acts.examples.arrow import ParquetReader

    reader = ParquetReader(
        level=acts.logging.INFO,
        inputDir=str(inputDir),
        collections=collections,
    )

    assert reader.availableEvents() == (0, expected_events)


def _add_arrow_writer(
    s: Sequencer,
    outputDir: Path,
    inputs_to_tables: dict[str, str],
) -> None:
    """Wire one ArrowParticleOutputConverter per (input, table) pair, and one
    ParquetWriter picking up all the resulting tables.

    @param inputs_to_tables: maps whiteboard key of SimParticleContainer to
        the desired whiteboard key / filename basename of the Arrow table.
        These must differ — the same key can't hold both an ACTS container
        and an arrow::Table.
    """
    from acts.examples.arrow import ArrowParticleOutputConverter, ParquetWriter

    field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
    for input_key, table_key in inputs_to_tables.items():
        assert (
            input_key != table_key
        ), "Arrow output key must differ from the SimParticleContainer key"
        s.addAlgorithm(
            ArrowParticleOutputConverter(
                level=acts.logging.INFO,
                inputParticles=input_key,
                outputTable=table_key,
                bField=field,
            )
        )

    s.addWriter(
        ParquetWriter(
            level=acts.logging.INFO,
            outputDir=str(outputDir),
            collections={
                table_key: table_key
                for table_key in inputs_to_tables.values()
            },
            eventsPerShard=2,
        )
    )


def test_particle_gun_generated(tmp_path, ptcl_gun):
    """Particle gun → generated particles → Parquet."""
    nevents = 5
    s = Sequencer(numThreads=1, events=nevents)
    ptcl_gun(s)
    _add_arrow_writer(s, tmp_path, {"particles_generated": "particles_generated_arrow"})
    s.run()

    _assert_particles_parquet(tmp_path / "particles_generated_arrow", nevents)
    _assert_parquet_reader_config(
        tmp_path,
        {"particles_generated_arrow": "particles_generated_arrow"},
        nevents,
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
