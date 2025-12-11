from pathlib import Path
import tempfile
import sys
import os
import time
import functools
import warnings
import numpy

import pytest

import acts
import acts.examples
from acts.examples import Sequencer
from acts import UnitConstants as u

ScopedFailureThreshold = acts.logging.ScopedFailureThreshold

from helpers import (
    dd4hepEnabled,
    hepmc3Enabled,
    geant4Enabled,
    AssertCollectionExistsAlg,
    isCI,
)


def with_pyhepmc(fn):
    """
    Some tests use pyhepmc to inspect the written output files.
    Locally, if pyhepmc is not present, we ignore these checks with a warning.
    In CI mode, they become a failure
    """

    try:
        import pyhepmc

        fn(pyhepmc)
    except ImportError:
        if isCI:
            raise
        else:
            warnings.warn("pyhepmc not available, skipping checks")
            pass


pytestmark = [
    pytest.mark.hepmc3,
    pytest.mark.skipif(not hepmc3Enabled, reason="HepMC3 plugin not available"),
]

cm = acts.examples.hepmc3.Compression
from acts.examples.hepmc3 import availableCompressionModes

_available_compression_modes = availableCompressionModes()


compression_modes = pytest.mark.parametrize(
    "compression",
    acts.examples.hepmc3.availableCompressionModes(),
    ids=[
        c.name if c != cm.none else "uncompressed"
        for c in acts.examples.hepmc3.availableCompressionModes()
    ],
)


def handle_path(out, compression):
    if compression == cm.none:
        actual_path = out
    else:
        assert not out.exists()
        actual_path = out.with_suffix(
            f".hepmc3{acts.examples.hepmc3.compressionExtension(compression)}"
        )
        assert actual_path.suffix == acts.examples.hepmc3.compressionExtension(
            compression
        )
    assert actual_path.exists()
    return actual_path


def check_sidecar(actual_path, expected_events):
    """Check that a sidecar metadata file exists and contains the expected event count"""
    sidecar_path = Path(str(actual_path) + ".json")
    assert sidecar_path.exists(), f"Sidecar file {sidecar_path} does not exist"

    import json

    with sidecar_path.open("r") as f:
        metadata = json.load(f)

    assert "num_events" in metadata, "Sidecar file missing 'num_events' field"
    assert (
        metadata["num_events"] == expected_events
    ), f"Expected {expected_events} events, but sidecar says {metadata['num_events']}"


all_formats = []
all_format_ids = []

for compression in acts.examples.hepmc3.availableCompressionModes():
    all_formats.append(("hepmc3", compression))
    all_format_ids.append(
        f"hepmc3-{compression.name if compression != cm.none else 'uncompressed'}"
    )

all_formats.append(("root", cm.none))
all_format_ids.append("root")

all_formats = pytest.mark.parametrize(
    "format,compression", all_formats, ids=all_format_ids
)

main_formats = []
main_format_ids = []

for compression in acts.examples.hepmc3.availableCompressionModes():
    main_formats.append(("hepmc3", compression))
    main_format_ids.append(
        f"hepmc3-{compression.name if compression != cm.none else 'uncompressed'}"
    )

main_formats.append(("root", cm.none))
main_format_ids.append("root")


main_formats = pytest.mark.parametrize(
    "format,compression", main_formats, ids=main_format_ids
)


@all_formats
def test_hemc3_writer(tmp_path, rng, compression, format):
    from acts.examples.hepmc3 import (
        HepMC3Writer,
    )

    s = Sequencer(numThreads=10, events=100)

    evGen = acts.examples.EventGenerator(
        level=acts.logging.DEBUG,
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=2),
                vertex=acts.examples.GaussianVertexGenerator(
                    stddev=acts.Vector4(50 * u.um, 50 * u.um, 150 * u.mm, 20 * u.ns),
                    mean=acts.Vector4(0, 0, 0, 0),
                ),
                particles=acts.examples.ParametricParticleGenerator(
                    p=(100 * u.GeV, 100 * u.GeV),
                    eta=(-2, 2),
                    phi=(0, 360 * u.degree),
                    randomizeCharge=True,
                    numParticles=2,
                ),
            )
        ],
        outputEvent="hepmc3_event",
        randomNumbers=rng,
    )
    s.addReader(evGen)

    s.addAlgorithm(
        acts.examples.hepmc3.HepMC3InputConverter(
            level=acts.logging.DEBUG,
            inputEvent=evGen.config.outputEvent,
            outputParticles="particles_generated",
            outputVertices="vertices_truth",
        )
    )

    alg = AssertCollectionExistsAlg(
        [
            "particles_generated",
            "vertices_truth",
        ],
        "check_alg",
        acts.logging.WARNING,
    )
    s.addAlgorithm(alg)

    out = tmp_path / "out" / f"pytest.{format}"
    out.parent.mkdir(parents=True, exist_ok=True)

    s.addWriter(
        HepMC3Writer(
            acts.logging.DEBUG,
            inputEvent="hepmc3_event",
            outputPath=out,
            compression=compression,
            # writeEventsInOrder=False,
        )
    )

    s.run()

    actual_path = handle_path(out, compression)

    # Check sidecar metadata file
    check_sidecar(actual_path, s.config.events)

    # pyhepmc does not support zstd
    if compression in (cm.none, cm.lzma, cm.bzip2, cm.zlib) and format != "root":

        @with_pyhepmc
        def check(pyhepmc):
            nevts = 0
            event_numbers = []
            with pyhepmc.open(actual_path) as f:
                for evt in f:
                    nevts += 1
                    event_numbers.append(evt.event_number)
                    assert len(evt.particles) == 4 + 2  # muons + beam particles
            assert nevts == s.config.events
            assert event_numbers == list(
                range(s.config.events)
            ), "Events are out of order"


class StallAlgorithm(acts.examples.IAlgorithm):
    sleep: float = 1

    def execute(self, ctx):

        if ctx.eventNumber == 50:
            print("BEGIN SLEEP")
            time.sleep(self.sleep)
            print("END OF SLEEP")
        else:
            time.sleep(0.01)

        return acts.examples.ProcessCode.SUCCESS


@pytest.fixture
def common_evgen(rng):
    def func(s: Sequencer):
        evGen = acts.examples.EventGenerator(
            level=acts.logging.INFO,
            generators=[
                acts.examples.EventGenerator.Generator(
                    multiplicity=acts.examples.FixedMultiplicityGenerator(n=2),
                    vertex=acts.examples.GaussianVertexGenerator(
                        stddev=acts.Vector4(
                            50 * u.um, 50 * u.um, 150 * u.mm, 20 * u.ns
                        ),
                        mean=acts.Vector4(0, 0, 0, 0),
                    ),
                    particles=acts.examples.ParametricParticleGenerator(
                        p=(100 * u.GeV, 100 * u.GeV),
                        eta=(-2, 2),
                        phi=(0, 360 * u.degree),
                        randomizeCharge=True,
                        numParticles=2,
                    ),
                )
            ],
            outputEvent="hepmc3_event",
            randomNumbers=rng,
        )
        s.addReader(evGen)
        return evGen

    return func


@pytest.fixture
def common_writer(tmp_path, common_evgen):
    def func(s: Sequencer, compression: acts.examples.hepmc3.Compression, format: str):
        from acts.examples.hepmc3 import HepMC3Writer

        evGen = common_evgen(s)

        out = tmp_path / "out" / f"events_pytest.{format}"
        out.parent.mkdir(parents=True, exist_ok=True)

        s.addWriter(
            HepMC3Writer(
                acts.logging.INFO,
                inputEvent=evGen.config.outputEvent,
                outputPath=out,
                compression=compression,
            )
        )

        return out

    return func


@pytest.mark.parametrize(
    "bufsize",
    [1, 5, 15, 50],
    ids=lambda v: f"buf{v}",
)
@pytest.mark.repeat(5)
@pytest.mark.timeout(5, method="thread")
def test_hepmc3_writer_stall(common_evgen, tmp_path, bufsize):
    """
    This test simulates that one event takes significantly longer than the
    other ones. In this case, the HepMC3 writer should keep a buffer of events
    and wait until the trailing event comes in.
    """
    from acts.examples.hepmc3 import (
        HepMC3Writer,
    )

    s = Sequencer(numThreads=10, events=150, logLevel=acts.logging.VERBOSE)

    evGen = common_evgen(s)

    out = tmp_path / "out" / "pytest.hepmc3"
    out.parent.mkdir(parents=True, exist_ok=True)

    stall = StallAlgorithm(name="stall_alg", level=acts.logging.INFO)
    s.addAlgorithm(stall)

    s.addWriter(
        HepMC3Writer(
            acts.logging.VERBOSE,
            inputEvent=evGen.config.outputEvent,
            outputPath=out,
            maxEventsPending=bufsize,
        )
    )

    s.run()

    # Check sidecar metadata file
    check_sidecar(out, s.config.events)

    @with_pyhepmc
    def check(pyhepmc):
        nevts = 0
        event_numbers = []
        with pyhepmc.open(out) as f:
            for evt in f:
                nevts += 1
                event_numbers.append(evt.event_number)
        assert nevts == s.config.events
        assert event_numbers == list(range(s.config.events)), "Events are out of order"


def test_hepmc3_writer_not_in_order(common_evgen, tmp_path):
    """
    Bypasses the event ordering. This test mainly checks that this code path completes
    """
    from acts.examples.hepmc3 import (
        HepMC3Writer,
    )

    s = Sequencer(numThreads=10, events=100)

    evGen = common_evgen(s)

    out = tmp_path / "out" / "pytest.hepmc3"
    out.parent.mkdir(parents=True, exist_ok=True)

    stall = StallAlgorithm(name="stall_alg", level=acts.logging.INFO)
    s.addAlgorithm(stall)

    s.addWriter(
        HepMC3Writer(
            acts.logging.VERBOSE,
            inputEvent=evGen.config.outputEvent,
            outputPath=out,
            maxEventsPending=5,
            writeEventsInOrder=False,
        )
    )

    s.run()

    # Check sidecar metadata file
    check_sidecar(out, s.config.events)

    @with_pyhepmc
    def check(pyhepmc):
        nevts = 0
        event_numbers = []
        with pyhepmc.open(out) as f:
            for evt in f:
                nevts += 1
                event_numbers.append(evt.event_number)
        assert nevts == s.config.events
        # We don't expect events to be in order, but they should be there
        assert set(event_numbers) == set(
            range(s.config.events)
        ), "Event numbers are different"


@main_formats
def test_hepmc3_writer_pythia8(tmp_path, rng, compression, format):
    from acts.examples.hepmc3 import (
        HepMC3Writer,
    )

    from acts.examples.simulation import addPythia8

    s = Sequencer(numThreads=1, events=3)

    vtxGen = acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(50 * u.um, 50 * u.um, 150 * u.mm, 20 * u.ns),
        mean=acts.Vector4(0, 0, 0, 0),
    )

    addPythia8(
        s,
        rnd=rng,
        hardProcess=["Top:qqbar2ttbar=on"],
        npileup=10,
        outputDirCsv=None,
        outputDirRoot=None,
        logLevel=acts.logging.INFO,
        vtxGen=vtxGen,
    )

    out = tmp_path / f"events.{format}"

    s.addWriter(
        HepMC3Writer(
            acts.logging.VERBOSE,
            inputEvent="pythia8-event",
            outputPath=out,
            compression=compression,
        )
    )

    # Assert particles and vertices are present
    alg = AssertCollectionExistsAlg(
        [
            "particles_generated",
            "vertices_truth",
            "pythia8-event",
        ],
        "check_alg",
        acts.logging.WARNING,
    )
    s.addAlgorithm(alg)

    s.run()

    actual_path = handle_path(out, compression)

    assert actual_path.exists(), f"File {actual_path} does not exist"

    # Check sidecar metadata file
    check_sidecar(actual_path, s.config.events)

    if compression in (cm.none, cm.lzma, cm.bzip2, cm.zlib) and format != "root":

        @with_pyhepmc
        def check(pyhepmc):
            nevts = 0
            with pyhepmc.open(actual_path) as f:
                for evt in f:
                    # we expect the event to be populated but the exact number is random
                    assert len(evt.particles) > 2000
                    nevts += 1
            assert nevts == s.config.events


@main_formats
def test_hepmc3_reader(common_writer, rng, compression, format):
    from acts.examples.hepmc3 import (
        HepMC3Reader,
    )

    nevents = 1200

    s = Sequencer(numThreads=10, events=nevents)

    out = common_writer(s, compression, format)

    s.run()

    actual_path = handle_path(out, compression)

    # Without external event number, we need to read all events
    # use multiple threads to test if seeking to the right event works
    s = Sequencer(numThreads=10)

    with acts.logging.ScopedFailureThreshold(acts.logging.ERROR):
        # We expect a warning about missing input event count
        s.addReader(
            HepMC3Reader(
                inputPath=actual_path,
                level=acts.logging.VERBOSE,
                outputEvent="hepmc3_event",
            )
        )

    s.addAlgorithm(
        acts.examples.hepmc3.HepMC3InputConverter(
            level=acts.logging.INFO,
            inputEvent="hepmc3_event",
            outputParticles="particles_read",
            outputVertices="vertices_read",
        )
    )

    alg = AssertCollectionExistsAlg(
        [
            "particles_read",
            "vertices_read",
            "hepmc3_event",
        ],
        "check_alg",
        acts.logging.WARNING,
    )
    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == nevents


@main_formats
def test_hepmc3_reader_explicit_num_events(common_writer, rng, compression, format):
    """
    The HepMC3Reader can be configured to read a specific number of events. If
    the explicit number of events is lower or equal to the actual number of
    events, reading works as expected.
    """
    from acts.examples.hepmc3 import (
        HepMC3Reader,
    )

    nevents = 1200
    s = Sequencer(numThreads=10, events=nevents)
    out = common_writer(s, compression, format)
    s.run()
    actual_path = handle_path(out, compression)

    # With external event number, we can read a specific event
    # Test 110 and 120 as both a lower number than is available and a higher number
    # than is available is valid
    for num in (nevents - 10, nevents):
        s = Sequencer(numThreads=10)

        s.addReader(
            HepMC3Reader(
                acts.logging.DEBUG,
                inputPath=actual_path,
                outputEvent="hepmc3_event",
                numEvents=num,
            )
        )

        alg = AssertCollectionExistsAlg(
            "hepmc3_event", "check_alg", acts.logging.WARNING
        )
        s.addAlgorithm(alg)

        s.run()

        assert alg.events_seen == num


@main_formats
def test_hepmc3_reader_explicit_num_events_too_large(
    common_writer, rng, compression, format
):
    """
    The HepMC3Reader can be configured to read a specific number of events. If
    the explicit number of events is larger than the actual number of events,
    the reader should throw an error.
    """
    from acts.examples.hepmc3 import (
        HepMC3Reader,
    )

    nevents = 1200
    s = Sequencer(numThreads=10, events=nevents)
    out = common_writer(s, compression, format)
    s.run()
    actual_path = handle_path(out, compression)

    s = Sequencer(numThreads=10)

    s.addReader(
        HepMC3Reader(
            acts.logging.DEBUG,
            inputPath=actual_path,
            outputEvent="hepmc3_event",
            numEvents=nevents + 10,
        )
    )

    with ScopedFailureThreshold(acts.logging.MAX):
        with pytest.raises(RuntimeError) as excinfo:
            s.run()

    assert "Failed to process event" in str(excinfo.value)


@main_formats
def test_hepmc3_reader_skip_events(common_writer, rng, compression, format):
    from acts.examples.hepmc3 import (
        HepMC3Reader,
    )

    nevents = 1200
    s = Sequencer(numThreads=10, events=nevents)
    out = common_writer(s, compression, format)
    s.run()
    actual_path = handle_path(out, compression)

    skip = 100
    s = Sequencer(numThreads=5, skip=skip, events=nevents - skip)

    with acts.logging.ScopedFailureThreshold(acts.logging.ERROR):
        # We expect a warning about missing input event count
        s.addReader(
            HepMC3Reader(
                inputPath=actual_path,
                level=acts.logging.DEBUG,
                outputEvent="hepmc3_event",
                maxEventBufferSize=10,
            )
        )

    class EventNumberCheckerAlg(acts.examples.IAlgorithm):
        events_seen = set()

        def __init__(
            self,
            name="check_alg",
            level=acts.logging.INFO,
            *args,
            **kwargs,
        ):
            acts.examples.IAlgorithm.__init__(
                self, name=name, level=level, *args, **kwargs
            )

        def execute(self, ctx):
            self.events_seen.add(ctx.eventNumber)
            return acts.examples.ProcessCode.SUCCESS

    alg = EventNumberCheckerAlg()
    s.addAlgorithm(alg)

    s.run()

    exp = set(range(skip, nevents))
    # We can't assert the correct number assignment without having to expose
    # more of the FW The HepMC3Reader internally asserts that the even number
    # read from the event is equal to the event number currently being
    # processed.
    assert alg.events_seen == exp


def test_hepmc3_compression_modes():
    assert cm.none in acts.examples.hepmc3.availableCompressionModes()


@compression_modes
def test_hepmc3_writer_compression_auto_detection(tmp_path, rng, compression):
    """Test that compression is automatically detected from the output path"""
    from acts.examples.hepmc3 import HepMC3Writer

    s = Sequencer(numThreads=1, events=10)

    evGen = acts.examples.EventGenerator(
        level=acts.logging.INFO,
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=1),
                vertex=acts.examples.GaussianVertexGenerator(
                    stddev=acts.Vector4(50 * u.um, 50 * u.um, 150 * u.mm, 20 * u.ns),
                    mean=acts.Vector4(0, 0, 0, 0),
                ),
                particles=acts.examples.ParametricParticleGenerator(
                    p=(100 * u.GeV, 100 * u.GeV),
                    eta=(-2, 2),
                    phi=(0, 360 * u.degree),
                    randomizeCharge=True,
                    numParticles=2,
                ),
            )
        ],
        outputEvent="hepmc3_event",
        randomNumbers=rng,
    )
    s.addReader(evGen)

    # Create output path WITH compression extension but WITHOUT specifying compression in config
    ext = acts.examples.hepmc3.compressionExtension(compression)
    out = tmp_path / f"events.hepmc3{ext}"
    out.parent.mkdir(parents=True, exist_ok=True)

    # Don't specify compression - it should be auto-detected from the path
    s.addWriter(
        HepMC3Writer(
            acts.logging.INFO,
            inputEvent="hepmc3_event",
            outputPath=out,
            # compression=None (default)
        )
    )

    s.run()

    # File should exist with the correct extension
    assert out.exists(), f"Output file {out} does not exist"

    # Check sidecar
    check_sidecar(out, s.config.events)

    # Verify we can read it back
    if compression in (cm.none, cm.lzma, cm.bzip2, cm.zlib):

        @with_pyhepmc
        def check(pyhepmc):
            nevts = 0
            with pyhepmc.open(out) as f:
                for evt in f:
                    nevts += 1
            assert nevts == s.config.events


def test_hepmc3_writer_compression_consistency(tmp_path):
    """Test that specifying both path extension and config compression must be consistent"""
    from acts.examples.hepmc3 import HepMC3Writer

    out = tmp_path / "out" / "events.hepmc3.gz"
    out.parent.mkdir(parents=True, exist_ok=True)

    # Path says .gz (zlib) but config says bzip2 - should fail
    with acts.logging.ScopedFailureThreshold(acts.logging.MAX), pytest.raises(
        ValueError
    ) as excinfo:
        HepMC3Writer(
            acts.logging.INFO,
            inputEvent="hepmc3_event",
            outputPath=out,
            compression=acts.examples.hepmc3.Compression.bzip2,
        )

    assert "Compression mismatch" in str(excinfo.value)
    assert "zlib" in str(excinfo.value)
    assert "bzip2" in str(excinfo.value)


def test_hepmc3_writer_compression_explicit_with_path(tmp_path, rng):
    """Test that you can specify compression in config and omit it from path"""
    from acts.examples.hepmc3 import HepMC3Writer

    s = Sequencer(numThreads=1, events=10)

    evGen = acts.examples.EventGenerator(
        level=acts.logging.INFO,
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=1),
                vertex=acts.examples.GaussianVertexGenerator(
                    stddev=acts.Vector4(50 * u.um, 50 * u.um, 150 * u.mm, 20 * u.ns),
                    mean=acts.Vector4(0, 0, 0, 0),
                ),
                particles=acts.examples.ParametricParticleGenerator(
                    p=(100 * u.GeV, 100 * u.GeV),
                    eta=(-2, 2),
                    phi=(0, 360 * u.degree),
                    randomizeCharge=True,
                    numParticles=2,
                ),
            )
        ],
        outputEvent="hepmc3_event",
        randomNumbers=rng,
    )
    s.addReader(evGen)

    # Specify path WITHOUT extension, but WITH compression in config
    out = tmp_path / "events.hepmc3"
    out.parent.mkdir(parents=True, exist_ok=True)

    s.addWriter(
        HepMC3Writer(
            acts.logging.INFO,
            inputEvent="hepmc3_event",
            outputPath=out,
            compression=acts.examples.hepmc3.Compression.zlib,
        )
    )

    s.run()

    # File should exist with compression extension added
    actual_path = tmp_path / "events.hepmc3.gz"
    assert actual_path.exists(), f"Output file {actual_path} does not exist"

    # Check sidecar
    check_sidecar(actual_path, s.config.events)


def test_hepmc3_writer_root_compression_error(tmp_path):
    """Test that an error is raised when trying to use compression with ROOT format"""
    from acts.examples.hepmc3 import HepMC3Writer

    out = tmp_path / "out" / "events.root"
    out.parent.mkdir(parents=True, exist_ok=True)

    with acts.logging.ScopedFailureThreshold(acts.logging.MAX), pytest.raises(
        ValueError
    ) as excinfo:
        HepMC3Writer(
            acts.logging.INFO,
            inputEvent="hepmc3_event",
            outputPath=out,
            compression=acts.examples.hepmc3.Compression.zlib,
        )

    assert "Compression not supported for ROOT format" in str(excinfo.value)


def test_hepmc3_reader_multiple_files(tmp_path, rng):
    from acts.examples.hepmc3 import HepMC3Writer, HepMC3Reader

    events = 100
    n_pileup = 10

    s = Sequencer(numThreads=10, events=events, logLevel=acts.logging.INFO)

    vtxGenZero = acts.examples.FixedVertexGenerator(
        fixed=acts.Vector4(0, 0, 0, 0),
    )

    vtxGen = acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(50 * u.um, 50 * u.um, 150 * u.mm, 20 * u.ns),
        mean=acts.Vector4(0, 0, 0, 0),
    )

    hard_scatter = acts.examples.EventGenerator(
        level=acts.logging.INFO,
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=1),
                vertex=vtxGenZero,
                particles=acts.examples.ParametricParticleGenerator(
                    p=(100 * u.GeV, 100 * u.GeV),
                    eta=(-2, 2),
                    phi=(0, 360 * u.degree),
                    randomizeCharge=True,
                    numParticles=2,
                ),
            )
        ],
        outputEvent="hard_scatter_event",
        randomNumbers=rng,
    )
    s.addReader(hard_scatter)

    out_hs = tmp_path / "out" / "events_pytest_hs.hepmc3"
    out_hs.parent.mkdir(parents=True, exist_ok=True)

    compression = acts.examples.hepmc3.Compression.bzip2

    s.addWriter(
        HepMC3Writer(
            acts.logging.INFO,
            inputEvent=hard_scatter.config.outputEvent,
            outputPath=out_hs,
            compression=compression,
        )
    )

    s.run()

    s = Sequencer(numThreads=10, events=events * n_pileup, logLevel=acts.logging.INFO)

    pileup = acts.examples.EventGenerator(
        level=acts.logging.INFO,
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=1),
                vertex=vtxGenZero,
                particles=acts.examples.ParametricParticleGenerator(
                    p=(100 * u.GeV, 100 * u.GeV),
                    eta=(-2, 2),
                    phi=(0, 360 * u.degree),
                    randomizeCharge=True,
                    numParticles=2,
                ),
            )
        ],
        outputEvent="pileup_event",
        randomNumbers=rng,
    )
    s.addReader(pileup)

    out_pu = tmp_path / "out" / "events_pytest_pu.hepmc3"
    out_pu.parent.mkdir(parents=True, exist_ok=True)

    s.addWriter(
        HepMC3Writer(
            acts.logging.INFO,
            inputEvent=pileup.config.outputEvent,
            outputPath=out_pu,
            compression=compression,
        )
    )

    s.run()

    act_hs = handle_path(out_hs, compression)
    act_pu = handle_path(out_pu, compression)

    # Check sidecar metadata files
    check_sidecar(act_hs, events)
    check_sidecar(act_pu, events * n_pileup)

    # do reading including merging and write combined file

    s = Sequencer(numThreads=10, logLevel=acts.logging.INFO)

    reader = HepMC3Reader(
        inputs=[
            HepMC3Reader.Input.Fixed(act_hs, 1),
            HepMC3Reader.Input.Fixed(act_pu, 10),
        ],
        level=acts.logging.VERBOSE,
        outputEvent="hepmc3_event",
        vertexGenerator=vtxGen,
        randomNumbers=rng,
    )
    s.addReader(reader)

    out_combined = tmp_path / "out" / "events_pytest_combined.hepmc3"
    out_combined.parent.mkdir(parents=True, exist_ok=True)

    s.addWriter(
        HepMC3Writer(
            acts.logging.INFO,
            inputEvent=reader.config.outputEvent,
            outputPath=out_combined,
            compression=compression,
        )
    )

    s.run()

    act_combined = handle_path(out_combined, compression)

    # Check sidecar metadata file for combined output
    check_sidecar(act_combined, events)

    @with_pyhepmc
    def check(pyhepmc):
        def get_vtx_pos(file: Path):
            with pyhepmc.open(file) as f:
                for evt in f:
                    for vtx in evt.vertices:
                        yield [vtx.position.x, vtx.position.y, vtx.position.z]

        hs = numpy.vstack(list(get_vtx_pos(act_hs))).T
        pu = numpy.vstack(list(get_vtx_pos(act_pu))).T
        combined = numpy.vstack(list(get_vtx_pos(act_combined))).T

        # NO smearing in hs and pu
        for arr in (hs, pu):
            vx, vy, vz = arr
            std = numpy.std(vx)
            assert std < 1 * u.um
            std = numpy.std(vy)
            assert std < 1 * u.um
            std = numpy.std(vz)
            assert std < 1 * u.um

        # Configured smearing in combined
        # Checked values are a bit lower than configured smearing due to limited stats
        vx, vy, vz = combined
        std = numpy.std(vx)
        assert std > 40 * u.um
        std = numpy.std(vy)
        assert std > 40 * u.um
        std = numpy.std(vz)
        assert std > 140 * u.mm
