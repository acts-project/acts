from pathlib import Path
import shutil
import tempfile
import sys
import os

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
)


pytestmark = [
    pytest.mark.hepmc3,
    pytest.mark.skipif(not hepmc3Enabled, reason="HepMC3 plugin not available"),
]

try:
    cm = acts.examples.hepmc3.Compression
    from acts.examples.hepmc3 import availableCompressionModes

    _available_compression_modes = availableCompressionModes()
except ImportError:
    _available_compression_modes = []


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


@pytest.mark.parametrize("per_event", [True, False], ids=["per_event", "combined"])
@compression_modes
def test_hepmc3_writer(tmp_path, rng, per_event, compression):
    from acts.examples.hepmc3 import (
        HepMC3Writer,
    )

    s = Sequencer(numThreads=1, events=10)

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

    out = tmp_path / "out" / "pytest.hepmc3"
    out.parent.mkdir(parents=True, exist_ok=True)

    s.addWriter(
        HepMC3Writer(
            acts.logging.DEBUG,
            inputEvent="hepmc3_event",
            outputPath=out,
            perEvent=per_event,
            compression=compression,
        )
    )

    s.run()

    if per_event:
        files = list(out.parent.iterdir())
        assert len(files) == 10
        assert all(f.stem.startswith("event") for f in files)
        if compression == cm.none:
            assert all(f.suffix == ".hepmc3" for f in files)
            for f in files:
                with f.open("r") as f:
                    assert len(f.readlines()) == 25
        else:
            ext = acts.examples.hepmc3.compressionExtension(compression)
            assert all(f.name.endswith(f".hepmc3{ext}") for f in files)

    else:
        actual_path = handle_path(out, compression)

        if compression == cm.none:
            with actual_path.open("r") as f:
                assert len(f.readlines()) == 214

        # pyhepmc does not support zstd
        if compression in (cm.lzma, cm.bzip2, cm.zlib):
            try:
                import pyhepmc

                nevts = 0
                with pyhepmc.open(actual_path) as f:
                    nevts += 1
                    for evt in f:
                        assert len(evt.particles) == 4 + 2  # muons + beam particles
                assert nevts == 1

            except ImportError:
                pass


@compression_modes
def test_hepmc3_writer_pythia8(tmp_path, rng, compression):
    from acts.examples.hepmc3 import (
        HepMC3Writer,
    )

    from pythia8 import addPythia8

    s = Sequencer(numThreads=1, events=1)

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

    out = tmp_path / "events.hepmc3"

    s.addWriter(
        HepMC3Writer(
            acts.logging.VERBOSE,
            inputEvent="pythia8-event",
            outputPath=out,
            perEvent=False,
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
    if compression in (cm.lzma, cm.bzip2, cm.zlib):
        try:
            import pyhepmc

            nevts = 0
            with pyhepmc.open(actual_path) as f:
                nevts += 1
                for evt in f:
                    assert len(evt.particles) == 7433
            assert nevts == 1

        except ImportError:
            pass


@compression_modes
def test_hepmc3_reader(tmp_path, rng, compression):
    from acts.examples.hepmc3 import (
        HepMC3Writer,
        HepMC3Reader,
    )

    s = Sequencer(numThreads=1, events=12)

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
        randomNumbers=rng,
    )

    s.addReader(evGen)

    out = tmp_path / "out" / "events_pytest.hepmc3"
    out.parent.mkdir(parents=True, exist_ok=True)

    s.addWriter(
        HepMC3Writer(
            acts.logging.DEBUG,
            inputEvent="hepmc3_event",
            outputPath=out,
            perEvent=False,
            compression=compression,
        )
    )

    s.run()

    actual_path = handle_path(out, compression)

    # Without external event number, we need to read all events
    s = Sequencer(numThreads=1)

    s.addReader(
        HepMC3Reader(
            acts.logging.DEBUG,
            inputPath=actual_path,
            perEvent=False,
            outputEvent="hepmc3_event",
        )
    )

    s.addAlgorithm(
        acts.examples.hepmc3.HepMC3InputConverter(
            level=acts.logging.DEBUG,
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

    assert alg.events_seen == 12

    # With external event number, we can read a specific event
    # Test 11 and 12 as both a lower number than is available and a higher number
    # than is available is valid
    for num in (11, 12):
        s = Sequencer(numThreads=1)

        s.addReader(
            HepMC3Reader(
                acts.logging.DEBUG,
                inputPath=actual_path,
                perEvent=False,
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

    s = Sequencer(numThreads=1)

    s.addReader(
        HepMC3Reader(
            acts.logging.DEBUG,
            inputPath=actual_path,
            perEvent=False,
            outputEvent="hepmc3_event",
            numEvents=13,
        )
    )

    with ScopedFailureThreshold(acts.logging.MAX):
        with pytest.raises(RuntimeError) as excinfo:
            s.run()

    assert "Failed to process event" in str(excinfo.value)


def test_hepmc3_reader_per_event(tmp_path, rng):
    from acts.examples.hepmc3 import (
        HepMC3Writer,
        HepMC3Reader,
    )

    s = Sequencer(numThreads=1, events=12)

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
            inputEvent="hepmc3_event",
            outputParticles="particles_generated",
            outputVertices="vertices_generated",
        )
    )
    alg = AssertCollectionExistsAlg(
        [
            "particles_generated",
            "vertices_generated",
            "hepmc3_event",
        ],
        "check_alg",
        acts.logging.WARNING,
    )
    s.addAlgorithm(alg)

    out = tmp_path / "out" / "pytest.hepmc3"
    out.parent.mkdir(parents=True, exist_ok=True)

    s.addWriter(
        HepMC3Writer(
            acts.logging.DEBUG,
            inputEvent="hepmc3_event",
            outputPath=out,
            perEvent=True,
        )
    )

    s.run()

    assert alg.events_seen == 12

    s = Sequencer(numThreads=1)

    s.addReader(
        HepMC3Reader(
            acts.logging.DEBUG,
            inputPath=out,
            perEvent=True,
            outputEvent="hepmc3_event",
        )
    )

    s.run()


def test_hepmc3_reader_per_event_with_num_events(tmp_path, rng):
    from acts.examples.hepmc3 import (
        HepMC3Reader,
    )

    # Test that using both perEvent=True and numEvents raises an error
    with pytest.raises(ValueError) as excinfo:
        HepMC3Reader(
            acts.logging.DEBUG,
            inputPath="dummy.hepmc3",
            perEvent=True,
            outputEvent="hepmc3_event",
            numEvents=5,
        )

    assert "perEvent and numEvents are mutually exclusive" in str(excinfo.value)


def test_hepmc3_compression_modes():
    assert cm.none in acts.examples.hepmc3.availableCompressionModes()
