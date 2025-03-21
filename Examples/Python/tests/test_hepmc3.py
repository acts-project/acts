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


@pytest.fixture(scope="session")
def hepmc_data_impl(tmp_path_factory):
    import subprocess

    script = (
        Path(__file__).parent.parent.parent.parent
        / "Examples"
        / "Scripts"
        / "Python"
        / "event_recording.py"
    )
    assert script.exists()

    with tempfile.TemporaryDirectory() as tmp_path:
        env = os.environ.copy()
        env["NEVENTS"] = "1"
        subprocess.check_call([sys.executable, str(script)], cwd=tmp_path, env=env)

        outfile = Path(tmp_path) / "hepmc3/event000000000-events.hepmc3"
        # fake = Path("/scratch/pagessin/acts/hepmc3/event000000000-events.hepmc3")

        # outfile.parent.mkdir()
        # shutil.copy(fake, outfile)

        assert outfile.exists()

        yield outfile


@pytest.fixture
def hepmc_data(hepmc_data_impl: Path, tmp_path):
    dest = tmp_path / hepmc_data_impl.name
    shutil.copy(hepmc_data_impl, dest)

    return dest


@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
@pytest.mark.odd
@pytest.mark.slow
def test_hepmc3_histogram(hepmc_data, tmp_path):
    from acts.examples.hepmc3 import (
        HepMC3AsciiReader,
        HepMCProcessExtractor,
    )

    s = Sequencer(numThreads=1)

    s.addReader(
        HepMC3AsciiReader(
            level=acts.logging.INFO,
            inputDir=str(hepmc_data.parent),
            inputStem="events",
            outputEvents="hepmc-events",
        )
    )

    s.addAlgorithm(
        HepMCProcessExtractor(
            level=acts.logging.INFO,
            inputEvents="hepmc-events",
            extractionProcess="Inelastic",
        )
    )

    # This segfaults, see https://github.com/acts-project/acts/issues/914
    # s.addWriter(
    #     RootNuclearInteractionParametersWriter(
    #         level=acts.logging.INFO, inputSimulationProcesses="event-fraction"
    #     )
    # )

    alg = AssertCollectionExistsAlg(
        "hepmc-events", name="check_alg", level=acts.logging.INFO
    )
    s.addAlgorithm(alg)

    s.run()


@pytest.mark.parametrize("per_event", [True, False])
def test_hepmc3_particle_writer(tmp_path, rng, per_event):
    from acts.examples.hepmc3 import (
        HepMC3AsciiWriter,
        HepMC3OutputConverter,
    )

    s = Sequencer(numThreads=4, events=10)

    evGen = acts.examples.EventGenerator(
        level=acts.logging.INFO,
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=4),
                vertex=acts.examples.GaussianVertexGenerator(
                    stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
                ),
                particles=acts.examples.ParametricParticleGenerator(
                    p=(100 * u.GeV, 100 * u.GeV),
                    eta=(-2, 2),
                    phi=(0, 360 * u.degree),
                    randomizeCharge=True,
                    numParticles=10,
                ),
            )
        ],
        outputParticles="particles_generated",
        outputVertices="vertices_input",
        randomNumbers=rng,
    )

    s.addReader(evGen)

    out = tmp_path / "out" / "events_pytest.hepmc3"
    out.parent.mkdir(parents=True, exist_ok=True)

    s.addAlgorithm(
        HepMC3OutputConverter(
            level=acts.logging.VERBOSE,
            inputParticles=evGen.config.outputParticles,
            inputVertices=evGen.config.outputVertices,
            outputEvents="hepmc-events",
        )
    )
    s.addWriter(
        HepMC3AsciiWriter(
            acts.logging.VERBOSE,
            inputEvents="hepmc-events",
            outputPath=out,
            perEvent=per_event,
        )
    )

    s.run()

    if per_event:
        files = list(out.parent.iterdir())
        assert len(files) == 10
        assert all(f.suffix == ".hepmc3" for f in files)
        assert all(f.stem.startswith("events_pytest") for f in files)
        for f in files:
            with f.open("r") as f:
                assert len(f.readlines()) == 46
    else:
        assert out.exists(), f"File {out} does not exist"
        with out.open("r") as f:
            assert len(f.readlines()) == 424


def test_hepmc3_particle_writer_pythia8(tmp_path, rng):
    from acts.examples.hepmc3 import (
        HepMC3AsciiWriter,
        HepMC3OutputConverter,
    )
    from pythia8 import addPythia8

    s = Sequencer(numThreads=1, events=10)

    addPythia8(
        s,
        rnd=rng,
        hardProcess=["Top:qqbar2ttbar=on"],
        npileup=0,
        outputDirCsv=None,
        outputDirRoot=None,
    )

    out = tmp_path / "events.hepmc3"

    s.addAlgorithm(
        HepMC3OutputConverter(
            level=acts.logging.VERBOSE,
            inputParticles="particles_generated",
            inputVertices="vertices_generated",
            outputEvents="hepmc-events",
        )
    )
    s.addWriter(
        HepMC3AsciiWriter(
            acts.logging.VERBOSE,
            inputEvents="hepmc-events",
            outputPath=out,
        )
    )

    s.run()

    assert out.exists(), f"File {out} does not exist"

    import shutil

    shutil.copy(out, Path.cwd() / "events.hepmc3")

    # with out.open("r") as f:
    #     assert len(f.readlines()) == 326237
