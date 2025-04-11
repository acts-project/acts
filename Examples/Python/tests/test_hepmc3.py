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


@pytest.mark.parametrize("per_event", [True, False], ids=["per_event", "combined"])
def test_hepmc3_particle_writer(tmp_path, rng, per_event):
    from acts.examples.hepmc3 import (
        HepMC3Writer,
        HepMC3OutputConverter,
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
        outputParticles="particles_generated",
        outputVertices="vertices_input",
        outputEvent="hepmc3_event",
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
                assert len(f.readlines()) == 25
    else:
        assert out.exists(), f"File {out} does not exist"
        with out.open("r") as f:
            assert len(f.readlines()) == 214

        try:
            import pyhepmc
            from pyhepmc.view import to_dot

            nevts = 0
            with pyhepmc.open(out) as f:
                nevts += 1
                for evt in f:
                    assert len(evt.particles) == 4 + 2  # muons + beam particles
            assert nevts == 1

        except ImportError:
            pass


def test_hepmc3_particle_writer_pythia8(tmp_path, rng):
    from acts.examples.hepmc3 import (
        HepMC3Writer,
        HepMC3OutputConverter,
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
        outputEvent="hepmc3_event",
        logLevel=acts.logging.INFO,
        vtxGen=vtxGen,
    )

    out = tmp_path / "events.hepmc3"

    # s.addAlgorithm(
    #     HepMC3OutputConverter(
    #         level=acts.logging.VERBOSE,
    #         inputParticles="particles_generated",
    #         inputVertices="vertices_generated",
    #         outputEvents="hepmc-events",
    #     )
    # )

    s.addWriter(
        HepMC3Writer(
            acts.logging.VERBOSE,
            inputEvent="hepmc3_event",
            outputPath=out,
            perEvent=False,
        )
    )

    s.run()

    assert out.exists(), f"File {out} does not exist"
    with out.open("r") as f:
        assert len(f.readlines()) == 18679

    try:
        import pyhepmc

        nevts = 0
        with pyhepmc.open(out) as f:
            nevts += 1
            for evt in f:
                assert len(evt.particles) == 7433
        assert nevts == 1

    except ImportError:
        pass
