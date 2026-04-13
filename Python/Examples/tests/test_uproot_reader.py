# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import pytest
from pathlib import Path

from helpers import AssertCollectionExistsAlg

import acts
import acts.examples
from acts.examples.uproot_reader import UprootParticleReader, UprootSimHitReader
from acts.examples import Sequencer
from acts.examples.root import RootParticleWriter, RootSimHitWriter


@pytest.fixture
def sim_output(tmp_path, fatras, conf_const):
    """Run one simulation and return paths to the two ROOT files."""
    particle_file = tmp_path / "particles_simulation.root"
    hit_file = tmp_path / "hits.root"

    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.WARNING)
    evGen, simAlg, _ = fatras(s)

    s.addWriter(
        conf_const(
            RootParticleWriter,
            acts.logging.WARNING,
            inputParticles=simAlg.config.outputParticles,
            filePath=str(particle_file),
        )
    )
    s.addWriter(
        conf_const(
            RootSimHitWriter,
            acts.logging.WARNING,
            inputSimHits=simAlg.config.outputSimHits,
            filePath=str(hit_file),
        )
    )
    s.run()

    assert particle_file.exists()
    assert hit_file.exists()
    return particle_file, hit_file


@pytest.mark.root
def test_uproot_particle_reader(sim_output, conf_const):
    """Read back particles."""
    particle_file, _ = sim_output

    s = Sequencer(numThreads=1, logLevel=acts.logging.WARNING)
    s.addReader(
        UprootParticleReader(
            filePath=particle_file,
            outputParticles="particles_read",
            level=acts.logging.WARNING,
        )
    )
    alg = AssertCollectionExistsAlg(
        ["particles_read"],
        "check_alg",
        acts.logging.WARNING,
    )
    s.addAlgorithm(alg)
    s.run()
    assert alg.events_seen == 10


@pytest.mark.root
def test_uproot_sim_hit_reader(sim_output, conf_const):
    """Read back simhits."""
    _, hit_file = sim_output

    s = Sequencer(numThreads=1, logLevel=acts.logging.WARNING)
    s.addReader(
        UprootSimHitReader(
            filePath=hit_file,
            outputSimHits="simhits_read",
            level=acts.logging.WARNING,
        )
    )
    alg = AssertCollectionExistsAlg(
        ["simhits_read"],
        "check_alg",
        acts.logging.WARNING,
    )
    s.addAlgorithm(alg)
    s.run()
    assert alg.events_seen == 10


@pytest.mark.root
def test_root_write_uproot_read(tmp_path, fatras, conf_const):
    """Full round-trip: simulate with fatras, write ROOT files, read back with uproot readers."""
    particle_file = tmp_path / "particles_simulation.root"
    hit_file = tmp_path / "hits.root"

    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.WARNING)
    evGen, simAlg, _ = fatras(s)
    s.addWriter(
        conf_const(
            RootParticleWriter,
            acts.logging.WARNING,
            inputParticles=simAlg.config.outputParticles,
            filePath=str(particle_file),
        )
    )
    s.addWriter(
        conf_const(
            RootSimHitWriter,
            acts.logging.WARNING,
            inputSimHits=simAlg.config.outputSimHits,
            filePath=str(hit_file),
        )
    )
    s.run()

    s2 = Sequencer(numThreads=1, logLevel=acts.logging.WARNING)
    s2.addReader(
        UprootParticleReader(
            filePath=particle_file,
            outputParticles="particles_read",
            level=acts.logging.WARNING,
        )
    )
    s2.addReader(
        UprootSimHitReader(
            filePath=hit_file,
            outputSimHits="simhits_read",
            level=acts.logging.WARNING,
        )
    )
    alg = AssertCollectionExistsAlg(
        ["particles_read", "simhits_read"],
        "check_alg",
        acts.logging.WARNING,
    )
    s2.addAlgorithm(alg)
    s2.run()
    assert alg.events_seen == 10
