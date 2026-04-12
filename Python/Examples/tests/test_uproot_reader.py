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
from acts.examples.uproot_reader import UprootReader
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
def test_uproot_reader_both(sim_output, conf_const):
    """Read back particles and simhits together."""
    particle_file, hit_file = sim_output

    s = Sequencer(numThreads=1, logLevel=acts.logging.WARNING)
    s.addReader(
        UprootReader(
            particleFilePath=particle_file,
            simHitFilePath=hit_file,
            outputParticles="particles_read",
            outputSimHits="simhits_read",
            level=acts.logging.WARNING,
        )
    )
    alg = AssertCollectionExistsAlg(
        ["particles_read", "simhits_read"],
        "check_alg",
        acts.logging.WARNING,
    )
    s.addAlgorithm(alg)
    s.run()
    assert alg.events_seen == 10


@pytest.mark.root
def test_uproot_reader_particles_only(sim_output, conf_const):
    """Read back particles without a simhit file."""
    particle_file, _ = sim_output

    s = Sequencer(numThreads=1, logLevel=acts.logging.WARNING)
    s.addReader(
        UprootReader(
            particleFilePath=particle_file,
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
def test_uproot_reader_hits_only(sim_output, conf_const):
    """Read back simhits without a particle file."""
    _, hit_file = sim_output

    s = Sequencer(numThreads=1, logLevel=acts.logging.WARNING)
    s.addReader(
        UprootReader(
            simHitFilePath=hit_file,
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
