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


@pytest.mark.root
def test_uproot_reader(tmp_path, fatras, conf_const):
    """Write particles and simhits with ROOT writers, read back with UprootReader."""
    particle_file = tmp_path / "particles_simulation.root"
    hit_file = tmp_path / "hits.root"

    # Step 1: Run simulation and write ROOT output
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

    assert particle_file.exists(), "particle ROOT file was not created"
    assert hit_file.exists(), "simhit ROOT file was not created"

    # Step 2: Read back with UprootReader (no ROOT plugin needed)
    s2 = Sequencer(numThreads=1, logLevel=acts.logging.WARNING)

    s2.addReader(
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
    s2.addAlgorithm(alg)

    s2.run()

    assert alg.events_seen == 10
