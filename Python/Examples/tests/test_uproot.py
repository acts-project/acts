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
from acts.examples.uproot import UprootParticleReader, UprootSimHitReader
from acts.examples import Sequencer
from acts.examples.root import RootParticleWriter, RootSimHitWriter


@pytest.mark.root
def test_root_write_uproot_read(tmp_path, fatras, conf_const):
    """Full round-trip: simulate with fatras, write ROOT files, read back with uproot readers."""
    particle_file = tmp_path / "particles_simulation.root"
    hit_file = tmp_path / "hits.root"

    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.INFO)
    evGen, simAlg, _ = fatras(s)
    s.addWriter(
        conf_const(
            RootParticleWriter,
            acts.logging.INFO,
            inputParticles=simAlg.config.outputParticles,
            filePath=str(particle_file),
        )
    )
    s.addWriter(
        conf_const(
            RootSimHitWriter,
            acts.logging.INFO,
            inputSimHits=simAlg.config.outputSimHits,
            filePath=str(hit_file),
        )
    )
    s.run()

    s2 = Sequencer(numThreads=1, logLevel=acts.logging.INFO)
    s2.addReader(
        UprootParticleReader(
            filePath=particle_file,
            outputParticles="particles_read",
            level=acts.logging.INFO,
        )
    )
    s2.addReader(
        UprootSimHitReader(
            filePath=hit_file,
            outputSimHits="simhits_read",
            level=acts.logging.INFO,
        )
    )
    alg = AssertCollectionExistsAlg(
        ["particles_read", "simhits_read"],
        "check_alg",
        acts.logging.INFO,
    )
    s2.addAlgorithm(alg)
    s2.run()
    assert alg.events_seen == 10
