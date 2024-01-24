import pytest
from pathlib import Path

import acts


def test_particle_hypothesis():
    muon = acts.ParticleHypothesis.muon
    pion = acts.ParticleHypothesis.pion
    electron = acts.ParticleHypothesis.electron
    proton = acts.ParticleHypothesis.proton
    geantino = acts.ParticleHypothesis.geantino
    chargedGeantino = acts.ParticleHypothesis.chargedGeantino

    # create new particle hypothesis
    kaon = acts.ParticleHypothesis(321, 0.493677, 1)

    # check pdg
    assert muon.pdg() == acts.PdgParticle.muon
    assert pion.pdg() == acts.PdgParticle.pion
    assert electron.pdg() == acts.PdgParticle.electron
    assert proton.pdg() == acts.PdgParticle.proton
    assert geantino.pdg() == acts.PdgParticle.geantino
    assert chargedGeantino.pdg() == acts.PdgParticle.chargedGeantino
    assert kaon.pdg() == 321

    # check mass
    assert electron.mass() != 0
    assert electron.mass() < muon.mass()
    assert muon.mass() < pion.mass()
    assert pion.mass() < kaon.mass()
    assert kaon.mass() < proton.mass()
    assert geantino.mass() == 0
    assert chargedGeantino.mass() == 0

    # check charge
    assert electron.absoluteCharge() == 1
    assert muon.absoluteCharge() == 1
    assert pion.absoluteCharge() == 1
    assert proton.absoluteCharge() == 1
    assert geantino.absoluteCharge() == 0
    assert chargedGeantino.absoluteCharge() == 1
    assert kaon.absoluteCharge() == 1

    # printing should show something sensible
    print(muon)
    print(pion)
    print(electron)
    print(proton)
    print(geantino)
    print(chargedGeantino)
