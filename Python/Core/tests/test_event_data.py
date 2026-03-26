import acts


def test_space_point_container():
    """Test SpacePointContainer2 and ConstSpacePointProxy2 bindings (read-only)."""
    columns = (
        acts.SpacePointColumns.X | acts.SpacePointColumns.Y | acts.SpacePointColumns.Z
    )
    container = acts.SpacePointContainer2(columns)
    assert container.empty
    assert len(container) == 0

    # Mutable bindings removed - container is read-only from Python
    # Can still test iteration and len on empty container
    assert list(container) == []


def test_seed_container():
    """Test SeedContainer2 and ConstSeedProxy2 bindings (read-only)."""
    seed_container = acts.SeedContainer2()
    assert seed_container.empty
    assert len(seed_container) == 0

    # Mutable bindings removed - container is read-only from Python
    assert list(seed_container) == []


def test_particle_hypothesis():
    muon = acts.ParticleHypothesis.muon
    pion = acts.ParticleHypothesis.pion
    electron = acts.ParticleHypothesis.electron
    proton = acts.ParticleHypothesis.proton
    kaon = acts.ParticleHypothesis.kaon
    geantino = acts.ParticleHypothesis.geantino
    chargedGeantino = acts.ParticleHypothesis.chargedGeantino

    # create new particle hypothesis

    # check pdg
    assert muon.absolutePdg == acts.PdgParticle.eMuon
    assert pion.absolutePdg == acts.PdgParticle.ePionPlus
    assert electron.absolutePdg == acts.PdgParticle.eElectron
    assert kaon.absolutePdg == acts.PdgParticle.eKaonPlus
    assert proton.absolutePdg == acts.PdgParticle.eProton
    assert geantino.absolutePdg == acts.PdgParticle.eInvalid
    assert chargedGeantino.absolutePdg == acts.PdgParticle.eInvalid

    # check mass
    assert electron.mass != 0
    assert electron.mass < muon.mass
    assert muon.mass < pion.mass
    assert pion.mass < kaon.mass
    assert kaon.mass < proton.mass
    assert geantino.mass == 0
    assert chargedGeantino.mass == 0

    # check charge
    assert electron.absoluteCharge == 1
    assert muon.absoluteCharge == 1
    assert pion.absoluteCharge == 1
    assert proton.absoluteCharge == 1
    assert geantino.absoluteCharge == 0
    assert chargedGeantino.absoluteCharge == 1
    assert kaon.absoluteCharge == 1

    # printing should show something sensible
    assert str(muon) == "ParticleHypothesis{absPdg=mu, mass=0.105658, absCharge=1}"
    assert str(pion) == "ParticleHypothesis{absPdg=pi, mass=0.13957, absCharge=1}"
    assert (
        str(electron) == "ParticleHypothesis{absPdg=e, mass=0.000510999, absCharge=1}"
    )
    assert str(kaon) == "ParticleHypothesis{absPdg=K, mass=0.493677, absCharge=1}"
    assert str(proton) == "ParticleHypothesis{absPdg=p, mass=0.938272, absCharge=1}"
    assert str(geantino) == "ParticleHypothesis{absPdg=0, mass=0, absCharge=0}"
    assert str(chargedGeantino) == "ParticleHypothesis{absPdg=0, mass=0, absCharge=1}"
