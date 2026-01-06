import pytest
import acts


def test_algebra():

    # Vector2 tests
    v2 = acts.Vector2(1, 2)
    with pytest.raises(TypeError):
        acts.Vector2(1)
    with pytest.raises(TypeError):
        acts.Vector2(1, 2, 3)

    v2 = acts.Vector2([1, 2]) + acts.Vector2([3, 4])
    assert v2[0] == 4 and v2[1] == 6
    with pytest.raises(TypeError):
        acts.Vector2([1, 2]) + acts.Vector3([1, 2, 3])

    v2 = acts.Vector2([5, 7]) - acts.Vector2([2, 3])
    assert v2[0] == 3 and v2[1] == 4
    with pytest.raises(TypeError):
        acts.Vector2([1, 2]) - acts.Vector3([1, 2, 3])

    v2 = acts.Vector2([2, 3]) * 3
    assert v2[0] == 6 and v2[1] == 9

    s2 = acts.Vector2([8, 4]) * acts.Vector2([0.5, 0.25])
    assert s2[0] == 4 and s2[1] == 1

    # Vector3 tests
    v3 = acts.Vector3(1, 2, 3)
    with pytest.raises(TypeError):
        acts.Vector3(1, 2, 3, 4)
    with pytest.raises(TypeError):
        acts.Vector3(1, 2)

    v3 = acts.Vector3([1, 2, 3])
    with pytest.raises(TypeError):
        acts.Vector3([1, 2, 3, 4])
    with pytest.raises(TypeError):
        acts.Vector3([1, 2])
    with pytest.raises(TypeError):
        acts.Vector3()

    unit_x = acts.Vector3.UnitX()
    assert unit_x[0] == 1 and unit_x[1] == 0 and unit_x[2] == 0

    unit_y = acts.Vector3.UnitY()
    assert unit_y[0] == 0 and unit_y[1] == 1 and unit_y[2] == 0

    unit_z = acts.Vector3.UnitZ()
    assert unit_z[0] == 0 and unit_z[1] == 0 and unit_z[2] == 1

    v3 = acts.Vector3([1, 0, 0]) + acts.Vector3([1, 2, 3])
    assert v3[0] == 2 and v3[1] == 2 and v3[2] == 3

    v3 = acts.Vector3([5, 7, 9]) - acts.Vector3([2, 3, 4])
    assert v3[0] == 3 and v3[1] == 4 and v3[2] == 5

    v3 = acts.Vector3([2, 3, 4]) * 2
    assert v3[0] == 4 and v3[1] == 6 and v3[2] == 8

    # Vector4 tests
    v4 = acts.Vector4(1, 2, 3, 4)
    with pytest.raises(TypeError):
        v4 = acts.Vector4(1, 2, 3)
    v4 = acts.Vector4([1, 2, 3, 4])
    with pytest.raises(TypeError):
        acts.Vector4([1, 2, 3])
    with pytest.raises(TypeError):
        acts.Vector4()

    v4 = acts.Vector4([1, 2, 3, 4]) + acts.Vector4([4, 3, 2, 1])
    assert v4[0] == 5 and v4[1] == 5 and v4[2] == 5 and v4[3] == 5

    v4 = acts.Vector4([5, 7, 9, 11]) - acts.Vector4([1, 2, 3, 4])
    assert v4[0] == 4 and v4[1] == 5 and v4[2] == 6 and v4[3] == 7

    v4 = acts.Vector4([2, 3, 4, 5]) * 3
    assert v4[0] == 6 and v4[1] == 9 and v4[2] == 12 and v4[3] == 15


def test_pgd_particle():
    assert len(acts.PdgParticle.__members__) == 32


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
    assert muon.absolutePdg() == acts.PdgParticle.eMuon
    assert pion.absolutePdg() == acts.PdgParticle.ePionPlus
    assert electron.absolutePdg() == acts.PdgParticle.eElectron
    assert kaon.absolutePdg() == acts.PdgParticle.eKaonPlus
    assert proton.absolutePdg() == acts.PdgParticle.eProton
    assert geantino.absolutePdg() == acts.PdgParticle.eInvalid
    assert chargedGeantino.absolutePdg() == acts.PdgParticle.eInvalid

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
    assert str(muon) == "ParticleHypothesis{absPdg=mu, mass=0.105658, absCharge=1}"
    assert str(pion) == "ParticleHypothesis{absPdg=pi, mass=0.13957, absCharge=1}"
    assert (
        str(electron) == "ParticleHypothesis{absPdg=e, mass=0.000510999, absCharge=1}"
    )
    assert str(kaon) == "ParticleHypothesis{absPdg=K, mass=0.493677, absCharge=1}"
    assert str(proton) == "ParticleHypothesis{absPdg=p, mass=0.938272, absCharge=1}"
    assert str(geantino) == "ParticleHypothesis{absPdg=0, mass=0, absCharge=0}"
    assert str(chargedGeantino) == "ParticleHypothesis{absPdg=0, mass=0, absCharge=1}"
