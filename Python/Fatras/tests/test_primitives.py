import acts
import acts.fatras


def test_barcode_roundtrip_properties():
    barcode = acts.fatras.Barcode()
    assert not barcode.isValid()

    barcode.vertexPrimary = 1
    barcode.vertexSecondary = 2
    barcode.particle = 3
    barcode.generation = 4
    barcode.subParticle = 5

    assert barcode.vertexPrimary == 1
    assert barcode.vertexSecondary == 2
    assert barcode.particle == 3
    assert barcode.generation == 4
    assert barcode.subParticle == 5
    assert barcode.isValid()
    assert "vp=" in repr(barcode)


def test_invalid_barcode_factory():
    invalid = acts.fatras.Barcode.Invalid()
    assert isinstance(invalid, acts.fatras.Barcode)
    assert not invalid.isValid()


def test_process_and_outcome_enums_available():
    assert (
        acts.fatras.GenerationProcess.eUndefined != acts.fatras.GenerationProcess.eDecay
    )
    assert (
        acts.fatras.GenerationProcess.ePhotonConversion
        == acts.fatras.GenerationProcess.ePhotonConversion
    )
    assert (
        acts.fatras.SimulationOutcome.Alive
        != acts.fatras.SimulationOutcome.KilledInteraction
    )


def test_particle_construction_and_properties():
    barcode = acts.fatras.Barcode()
    barcode.vertexPrimary = 1
    barcode.particle = 42

    particle = acts.fatras.Particle(barcode, acts.PdgParticle.eMuon, -1.0, 0.105658)

    assert isinstance(particle.particleId, acts.fatras.Barcode)
    assert particle.particleId.vertexPrimary == 1
    assert particle.particleId.particle == 42
    assert particle.pdg == acts.PdgParticle.eMuon
    assert particle.absolutePdg == acts.PdgParticle.eMuon
    assert particle.charge == -1.0
    assert particle.mass == 0.105658


def test_particle_construction_from_pdg_table():
    barcode = acts.fatras.Barcode()
    barcode.vertexPrimary = 1
    barcode.particle = 7

    particle = acts.fatras.Particle(barcode, acts.PdgParticle.eMuon)
    assert particle.particleId.particle == 7
    assert particle.pdg == acts.PdgParticle.eMuon
    assert particle.absolutePdg == acts.PdgParticle.eMuon
    assert particle.mass > 0
