import gc
import pytest
import acts


def test_space_point_container():
    """Test SpacePointContainer2 and MutableSpacePointProxy2 bindings."""
    columns = (
        acts.SpacePointColumns.X | acts.SpacePointColumns.Y | acts.SpacePointColumns.Z
    )
    container = acts.SpacePointContainer2(columns)
    assert container.empty
    assert len(container) == 0
    assert list(container) == []

    # createSpacePoint returns a MutableSpacePointProxy2
    sp0 = container.createSpacePoint()
    assert isinstance(sp0, acts.MutableSpacePointProxy2)
    assert not container.empty
    assert len(container) == 1
    assert sp0.index == 0

    # Write and read back scalar fields
    sp0.x = 1.0
    sp0.y = 2.0
    sp0.z = 3.0
    assert sp0.x == pytest.approx(1.0)
    assert sp0.y == pytest.approx(2.0)
    assert sp0.z == pytest.approx(3.0)

    # __getitem__ returns MutableSpacePointProxy2
    via_index = container[0]
    assert isinstance(via_index, acts.MutableSpacePointProxy2)
    assert via_index.x == pytest.approx(1.0)

    # Mutation via __getitem__ proxy is reflected in the container
    via_index.x = 99.0
    assert sp0.x == pytest.approx(99.0)

    # __iter__ yields MutableSpacePointProxy2
    sp1 = container.createSpacePoint()
    sp1.x = 5.0
    items = list(container)
    assert len(items) == 2
    assert isinstance(items[0], acts.MutableSpacePointProxy2)
    assert items[0].x == pytest.approx(99.0)
    assert items[1].x == pytest.approx(5.0)


def test_space_point_all_columns_round_trip():
    """Read/write every guarded column when the columns are present."""
    Cols = acts.SpacePointColumns
    columns = (
        Cols.X
        | Cols.Y
        | Cols.Z
        | Cols.R
        | Cols.Phi
        | Cols.Time
        | Cols.VarianceZ
        | Cols.VarianceR
        | Cols.VarianceT
        | Cols.CopiedFromIndex
    )
    container = acts.SpacePointContainer2(columns)
    sp = container.createSpacePoint()

    # Scalar fields, including varianceT (writeable in C++ but previously not
    # exposed on the mutable proxy).
    sp.x = 1.0
    sp.y = 2.0
    sp.z = 3.0
    sp.r = 4.0
    sp.phi = 0.5
    sp.time = 6.0
    sp.varianceZ = 7.0
    sp.varianceR = 8.0
    sp.varianceT = 9.0
    assert sp.x == pytest.approx(1.0)
    assert sp.time == pytest.approx(6.0)
    assert sp.varianceT == pytest.approx(9.0)

    sp.copiedFromIndex = 0
    assert sp.copiedFromIndex == 0

    # Container-level numpy column views work when the column is present.
    assert container.x.shape == (1,)
    assert container.varianceT.shape == (1,)

    # Const proxy reads the same values.
    const_sp = container[0]
    assert const_sp.varianceT == pytest.approx(9.0)


def test_space_point_missing_column_raises():
    """Accessing a column that was not requested raises AttributeError
    instead of segfaulting on a disengaged optional column."""
    Cols = acts.SpacePointColumns

    # Proxy read of a column absent from an otherwise-populated container.
    container = acts.SpacePointContainer2(Cols.X)
    sp = container.createSpacePoint()
    with pytest.raises(AttributeError, match="time"):
        _ = sp.time

    # Proxy read on a fully empty (None) container: even x is absent.
    none_container = acts.SpacePointContainer2()
    none_sp = none_container.createSpacePoint()
    with pytest.raises(AttributeError, match="x"):
        _ = none_sp.x

    # Proxy write to an absent scalar column.
    with pytest.raises(AttributeError, match="time"):
        sp.time = 1.0

    # Proxy write to an absent array column.
    with pytest.raises(AttributeError, match="topStripVector"):
        sp.topStripVector = [0.0, 0.0, 0.0]

    # sourceLinks read without the SourceLinks column.
    with pytest.raises(AttributeError, match="sourceLinks"):
        _ = sp.sourceLinks

    # Const proxy read of an absent column.
    const_sp = container[0]
    with pytest.raises(AttributeError, match="varianceZ"):
        _ = const_sp.varianceZ

    # Container-level numpy column property of an absent column.
    with pytest.raises(AttributeError, match="time"):
        _ = container.time


def test_seed_container():
    """Test SeedContainer2, MutableSeedProxy2 and ConstSeedProxy2 bindings."""
    seed_container = acts.SeedContainer2()
    assert seed_container.empty
    assert len(seed_container) == 0
    assert list(seed_container) == []

    # createSeed returns a MutableSeedProxy2
    seed0 = seed_container.createSeed()
    assert isinstance(seed0, acts.MutableSeedProxy2)
    assert not seed_container.empty
    assert len(seed_container) == 1

    # Newly created seed is empty and has default values
    assert seed0.empty
    assert seed0.index == 0
    assert seed0.size == 0

    # Write quality and vertexZ
    seed0.quality = 0.9
    seed0.vertexZ = 12.5
    assert seed0.quality == pytest.approx(0.9)
    assert seed0.vertexZ == pytest.approx(12.5)

    # Assign space point indices and verify size
    seed0.assignSpacePointIndices([0, 1, 2])
    assert not seed0.empty
    assert seed0.size == 3

    # Create a second seed
    seed1 = seed_container.createSeed()
    assert seed1.index == 1
    assert len(seed_container) == 2

    # Read back via ConstSeedProxy2 (__getitem__)
    const0 = seed_container[0]
    assert isinstance(const0, acts.ConstSeedProxy2)
    assert const0.quality == pytest.approx(0.9)
    assert const0.vertexZ == pytest.approx(12.5)
    assert const0.size == 3
    assert list(const0.spacePointIndices) == [0, 1, 2]

    # Iteration yields ConstSeedProxy2 objects
    seeds = list(seed_container)
    assert len(seeds) == 2
    assert isinstance(seeds[0], acts.ConstSeedProxy2)

    # assignSpacePointContainer links the two containers (shared ownership).
    # After the call, spacePoints() on a seed proxy is resolvable.
    sp_container = acts.SpacePointContainer2(
        acts.SpacePointColumns.X | acts.SpacePointColumns.Y | acts.SpacePointColumns.Z
    )
    seed_container.assignSpacePointContainer(sp_container)
    del sp_container
    gc.collect()
    assert list(const0.spacePointIndices) == [0, 1, 2]


def test_space_point_proxy_fails_loud_after_disown():
    """A space point proxy/iterator whose container was transferred away (the
    whiteboard takes it as a unique_ptr, disowning the Python wrapper) must raise
    ValueError on access instead of dereferencing freed memory (SIGSEGV)."""
    tb = acts._testing

    Cols = acts.SpacePointColumns
    container = acts.SpacePointContainer2(Cols.X | Cols.Y | Cols.Z)
    sp = container.createSpacePoint()
    sp.x = 1.0
    via_index = container[0]
    it = iter(container)

    # Positive control: everything works while the container is owned.
    assert sp.x == pytest.approx(1.0)
    assert via_index.x == pytest.approx(1.0)

    # Transfer the container away (same mechanism as a whiteboard write).
    tb.consume_spacepoints(container)

    # Every previously-obtained handle now fails loud rather than segfaulting.
    with pytest.raises(ValueError, match="consumed"):
        _ = sp.x
    with pytest.raises(ValueError, match="consumed"):
        _ = via_index.x
    with pytest.raises(ValueError, match="consumed"):
        next(it)


def test_seed_proxy_fails_loud_after_disown():
    """Same fail-loud contract for seed proxies/iterators after the container is
    disowned by a unique_ptr transfer."""
    tb = acts._testing

    container = acts.SeedContainer2()
    seed = container.createSeed()
    seed.quality = 0.9
    via_index = container[0]
    it = iter(container)

    assert seed.quality == pytest.approx(0.9)
    assert via_index.quality == pytest.approx(0.9)

    tb.consume_seeds(container)

    with pytest.raises(ValueError, match="consumed"):
        _ = seed.quality
    with pytest.raises(ValueError, match="consumed"):
        _ = via_index.quality
    with pytest.raises(ValueError, match="consumed"):
        next(it)


def test_proxy_survives_container_gc():
    """The tether holds a strong reference to its container, so a proxy keeps
    working after the local container handle is dropped and collected (the
    container is only truly gone once transferred, not merely unreferenced)."""
    Cols = acts.SpacePointColumns
    sp_container = acts.SpacePointContainer2(Cols.X)
    sp = sp_container.createSpacePoint()
    sp.x = 7.0
    del sp_container
    gc.collect()
    assert sp.x == pytest.approx(7.0)
    assert isinstance(sp, acts.MutableSpacePointProxy2)

    seed_container = acts.SeedContainer2()
    seed = seed_container.createSeed()
    seed.quality = 3.0
    del seed_container
    gc.collect()
    assert seed.quality == pytest.approx(3.0)
    assert isinstance(seed, acts.MutableSeedProxy2)


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
