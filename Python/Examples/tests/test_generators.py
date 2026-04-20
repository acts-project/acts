import acts

u = acts.UnitConstants


def test_uniform_vertex_generator():
    vertex_generator = acts.examples.UniformVertexGenerator(
        min=acts.Vector4(
            -1 * u.um,
            -2 * u.um,
            -3 * u.mm,
            -4 * u.ns,
        ),
        max=acts.Vector4(
            1 * u.um,
            2 * u.um,
            3 * u.mm,
            4 * u.ns,
        ),
    )

    assert vertex_generator.min[0] == -1 * u.um
    assert vertex_generator.min[1] == -2 * u.um
    assert vertex_generator.min[2] == -3 * u.mm
    assert vertex_generator.min[3] == -4 * u.ns
    assert vertex_generator.max[0] == 1 * u.um
    assert vertex_generator.max[1] == 2 * u.um
    assert vertex_generator.max[2] == 3 * u.mm
    assert vertex_generator.max[3] == 4 * u.ns


def test_fixed_vertex_generator():
    rng = acts.examples.RandomEngine()
    gen = acts.examples.FixedVertexGenerator(
        fixed=acts.Vector4(1 * u.mm, 2 * u.mm, 3 * u.mm, 4 * u.ns)
    )
    pos = gen(rng, 0)
    assert pos[0] == 1 * u.mm
    assert pos[1] == 2 * u.mm
    assert pos[2] == 3 * u.mm
    assert pos[3] == 4 * u.ns


def test_gaussian_vertex_generator():
    rng = acts.examples.RandomEngine()
    gen = acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(1 * u.mm, 1 * u.mm, 1 * u.mm, 0),
        mean=acts.Vector4(0, 0, 0, 0),
    )
    # Generate several vertices and check they are not all identical
    positions = [gen(rng, i) for i in range(10)]
    xs = [p[0] for p in positions]
    assert len(set(xs)) > 1, "Gaussian generator should produce varying positions"


def test_additive_vertex_generator():
    rng = acts.examples.RandomEngine()
    g1 = acts.examples.FixedVertexGenerator(
        fixed=acts.Vector4(1 * u.mm, 2 * u.mm, 3 * u.mm, 0)
    )
    g2 = acts.examples.FixedVertexGenerator(
        fixed=acts.Vector4(10 * u.mm, 20 * u.mm, 30 * u.mm, 0)
    )
    additive = acts.examples.AdditiveVertexGenerator(generators=[g1, g2])

    assert len(additive.generators) == 2

    pos = additive(rng, 0)
    assert pos[0] == 11 * u.mm
    assert pos[1] == 22 * u.mm
    assert pos[2] == 33 * u.mm


def test_lumi_block_vertex_generator():
    rng = acts.examples.RandomEngine()
    gen = acts.examples.LumiBlockVertexGenerator(
        blockSize=1000,
        stddev=acts.Vector4(0.01 * u.mm, 0.01 * u.mm, 1.0 * u.mm, 0),
    )

    assert gen.blockSize == 1000
    assert gen.stddev[0] == 0.01 * u.mm
    assert gen.stddev[1] == 0.01 * u.mm
    assert gen.stddev[2] == 1.0 * u.mm
    assert gen.stddev[3] == 0

    # Events in the same block must produce the same offset
    pos_a = gen(rng, 0)
    pos_b = gen(rng, 999)
    for i in range(4):
        assert pos_a[i] == pos_b[i], f"Same block, component {i} should match"

    # Events in different blocks must produce different offsets
    pos_c = gen(rng, 1000)
    differs = any(pos_a[i] != pos_c[i] for i in range(3))
    assert differs, "Different blocks should produce different offsets"


def test_lumi_block_with_gaussian_composition():
    """Integration test: AdditiveVertexGenerator composing Gaussian + LumiBlock."""
    rng = acts.examples.RandomEngine()
    gaussian = acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0.01 * u.mm, 0.01 * u.mm, 50 * u.mm, 0),
        mean=acts.Vector4(0, 0, 0, 0),
    )
    lumi = acts.examples.LumiBlockVertexGenerator(
        blockSize=100,
        stddev=acts.Vector4(1 * u.mm, 1 * u.mm, 10 * u.mm, 0),
    )
    combined = acts.examples.AdditiveVertexGenerator(generators=[gaussian, lumi])

    # The lumi block contribution should be constant within a block,
    # so two calls within the same block should differ only by the gaussian part
    lumi_offset_0 = lumi(rng, 0)
    lumi_offset_50 = lumi(rng, 50)
    for i in range(4):
        assert lumi_offset_0[i] == lumi_offset_50[i]

    # The combined generator should produce varying results even within a block
    positions = [combined(rng, evt) for evt in range(10)]
    xs = [p[0] for p in positions]
    assert len(set(xs)) > 1, "Combined generator should vary within a block"


def test_lumi_block_rotation_vertex_generator():
    """Test that the rotation wrapper applies a consistent per-block tilt."""
    rng = acts.examples.RandomEngine()
    base = acts.examples.FixedVertexGenerator(fixed=acts.Vector4(0, 0, 100 * u.mm, 0))
    gen = acts.examples.LumiBlockRotationVertexGenerator(
        base=base,
        blockSize=1000,
        xAngleStddev=0.01,
        yAngleStddev=0.01,
    )

    # Same block must produce identical results
    pos_a = gen(rng, 0)
    pos_b = gen(rng, 500)
    for i in range(4):
        assert pos_a[i] == pos_b[i], f"Same block, component {i} should match"

    # Different block should produce a different rotation
    pos_c = gen(rng, 1000)
    differs = any(pos_a[i] != pos_c[i] for i in range(3))
    assert differs, "Different blocks should produce different rotations"

    # Time component must be unchanged
    assert pos_a[3] == 0


def test_lumi_block_rotation_preserves_length():
    """Rotation should preserve the 3D length of the position vector."""
    import math

    rng = acts.examples.RandomEngine()
    base = acts.examples.FixedVertexGenerator(
        fixed=acts.Vector4(1 * u.mm, 2 * u.mm, 100 * u.mm, 5 * u.ns)
    )
    gen = acts.examples.LumiBlockRotationVertexGenerator(
        base=base,
        blockSize=100,
        xAngleStddev=0.1,
        yAngleStddev=0.1,
    )

    original_len = math.sqrt((1 * u.mm) ** 2 + (2 * u.mm) ** 2 + (100 * u.mm) ** 2)
    for evt in [0, 100, 200, 999]:
        pos = gen(rng, evt)
        rotated_len = math.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2)
        assert abs(rotated_len - original_len) < 1e-10 * u.mm


def test_rotation_uncorrelated_with_position_shift():
    """Rotation and position shift use different seeds for the same block size."""
    rng = acts.examples.RandomEngine()
    base = acts.examples.FixedVertexGenerator(fixed=acts.Vector4(0, 0, 100 * u.mm, 0))
    rot_gen = acts.examples.LumiBlockRotationVertexGenerator(
        base=base,
        blockSize=100,
        xAngleStddev=0.01,
        yAngleStddev=0.01,
    )
    pos_gen = acts.examples.LumiBlockVertexGenerator(
        blockSize=100,
        stddev=acts.Vector4(1 * u.mm, 1 * u.mm, 1 * u.mm, 0),
    )

    # Collect which blocks change for rotation vs position shift
    rot_blocks_differ = []
    pos_blocks_differ = []
    ref_rot = rot_gen(rng, 0)
    ref_pos = pos_gen(rng, 0)
    for block in range(1, 10):
        evt = block * 100
        r = rot_gen(rng, evt)
        p = pos_gen(rng, evt)
        rot_blocks_differ.append(any(ref_rot[i] != r[i] for i in range(3)))
        pos_blocks_differ.append(any(ref_pos[i] != p[i] for i in range(3)))

    # Both should vary across blocks
    assert any(rot_blocks_differ), "Rotation should vary across blocks"
    assert any(pos_blocks_differ), "Position shift should vary across blocks"


def test_full_composition_tilt_shift_gaussian():
    """Full integration: Gaussian smearing + lumi-block tilt + lumi-block shift."""
    rng = acts.examples.RandomEngine()
    gaussian = acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0.01 * u.mm, 0.01 * u.mm, 50 * u.mm, 1 * u.ns),
        mean=acts.Vector4(0, 0, 0, 0),
    )
    tilted_gaussian = acts.examples.LumiBlockRotationVertexGenerator(
        base=gaussian,
        blockSize=1000,
        xAngleStddev=0.001,
        yAngleStddev=0.001,
    )
    lumi_shift = acts.examples.LumiBlockVertexGenerator(
        blockSize=1000,
        stddev=acts.Vector4(0.01 * u.mm, 0.01 * u.mm, 1.0 * u.mm, 0),
    )
    combined = acts.examples.AdditiveVertexGenerator(
        generators=[tilted_gaussian, lumi_shift]
    )

    # Should produce varying positions within a block (gaussian component)
    positions = [combined(rng, evt) for evt in range(10)]
    xs = [p[0] for p in positions]
    assert len(set(xs)) > 1

    # Shift component should be constant within the block
    shift_a = lumi_shift(rng, 0)
    shift_b = lumi_shift(rng, 999)
    for i in range(4):
        assert shift_a[i] == shift_b[i]


def test_random_engine_seed():
    """RandomEngine wrapper exposes the seed it was constructed with."""
    rng_default = acts.examples.RandomEngine()
    assert rng_default.seed() == 0

    rng_seeded = acts.examples.RandomEngine(seed=42)
    assert rng_seeded.seed() == 42


def test_lumi_block_seed_dependence():
    """Different seeds must produce different lumi-block offsets."""
    gen = acts.examples.LumiBlockVertexGenerator(
        blockSize=100,
        stddev=acts.Vector4(1 * u.mm, 1 * u.mm, 1 * u.mm, 0),
    )

    rng_a = acts.examples.RandomEngine(seed=1)
    rng_b = acts.examples.RandomEngine(seed=2)

    pos_a = gen(rng_a, 0)
    pos_b = gen(rng_b, 0)
    differs = any(pos_a[i] != pos_b[i] for i in range(3))
    assert differs, "Different seeds should produce different lumi-block offsets"


def test_lumi_block_seed_reproducibility():
    """Same seed + same block must always give the same offset."""
    gen = acts.examples.LumiBlockVertexGenerator(
        blockSize=100,
        stddev=acts.Vector4(1 * u.mm, 1 * u.mm, 1 * u.mm, 0),
    )

    rng1 = acts.examples.RandomEngine(seed=42)
    rng2 = acts.examples.RandomEngine(seed=42)

    pos1 = gen(rng1, 50)
    pos2 = gen(rng2, 50)
    for i in range(4):
        assert pos1[i] == pos2[i], f"Same seed, same block: component {i} should match"


def test_lumi_block_rotation_seed_dependence():
    """Different seeds must produce different rotations."""
    base = acts.examples.FixedVertexGenerator(fixed=acts.Vector4(0, 0, 100 * u.mm, 0))
    gen = acts.examples.LumiBlockRotationVertexGenerator(
        base=base,
        blockSize=100,
        xAngleStddev=0.01,
        yAngleStddev=0.01,
    )

    rng_a = acts.examples.RandomEngine(seed=1)
    rng_b = acts.examples.RandomEngine(seed=2)

    pos_a = gen(rng_a, 0)
    pos_b = gen(rng_b, 0)
    differs = any(pos_a[i] != pos_b[i] for i in range(3))
    assert differs, "Different seeds should produce different rotations"


def test_examples_fatras_aliases_present():
    assert hasattr(acts.examples, "SimBarcode")
    assert hasattr(acts.examples, "GenerationProcess")
    assert hasattr(acts.examples, "SimulationOutcome")
    assert hasattr(acts.examples, "SimParticleState")

    assert acts.examples.SimBarcode is acts.fatras.Barcode
    assert acts.examples.GenerationProcess is acts.fatras.GenerationProcess
    assert acts.examples.SimulationOutcome is acts.fatras.SimulationOutcome
    assert acts.examples.SimParticleState is acts.fatras.Particle
