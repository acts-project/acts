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


def test_examples_fatras_aliases_present():
    assert hasattr(acts.examples, "SimBarcode")
    assert hasattr(acts.examples, "GenerationProcess")
    assert hasattr(acts.examples, "SimulationOutcome")
    assert hasattr(acts.examples, "SimParticleState")

    assert acts.examples.SimBarcode is acts.fatras.Barcode
    assert acts.examples.GenerationProcess is acts.fatras.GenerationProcess
    assert acts.examples.SimulationOutcome is acts.fatras.SimulationOutcome
    assert acts.examples.SimParticleState is acts.fatras.Particle
