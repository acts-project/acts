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
