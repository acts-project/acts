import acts
import pytest

def test_pgd_particle():
    assert len(acts.PdgParticle.__members__) == 19

def test_algebra():
    # Vector2 testing
    v2 = acts.Vector2(1, 2)
    v2 = acts.Vector2([1, 2])
    assert v2[0] == 1
    assert v2[1] == 2
    print(str(v2))
    with pytest.raises(TypeError):
        acts.Vector2(1, 2, 3)

    # Vector3 testing
    v3 = acts.Vector3(1, 2, 3)
    with pytest.raises(TypeError):
        acts.Vector3(1, 2, 3, 4)
    with pytest.raises(TypeError):
        acts.Vector3(1, 2)
    assert v3[0] == 1
    assert v3[1] == 2
    assert v3[2] == 3
    print(str(v3))

    v3 = acts.Vector3([1, 2, 3])
    with pytest.raises(TypeError):
        acts.Vector3([1, 2, 3, 4])
    with pytest.raises(TypeError):
        acts.Vector3([1, 2])
    with pytest.raises(TypeError):
        acts.Vector3()

    vX = acts.Vector3.UnitX()
    vY = acts.Vector3.UnitY()
    vZ = acts.Vector3.UnitZ()
    assert vX[0] == 1 and vX[1] == 0 and vX[2] == 0
    assert vY[0] == 0 and vY[1] == 1 and vY[2] == 0
    assert vZ[0] == 0 and vZ[1] == 0 and vZ[2] == 1

    # Vector4 testing
    v4 = acts.Vector4(1, 2, 3, 4)
    assert v4[0] == 1
    assert v4[1] == 2
    assert v4[2] == 3
    assert v4[3] == 4
    print(str(v4))
    with pytest.raises(TypeError):
        v4 = acts.Vector4(1, 2, 3)
    v4 = acts.Vector4([1, 2, 3, 4])
    with pytest.raises(TypeError):
        acts.Vector4([1, 2, 3])
    with pytest.raises(TypeError):
        acts.Vector4()

    # AngleAxis3 testing
    aa = acts.AngleAxis3(0.5, acts.Vector3(1, 0, 0))
    print(str(aa))

    # Translation3 testing
    t3 = acts.Translation3(acts.Vector3(1, 2, 3))
    assert t3[0] == 1
    assert t3[1] == 2
    assert t3[2] == 3
    print(str(t3))

    # Transform3 testing
    t3 = acts.Transform3(acts.Vector3(1, 2, 3))

def test_units() :
    assert len(dir(acts.UnitConstants)) == 44