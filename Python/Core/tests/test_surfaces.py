import pytest

import acts


def test_surface_bounds_base_api():
    bounds = acts.RectangleBounds(10.0, 5.0)

    assert bounds.type() == acts.SurfaceBoundsType.Rectangle
    assert bounds.isCartesian() is True
    assert bounds.values() == [10.0, 5.0]
    assert bounds.inside(acts.Vector2(0.0, 0.0)) is True
    assert bounds.inside(acts.Vector2(100.0, 100.0)) is False
    assert bounds.distance(acts.Vector2(0.0, 0.0)) == pytest.approx(5.0)
    center = bounds.center()
    assert center[0] == pytest.approx(0.0)
    assert center[1] == pytest.approx(0.0)
    assert "RectangleBounds" in str(bounds)


@pytest.mark.parametrize(
    "bounds, expected_size",
    [
        (acts.CylinderBounds(10.0, 20.0), 4),
        (acts.AnnulusBounds(10.0, 20.0, -0.2, 0.2), 6),
        (acts.RadialBounds(10.0, 20.0), 4),
        (acts.LineBounds(1.0, 100.0), 2),
        (acts.RectangleBounds(10.0, 20.0), 2),
        (acts.TrapezoidBounds(8.0, 12.0, 20.0, 0.1), 4),
    ],
)
def test_surface_bounds_indexing(bounds, expected_size):
    assert len(bounds) == expected_size
    for i in range(expected_size):
        assert bounds[i] == pytest.approx(bounds.values()[i])

    with pytest.raises(IndexError):
        _ = bounds[expected_size]

    with pytest.raises(IndexError):
        _ = bounds[-1]


def test_surface_factory_and_surface_api():
    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
    transform = acts.Transform3.Identity()

    plane_bounds = acts.RectangleBounds(10.0, 5.0)
    plane = acts.Surface.createPlane(transform, plane_bounds)
    assert isinstance(plane, acts.PlaneSurface)
    assert plane.type == acts.SurfaceType.Plane
    assert plane.bounds == plane_bounds
    assert plane.thickness == pytest.approx(0.0)
    assert isinstance(plane.isSensitive, bool)
    assert isinstance(plane.isAlignable, bool)

    local = acts.Vector2(1.0, 1.0)
    direction = acts.Vector3(0.0, 0.0, 1.0)
    global_position = plane.localToGlobal(gctx, local, direction)
    assert plane.insideBounds(local) is True
    assert plane.isOnSurface(gctx, global_position, direction) is True
    assert "PlaneSurface" in plane.name
    assert "PlaneSurface" in plane.toString(gctx)

    cylinder = acts.Surface.createCylinder(transform, acts.CylinderBounds(20.0, 50.0))
    assert isinstance(cylinder, acts.CylinderSurface)
    assert cylinder.type == acts.SurfaceType.Cylinder

    perigee = acts.Surface.createPerigee(acts.Vector3(0.0, 0.0, 0.0))
    assert isinstance(perigee, acts.PerigeeSurface)
    assert perigee.type == acts.SurfaceType.Perigee
