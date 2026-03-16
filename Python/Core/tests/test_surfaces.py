import pytest

import acts


def test_surface_bounds_base_api():
    bounds = acts.RectangleBounds(10.0, 5.0)

    assert bounds.type == acts.SurfaceBoundsType.Rectangle
    assert bounds.isCartesian() is True
    assert bounds.values() == [-10.0, -5.0, 10.0, 5.0]
    assert bounds.inside(acts.Vector2(0.0, 0.0)) is True
    assert bounds.inside(acts.Vector2(100.0, 100.0)) is False
    assert bounds.distance(acts.Vector2(0.0, 0.0)) == pytest.approx(5.0)
    center = bounds.center()
    assert center[0] == pytest.approx(0.0)
    assert center[1] == pytest.approx(0.0)
    assert "RectangleBounds" in str(bounds)


def test_boundary_tolerance_binding_and_inside_overload():
    none_tol = acts.BoundaryTolerance.none()
    inf_tol = acts.BoundaryTolerance.infinite()
    abs_tol = acts.BoundaryTolerance.absoluteEuclidean(0.5)

    assert none_tol.isNone() is True
    assert none_tol.isInfinite() is False

    assert inf_tol.isInfinite() is True
    assert inf_tol.isNone() is False

    assert abs_tol.hasAbsoluteEuclidean() is True
    assert abs_tol.hasChi2Bound() is False
    assert abs_tol.hasChi2Cartesian() is False

    bounds = acts.RectangleBounds(10.0, 5.0)
    assert bounds.inside(acts.Vector2(0.0, 0.0), none_tol) is True
    assert bounds.inside(acts.Vector2(100.0, 100.0), inf_tol) is True


def test_bound_value_enums_exposed():
    assert acts.CylinderBoundsValue.R is not None
    assert acts.AnnulusBoundsValue.MinR is not None
    assert acts.RadialBoundsValue.MinR is not None
    assert acts.LineBoundsValue.R is not None
    assert acts.RectangleBoundsValue.MinX is not None
    assert acts.TrapezoidBoundsValue.HalfLengthXnegY is not None


@pytest.mark.parametrize(
    "bounds, expected_size",
    [
        (acts.CylinderBounds(10.0, 20.0), 6),
        (acts.AnnulusBounds(10.0, 20.0, -0.2, 0.2), 7),
        (acts.RadialBounds(10.0, 20.0), 4),
        (acts.LineBounds(1.0, 100.0), 2),
        (acts.RectangleBounds(10.0, 20.0), 4),
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


def test_surface_bounds_get_method():
    cylinder = acts.CylinderBounds(10.0, 20.0)
    assert cylinder.get(acts.CylinderBoundsValue.R) == pytest.approx(10.0)
    assert cylinder.get(acts.CylinderBoundsValue.HalfLengthZ) == pytest.approx(20.0)
    assert cylinder.get(acts.CylinderBoundsValue.HalfPhiSector) == pytest.approx(
        cylinder.values()[2]
    )
    assert cylinder.get(acts.CylinderBoundsValue.AveragePhi) == pytest.approx(
        cylinder.values()[3]
    )
    assert cylinder.get(acts.CylinderBoundsValue.BevelMinZ) == pytest.approx(
        cylinder.values()[4]
    )
    assert cylinder.get(acts.CylinderBoundsValue.BevelMaxZ) == pytest.approx(
        cylinder.values()[5]
    )

    annulus = acts.AnnulusBounds(10.0, 20.0, -0.2, 0.2)
    assert annulus.get(acts.AnnulusBoundsValue.MinR) == pytest.approx(10.0)
    assert annulus.get(acts.AnnulusBoundsValue.MaxR) == pytest.approx(20.0)
    assert annulus.get(acts.AnnulusBoundsValue.MinPhiRel) == pytest.approx(-0.2)
    assert annulus.get(acts.AnnulusBoundsValue.MaxPhiRel) == pytest.approx(0.2)
    assert annulus.get(acts.AnnulusBoundsValue.OriginX) == pytest.approx(
        annulus.values()[4]
    )
    assert annulus.get(acts.AnnulusBoundsValue.OriginY) == pytest.approx(
        annulus.values()[5]
    )

    radial = acts.RadialBounds(11.0, 21.0)
    assert radial.get(acts.RadialBoundsValue.MinR) == pytest.approx(11.0)
    assert radial.get(acts.RadialBoundsValue.MaxR) == pytest.approx(21.0)
    assert radial.get(acts.RadialBoundsValue.HalfPhiSector) == pytest.approx(
        radial.values()[2]
    )
    assert radial.get(acts.RadialBoundsValue.AveragePhi) == pytest.approx(
        radial.values()[3]
    )

    line = acts.LineBounds(1.0, 100.0)
    assert line.get(acts.LineBoundsValue.R) == pytest.approx(1.0)
    assert line.get(acts.LineBoundsValue.HalfLengthZ) == pytest.approx(100.0)

    rectangle = acts.RectangleBounds(10.0, 20.0)
    assert rectangle.get(acts.RectangleBoundsValue.MinX) == pytest.approx(-10.0)
    assert rectangle.get(acts.RectangleBoundsValue.MinY) == pytest.approx(-20.0)
    assert rectangle.get(acts.RectangleBoundsValue.MaxX) == pytest.approx(10.0)
    assert rectangle.get(acts.RectangleBoundsValue.MaxY) == pytest.approx(20.0)

    trapezoid = acts.TrapezoidBounds(8.0, 12.0, 20.0, 0.1)
    assert trapezoid.get(acts.TrapezoidBoundsValue.HalfLengthXnegY) == pytest.approx(
        8.0
    )
    assert trapezoid.get(acts.TrapezoidBoundsValue.HalfLengthXposY) == pytest.approx(
        12.0
    )
    assert trapezoid.get(acts.TrapezoidBoundsValue.HalfLengthY) == pytest.approx(20.0)
    assert trapezoid.get(acts.TrapezoidBoundsValue.RotationAngle) == pytest.approx(0.1)


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
    assert plane.insideBounds(local, acts.BoundaryTolerance.none()) is True
    assert (
        plane.isOnSurface(
            gctx, global_position, direction, acts.BoundaryTolerance.none(), 0.0
        )
        is True
    )
    assert "PlaneSurface" in plane.name
    assert "PlaneSurface" in plane.toString(gctx)

    cylinder = acts.Surface.createCylinder(transform, acts.CylinderBounds(20.0, 50.0))
    assert isinstance(cylinder, acts.CylinderSurface)
    assert cylinder.type == acts.SurfaceType.Cylinder

    perigee = acts.Surface.createPerigee(acts.Vector3(0.0, 0.0, 0.0))
    assert isinstance(perigee, acts.PerigeeSurface)
    assert perigee.type == acts.SurfaceType.Perigee
