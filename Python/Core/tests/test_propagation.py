import sys

import numpy as np
import pytest

import acts


def position_on_target(result, target, context):
    return target.localToGlobal(
        context,
        acts.Vector2(result.parameters[0], result.parameters[1]),
        acts.Vector3(0.0, 0.0, 1.0),
    )


def test_propagator_max_step_size():
    """The bound step limit controls actual propagation, not only a Python value."""
    geo_context = acts.GeometryContext.dangerouslyDefaultConstruct()
    field_context = acts.MagneticFieldContext()
    options = acts.PropagatorPlainOptions(geo_context, field_context)
    assert options.stepping.maxStepSize == sys.float_info.max

    propagator = acts.EigenVoidPropagator(
        acts.EigenStepper(acts.ConstantBField(acts.Vector3(0.0, 0.0, 0.0))),
        acts.VoidNavigator(),
    )
    start = acts.BoundTrackParameters.createCurvilinear(
        acts.Vector4(0.0, 0.0, 0.0, 0.0),
        acts.Vector3(0.0, 0.0, 1.0),
        1.0 / acts.UnitConstants.GeV,
        None,
        acts.ParticleHypothesis.pion,
    )
    distance = 20.0 * acts.UnitConstants.mm
    target = acts.Surface.createPlane(
        acts.Transform3(acts.Vector3(0.0, 0.0, distance)),
        acts.RectangleBounds(
            10.0 * acts.UnitConstants.mm, 10.0 * acts.UnitConstants.mm
        ),
    )
    options.maxSteps = 2
    result = propagator.propagateToSurface(start, target, options)
    np.testing.assert_allclose(
        position_on_target(result, target, geo_context), [0.0, 0.0, distance], atol=1e-6
    )

    options.stepping.maxStepSize = 1.0 * acts.UnitConstants.mm
    assert options.stepping.maxStepSize == 1.0 * acts.UnitConstants.mm
    with pytest.raises(RuntimeError, match="Propagation to surface failed"):
        propagator.propagateToSurface(start, target, options)

    options.maxSteps = 100
    result = propagator.propagateToSurface(start, target, options)
    np.testing.assert_allclose(
        position_on_target(result, target, geo_context), [0.0, 0.0, distance], atol=1e-6
    )


def test_limited_step_reaches_inner_barrel_plane():
    """A short step resolves a finite target close to a curved trajectory start."""
    geo_context = acts.GeometryContext.dangerouslyDefaultConstruct()
    field_context = acts.MagneticFieldContext()
    options = acts.PropagatorPlainOptions(geo_context, field_context)
    options.stepping.maxStepSize = 10.0 * acts.UnitConstants.mm
    field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 3.0 * acts.UnitConstants.T))
    propagator = acts.EigenVoidPropagator(
        acts.EigenStepper(field), acts.VoidNavigator()
    )
    eta = 1.4
    phi = np.deg2rad(292.5)
    start = acts.BoundTrackParameters.createCurvilinear(
        acts.Vector4(0.5, 0.5, 150.0, 0.0),
        phi,
        2.0 * np.arctan(np.exp(-eta)),
        1.0 / (np.cosh(eta) * acts.UnitConstants.GeV),
        None,
        acts.ParticleHypothesis.pion,
    )
    angle = np.deg2rad(280.0)
    rotation = acts.RotationMatrix3(
        acts.Vector3(-np.sin(angle), np.cos(angle), 0.0),
        acts.Vector3(0.0, 0.0, 1.0),
        acts.Vector3(np.cos(angle), np.sin(angle), 0.0),
    )
    target = acts.Surface.createPlane(
        acts.Transform3(
            acts.Vector3(34.0 * np.cos(angle), 34.0 * np.sin(angle), 208.8),
            rotation,
        ),
        acts.RectangleBounds(10.0, 9.6),
    )
    result = propagator.propagateToSurface(start, target, options)
    np.testing.assert_allclose(
        position_on_target(result, target, geo_context),
        [13.42655778, -32.15704041, 216.88629003],
        atol=1e-5,
    )
    assert abs(result.parameters[0]) < 10.0
    assert abs(result.parameters[1]) < 9.6
