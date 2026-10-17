// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/detail/AlignmentHelper.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <utility>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(SurfacesSuite)

/// Test for rotation matrix and calculation of derivative of rotated x/y/z axis
/// w.r.t. rotation parameters
BOOST_AUTO_TEST_CASE(alignment_helper_test) {
  // (a) Test with non-identity rotation matrix
  // Rotation angle parameters
  const double alpha = std::numbers::pi;
  const double beta = 0.;
  const double gamma = std::numbers::pi / 2.;
  // rotation around x axis
  AngleAxis3 rotX(alpha, Vector3(1., 0., 0.));
  // rotation around y axis
  AngleAxis3 rotY(beta, Vector3(0., 1., 0.));
  // rotation around z axis
  AngleAxis3 rotZ(gamma, Vector3(0., 0., 1.));
  double sz = std::sin(gamma);
  double cz = std::cos(gamma);
  double sy = std::sin(beta);
  double cy = std::cos(beta);
  double sx = std::sin(alpha);
  double cx = std::cos(alpha);

  // Calculate the expected rotation matrix for rotZ * rotY * rotX,
  // (i.e. first rotation around x axis, then y axis, last z axis):
  // [ cz*cy  cz*sy*sx-cx*sz  sz*sx+cz*cx*sy ]
  // [ cy*sz  cz*cx+sz*sy*sx  cx*sz*sy-cz*sx ]
  // [ -sy    cy*sx           cy*cx          ]
  RotationMatrix3 expRot = RotationMatrix3::Zero();
  expRot.col(0) = Vector3(cz * cy, cy * sz, -sy);
  expRot.col(1) =
      Vector3(cz * sy * sx - cx * sz, cz * cx + sz * sy * sx, cy * sx);
  expRot.col(2) =
      Vector3(sz * sx + cz * cx * sy, cx * sz * sy - cz * sx, cy * cx);

  // Calculate the expected derivative of local x axis to its rotation
  RotationMatrix3 expRotToXAxis = RotationMatrix3::Zero();
  expRotToXAxis.col(0) = expRot * Vector3(0, 0, 0);
  expRotToXAxis.col(1) = expRot * Vector3(0, 0, -1);
  expRotToXAxis.col(2) = expRot * Vector3(0, 1, 0);

  // Calculate the expected derivative of local y axis to its rotation
  RotationMatrix3 expRotToYAxis = RotationMatrix3::Zero();
  expRotToYAxis.col(0) = expRot * Vector3(0, 0, 1);
  expRotToYAxis.col(1) = expRot * Vector3(0, 0, 0);
  expRotToYAxis.col(2) = expRot * Vector3(-1, 0, 0);

  // Calculate the expected derivative of local z axis to its rotation
  RotationMatrix3 expRotToZAxis = RotationMatrix3::Zero();
  expRotToZAxis.col(0) = expRot * Vector3(0, -1, 0);
  expRotToZAxis.col(1) = expRot * Vector3(1, 0, 0);
  expRotToZAxis.col(2) = expRot * Vector3(0, 0, 0);

  // Construct a transform
  Translation3 translation(Vector3(0., 0., 0.));
  Transform3 transform(translation);
  // Rotation with rotZ * rotY * rotX
  transform *= rotZ;
  transform *= rotY;
  transform *= rotX;
  // Get the rotation of the transform
  const auto rotation = transform.rotation();

  // Check if the rotation matrix is as expected
  CHECK_CLOSE_ABS(rotation, expRot, 1e-15);

  // Call the alignment helper to calculate the derivative of local frame axes
  // w.r.t its rotation
  const auto& [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] =
      detail::rotationToLocalAxesDerivative(rotation);

  // Check if the derivative for local x axis is as expected
  CHECK_CLOSE_ABS(rotToLocalXAxis, expRotToXAxis, 1e-15);

  // Check if the derivative for local y axis is as expected
  CHECK_CLOSE_ABS(rotToLocalYAxis, expRotToYAxis, 1e-15);

  // Check if the derivative for local z axis is as expected
  CHECK_CLOSE_ABS(rotToLocalZAxis, expRotToZAxis, 1e-15);

  // (b) Test with identity rotation matrix
  RotationMatrix3 iRotation = RotationMatrix3::Identity();

  // Call the alignment helper to calculate the derivative of local frame axes
  // w.r.t its rotation
  const auto& [irotToLocalXAxis, irotToLocalYAxis, irotToLocalZAxis] =
      detail::rotationToLocalAxesDerivative(iRotation);

  // The expected derivatives
  expRotToXAxis << 0, 0, 0, 0, 0, 1, 0, -1, 0;
  expRotToYAxis << 0, 0, -1, 0, 0, 0, 1, 0, 0;
  expRotToZAxis << 0, 1, 0, -1, 0, 0, 0, 0, 0;

  // Check if the derivative for local x axis is as expected
  CHECK_CLOSE_ABS(irotToLocalXAxis, expRotToXAxis, 1e-15);

  // Check if the derivative for local y axis is as expected
  CHECK_CLOSE_ABS(irotToLocalYAxis, expRotToYAxis, 1e-15);

  // Check if the derivative for local z axis is as expected
  CHECK_CLOSE_ABS(irotToLocalZAxis, expRotToZAxis, 1e-15);
}

namespace {

/// Move an object by local-frame alignment parameters (dt, dw)
Transform3 moveInLocalFrame(const Transform3& transform,
                            const AlignmentVector& params) {
  Transform3 delta = Transform3::Identity();
  delta.linear() = (AngleAxis3(params[eAlignmentRotation2], Vector3::UnitZ()) *
                    AngleAxis3(params[eAlignmentRotation1], Vector3::UnitY()) *
                    AngleAxis3(params[eAlignmentRotation0], Vector3::UnitX()))
                       .toRotationMatrix();
  delta.translation() = params.head<3>();
  return transform * delta;
}

/// Local-frame alignment parameters (dt, dw) moving an object from
/// @p nominal to @p moved
AlignmentVector localFrameParameters(const Transform3& nominal,
                                     const Transform3& moved) {
  const RotationMatrix3 rotation = nominal.rotation();
  const AngleAxis3 deltaRotation(rotation.transpose() * moved.rotation());
  AlignmentVector params;
  params.head<3>() =
      rotation.transpose() * (moved.translation() - nominal.translation());
  params.tail<3>() = deltaRotation.angle() * deltaRotation.axis();
  return params;
}

/// A composite structure and one of its components, both tilted and away from
/// the global origin
struct CompositeFixture {
  Transform3 composite = Translation3(10., -20., 300.) *
                         AngleAxis3(0.3, Vector3(1., 2., 3.).normalized());
  Transform3 component = composite * Translation3(50., 30., -40.) *
                         AngleAxis3(-1.1, Vector3(-2., 1., 0.5).normalized());

  /// The component transform after moving the composite as a rigid body
  Transform3 movedComponent(const AlignmentVector& compositeParams) const {
    return moveInLocalFrame(composite, compositeParams) * composite.inverse() *
           component;
  }
};

}  // namespace

/// Check the composite Jacobian against the finite-difference rigid motion of
/// a component, in local-frame and ACTS parameters
BOOST_AUTO_TEST_CASE(composite_jacobian_rigid_motion) {
  CompositeFixture fix;
  const AlignmentMatrix jacobian =
      detail::compositeToComponentJacobian(fix.composite, fix.component);
  const AlignmentMatrix jacobianActs =
      detail::localFrameToAlignmentParametersJacobian(fix.component) * jacobian;

  const double step = 1e-5;
  for (std::size_t iPar = 0; iPar < eAlignmentSize; ++iPar) {
    AlignmentVector delta = AlignmentVector::Zero();
    delta[iPar] = step;
    const Transform3 plus = fix.movedComponent(delta);
    const Transform3 minus = fix.movedComponent(-delta);

    const AlignmentVector numLocal =
        (localFrameParameters(fix.component, plus) -
         localFrameParameters(fix.component, minus)) /
        (2. * step);
    CHECK_CLOSE_ABS(numLocal, jacobian.col(iPar), 1e-6);

    // ACTS parameters: global center shift, local rotations
    AlignmentVector numActs = numLocal;
    numActs.head<3>() =
        (plus.translation() - minus.translation()) / (2. * step);
    CHECK_CLOSE_ABS(numActs, jacobianActs.col(iPar), 1e-6);
  }
}

/// Check the chained derivative of the bound local position on a component
/// plane w.r.t. the composite parameters, against moving the composite and
/// intersecting a straight track with the moved plane
BOOST_AUTO_TEST_CASE(composite_jacobian_plane_derivative) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  CompositeFixture fix;
  auto bounds = std::make_shared<const RectangleBounds>(1000., 1000.);
  auto plane = Surface::makeShared<PlaneSurface>(fix.component, bounds);

  // straight track crossing the plane at an angle
  const Vector3 direction =
      (fix.component.rotation() * Vector3(0.3, -0.2, 1.)).normalized();
  const Vector3 position =
      plane->localToGlobal(gctx, Vector2(3., -7.), direction);

  FreeVector pathDerivative = FreeVector::Zero();
  pathDerivative.head<3>() = direction;
  const AlignmentToBoundMatrix alignToBound = plane->alignmentToBoundDerivative(
      gctx, position, direction, pathDerivative);
  const Matrix<2, eAlignmentSize> expected =
      alignToBound.topRows<2>() *
      detail::localFrameToAlignmentParametersJacobian(fix.component) *
      detail::compositeToComponentJacobian(fix.composite, fix.component);

  auto boundLocal = [&](const AlignmentVector& compositeParams) {
    const Transform3 moved = fix.movedComponent(compositeParams);
    auto movedPlane = Surface::makeShared<PlaneSurface>(moved, bounds);
    const Vector3 normal = moved.rotation().col(2);
    const double path =
        normal.dot(moved.translation() - position) / normal.dot(direction);
    return movedPlane
        ->globalToLocal(gctx, position + path * direction, direction)
        .value();
  };

  const double step = 1e-5;
  for (std::size_t iPar = 0; iPar < eAlignmentSize; ++iPar) {
    AlignmentVector delta = AlignmentVector::Zero();
    delta[iPar] = step;
    const Vector2 numerical =
        (boundLocal(delta) - boundLocal(-delta)) / (2. * step);
    CHECK_CLOSE_ABS(numerical, expected.col(iPar), 1e-6);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
