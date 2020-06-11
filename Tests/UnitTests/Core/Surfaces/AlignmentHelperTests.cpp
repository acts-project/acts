// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>
#include <cmath>

#include "Acts/Surfaces/detail/AlignmentHelper.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

/// Test for rotation matrix and calculation of derivative of rotated x/y/z axis
/// w.r.t. rotation parameters
BOOST_AUTO_TEST_CASE(alignment_helper_test) {
  // (a) Test with non-identity rotation matrix
  // Rotation angle parameters
  const double alpha = M_PI;
  const double beta = 0;
  const double gamma = M_PI / 2;
  // rotation around x axis
  AngleAxis3D rotX(alpha, Vector3D(1., 0., 0.));
  // rotation around y axis
  AngleAxis3D rotY(beta, Vector3D(0., 1., 0.));
  // rotation around z axis
  AngleAxis3D rotZ(gamma, Vector3D(0., 0., 1.));
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
  RotationMatrix3D expRot = RotationMatrix3D::Zero();
  expRot.col(0) = Vector3D(cz * cy, cy * sz, -sy);
  expRot.col(1) =
      Vector3D(cz * sy * sx - cx * sz, cz * cx + sz * sy * sx, cy * sx);
  expRot.col(2) =
      Vector3D(sz * sx + cz * cx * sy, cx * sz * sy - cz * sx, cy * cx);

  // Calculate the expected derivative of local x axis to its rotation
  RotationMatrix3D expRotToXAxis = RotationMatrix3D::Zero();
  expRotToXAxis.col(0) = Vector3D(0, 0, 0);
  expRotToXAxis.col(1) = Vector3D(-cz * sy, -sz * sy, -cy);
  expRotToXAxis.col(2) = Vector3D(-sz * cy, cz * cy, 0);

  // Calculate the expected derivative of local y axis to its rotation
  RotationMatrix3D expRotToYAxis = RotationMatrix3D::Zero();
  expRotToYAxis.col(0) =
      Vector3D(cz * sy * cx + sz * sx, sz * sy * cx - cz * sx, cy * cx);
  expRotToYAxis.col(1) = Vector3D(cz * cy * sx, sz * cy * sx, -sy * sx);
  expRotToYAxis.col(2) =
      Vector3D(-sz * sy * sx - cz * cx, cz * sy * sx - sz * cx, 0);

  // Calculate the expected derivative of local z axis to its rotation
  RotationMatrix3D expRotToZAxis = RotationMatrix3D::Zero();
  expRotToZAxis.col(0) =
      Vector3D(sz * cx - cz * sy * sx, -sz * sy * sx - cz * cx, -cy * sx);
  expRotToZAxis.col(1) = Vector3D(cz * cy * cx, sz * cy * cx, -sy * cx);
  expRotToZAxis.col(2) =
      Vector3D(cz * sx - sz * sy * cx, cz * sy * cx + sz * sx, 0);

  // Construct a transform
  Vector3D translation(0., 0., 0.);
  auto transform = std::make_shared<Transform3D>(Translation3D(translation));
  // Rotation with rotZ * rotY * rotX
  (*transform) *= rotZ;
  (*transform) *= rotY;
  (*transform) *= rotX;
  // Get the rotation of the transform
  const auto rotation = transform->rotation();

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
  RotationMatrix3D iRotation = RotationMatrix3D::Identity();

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
}  // namespace Test
}  // namespace Acts
