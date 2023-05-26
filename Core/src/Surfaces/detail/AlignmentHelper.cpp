// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/detail/AlignmentHelper.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

Acts::detail::RotationToAxes Acts::detail::rotationToLocalAxesDerivative(
    const RotationMatrix3& rotation) {
  // Get Euler angles for rotation representated by rotZ * rotY * rotX, i.e.
  // first rotation around x axis, then y axis, last z axis
  // The elements stored in rotAngles is (rotZ, rotY, rotX)
  const Vector3 rotAngles = rotation.eulerAngles(2, 1, 0);
  double sx = std::sin(rotAngles(2));
  double cx = std::cos(rotAngles(2));
  double sy = std::sin(rotAngles(1));
  double cy = std::cos(rotAngles(1));
  double sz = std::sin(rotAngles(0));
  double cz = std::cos(rotAngles(0));
  // rotZ * rotY * rotX =
  // [ cz*cy  cz*sy*sx-cx*sz  sz*sx+cz*cx*sy ]
  // [ cy*sz  cz*cx+sz*sy*sx  cx*sz*sy-cz*sx ]
  // [ -sy   cy*sx         cy*cx        ]

  // Derivative of local x axis w.r.t. (rotX, rotY, rotZ)
  RotationMatrix3 rotToLocalXAxis = RotationMatrix3::Zero();
  rotToLocalXAxis.col(0) = Vector3(0, 0, 0);
  rotToLocalXAxis.col(1) = Vector3(-cz * sy, -sz * sy, -cy);
  rotToLocalXAxis.col(2) = Vector3(-sz * cy, cz * cy, 0);
  // Derivative of local y axis w.r.t. (rotX, rotY, rotZ)
  RotationMatrix3 rotToLocalYAxis = RotationMatrix3::Zero();
  rotToLocalYAxis.col(0) =
      Vector3(cz * sy * cx + sz * sx, sz * sy * cx - cz * sx, cy * cx);
  rotToLocalYAxis.col(1) = Vector3(cz * cy * sx, sz * cy * sx, -sy * sx);
  rotToLocalYAxis.col(2) =
      Vector3(-sz * sy * sx - cz * cx, cz * sy * sx - sz * cx, 0);
  // Derivative of local z axis w.r.t. (rotX, rotY, rotZ)
  RotationMatrix3 rotToLocalZAxis = RotationMatrix3::Zero();
  rotToLocalZAxis.col(0) =
      Vector3(sz * cx - cz * sy * sx, -sz * sy * sx - cz * cx, -cy * sx);
  rotToLocalZAxis.col(1) = Vector3(cz * cy * cx, sz * cy * cx, -sy * cx);
  rotToLocalZAxis.col(2) =
      Vector3(cz * sx - sz * sy * cx, cz * sy * cx + sz * sx, 0);

  return std::make_tuple(std::move(rotToLocalXAxis), std::move(rotToLocalYAxis),
                         std::move(rotToLocalZAxis));
}
