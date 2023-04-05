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
    const RotationMatrix3& compositeRotation,
    const RotationMatrix3& relRotation) {
  // Get Euler angles for rotation represented by rotZ * rotY * rotX, i.e.
  // first rotation around x axis, then y axis, last z axis
  // The elements stored in rotAngles is (rotZ, rotY, rotX)
  // const Vector3 rotAngles = rotation.eulerAngles(2, 1, 0);
  // double sx = std::sin(rotAngles(2));
  // double cx = std::cos(rotAngles(2));
  // double sy = std::sin(rotAngles(1));
  // double cy = std::cos(rotAngles(1));
  // double sz = std::sin(rotAngles(0));
  // double cz = std::cos(rotAngles(0));
  // rotZ * rotY * rotX =
  // [ cz*cy  cz*sy*sx-cx*sz  sz*sx+cz*cx*sy ]
  // [ cy*sz  cz*cx+sz*sy*sx  cx*sz*sy-cz*sx ]
  // [ -sy   cy*sx         cy*cx        ]

  // Derivative of local x axis w.r.t. (rotX, rotY, rotZ)
  RotationMatrix3 rotToCompositeLocalXAxis = RotationMatrix3::Zero();
  rotToCompositeLocalXAxis.col(0) = compositeRotation * Vector3(0, 0, 0);
  rotToCompositeLocalXAxis.col(1) = compositeRotation * Vector3(0, 0, -1);
  rotToCompositeLocalXAxis.col(2) = compositeRotation * Vector3(0, 1, 0);
  // Derivative of local y axis w.r.t. (rotX, rotY, rotZ)
  RotationMatrix3 rotToCompositeLocalYAxis = RotationMatrix3::Zero();
  rotToCompositeLocalYAxis.col(0) = compositeRotation * Vector3(0, 0, 1);
  rotToCompositeLocalYAxis.col(1) = compositeRotation * Vector3(0, 0, 0);
  rotToCompositeLocalYAxis.col(2) = compositeRotation * Vector3(-1, 0, 0);
  // Derivative of local z axis w.r.t. (rotX, rotY, rotZ)
  RotationMatrix3 rotToCompositeLocalZAxis = RotationMatrix3::Zero();
  rotToCompositeLocalZAxis.col(0) = compositeRotation * Vector3(0, -1, 0);
  rotToCompositeLocalZAxis.col(1) = compositeRotation * Vector3(1, 0, 0);
  rotToCompositeLocalZAxis.col(2) = compositeRotation * Vector3(0, 0, 0);

  RotationMatrix3 rotToLocalXAxis = RotationMatrix3::Zero();
  RotationMatrix3 rotToLocalYAxis = RotationMatrix3::Zero();
  RotationMatrix3 rotToLocalZAxis = RotationMatrix3::Zero();
  rotToLocalXAxis = rotToCompositeLocalXAxis * relRotation(0, 0) +
                    rotToCompositeLocalYAxis * relRotation(1, 0) +
                    rotToCompositeLocalZAxis * relRotation(2, 0);
  rotToLocalYAxis = rotToCompositeLocalXAxis * relRotation(0, 1) +
                    rotToCompositeLocalYAxis * relRotation(1, 1) +
                    rotToCompositeLocalZAxis * relRotation(2, 1);
  rotToLocalZAxis = rotToCompositeLocalXAxis * relRotation(0, 2) +
                    rotToCompositeLocalYAxis * relRotation(1, 2) +
                    rotToCompositeLocalZAxis * relRotation(2, 2);

  return std::make_tuple(std::move(rotToLocalXAxis), std::move(rotToLocalYAxis),
                         std::move(rotToLocalZAxis));
}
