// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"

#include <cassert>
#include <stdexcept>

namespace Acts {

/// Check orthogonality within the given tolerance.
/// @param rotation Matrix to check
/// @param tolerance Comparison tolerance (default: s_transformEquivalentTolerance)
inline bool isOrthogonal(const RotationMatrix3& rotation,
                         double tolerance = s_transformEquivalentTolerance) {
  return (rotation * rotation.transpose())
      .isApprox(RotationMatrix3::Identity(), tolerance);
}

/// Build a rigid transform: `p -> rotation * p + translation`.
/// @pre @p rotation is orthogonal (asserted).
/// @param rotation Local axes in the target frame
/// @param translation Local origin in the target frame
inline Transform3 makeTransform3(const RotationMatrix3& rotation,
                                 const Vector3& translation = Vector3::Zero()) {
  assert(isOrthogonal(rotation) &&
         "Transform3 requires an orthogonal rotation part");
  Transform3 transform = Transform3::Identity();
  transform.linear() = rotation;
  transform.translation() = translation;
  return transform;
}

/// Convert an affine transform, rejecting a non-orthogonal linear part.
/// @throws std::invalid_argument if the linear part is not orthogonal
inline Transform3 makeTransform3(const AffineTransform3& transform) {
  if (!isOrthogonal(transform.linear())) {
    throw std::invalid_argument(
        "Affine transform is not a rigid transformation");
  }
  return makeTransform3(transform.linear(), transform.translation());
}

}  // namespace Acts
