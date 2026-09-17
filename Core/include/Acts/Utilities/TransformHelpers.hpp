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

/// @brief Check whether a matrix is orthogonal, i.e. a rotation or reflection
///
/// This is the invariant @ref Acts::Transform3 relies on. Internal call sites
/// use @ref Acts::makeTransform3, which asserts it. Call this directly to
/// reject a matrix that comes from outside ACTS.
///
/// @param rotation The matrix to check
/// @return Whether the matrix is orthogonal within
///         @ref Acts::s_transformEquivalentTolerance
inline bool isOrthogonal(const RotationMatrix3& rotation) {
  return (rotation * rotation.transpose())
      .isApprox(RotationMatrix3::Identity(), s_transformEquivalentTolerance);
}

/// @brief Build a @ref Acts::Transform3 from a rotation and a translation
///
/// Maps `p -> rotation * p + translation`, i.e. the columns of @p rotation are
/// the local frame axes and @p translation its origin, both in the target
/// frame. Replaces the Eigen product `Translation3(translation) * rotation`,
/// which is affine and cannot be assigned to a @ref Acts::Transform3.
///
/// For the reverse order, `rotation * translation`, build the pure rotation
/// first and apply the translation to it: `makeTransform3(rotation) *
/// translation`.
///
/// Only a rotation that is a plain matrix needs this function. Eigen converts
/// its rotation types to a rigid transform directly, so write
/// `Transform3{AngleAxis3{angle, Vector3::UnitX()}}`,
/// `Transform3{Translation3{x, y, z}}` or
/// `Translation3{translation} * AngleAxis3{angle, axis}`. Products of
/// @ref Acts::Transform3 with @ref Acts::Translation3 or
/// @ref Acts::AngleAxis3 also stay rigid.
///
/// @param rotation The orthogonal linear part, i.e. the local frame axes
/// @param translation The local frame origin, in the target frame
/// @return The combined rigid transformation
inline Transform3 makeTransform3(const RotationMatrix3& rotation,
                                 const Vector3& translation = Vector3::Zero()) {
  assert(isOrthogonal(rotation) &&
         "Transform3 requires an orthogonal rotation part");
  Transform3 transform = Transform3::Identity();
  transform.linear() = rotation;
  transform.translation() = translation;
  return transform;
}

/// @brief Build a @ref Acts::Transform3 from a general affine transform
///
/// Use this for transforms that come from outside ACTS, e.g. an
/// @c Eigen::Affine3d from an external geometry source. Unlike the overload
/// above, this one checks the linear part in every build type.
///
/// @param transform The affine transform to convert
/// @throws std::invalid_argument if the linear part is not orthogonal
/// @return The equivalent rigid transformation
inline Transform3 makeTransform3(const AffineTransform3& transform) {
  if (!isOrthogonal(transform.linear())) {
    throw std::invalid_argument(
        "Affine transform is not a rigid transformation");
  }
  return makeTransform3(transform.linear(), transform.translation());
}

}  // namespace Acts
