// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/detail/TransformComparator.hpp"

#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/detail/EigenCompat.hpp"

namespace Acts::detail {
TransformComparator::TransformComparator(const double transTolerance,
                                         const double rotTolerance)
    : m_tolTrans{transTolerance}, m_tolRot{rotTolerance} {}
int TransformComparator::compare(const Acts::RotationMatrix3& a,
                                 const Acts::RotationMatrix3& b) const {
  const Acts::Vector3 anglesA =
      detail::EigenCompat::canonicalEulerAngles(a, 2, 1, 0);
  const Acts::Vector3 anglesB =
      detail::EigenCompat::canonicalEulerAngles(b, 2, 1, 0);
  for (int i = 0; i < 3; ++i) {
    const double diff = anglesA[i] - anglesB[i];
    if (Acts::abs(diff) > m_tolRot) {
      return copySign(1, diff);
    }
  }
  return 0;
}
int TransformComparator::compare(const Acts::Transform3& a,
                                 const Acts::Transform3& b) const {
  if (const int tCmp = compare<3>(a.translation(), b.translation());
      tCmp != 0) {
    return tCmp;
  }
  return compare(a.rotation(), b.rotation());
}
bool TransformComparator::operator()(const Acts::Transform3& a,
                                     const Acts::Transform3& b) const {
  return compare(a, b) < 0;
}
bool TransformComparator::operator()(const Acts::RotationMatrix3& a,
                                     const RotationMatrix3& b) const {
  return compare(a, b) < 0;
}
bool TransformComparator::operator()(const Acts::Vector3& a,
                                     const Acts::Vector3& b) const {
  return compare<3>(a, b) < 0;
}
bool TransformComparator::operator()(const Acts::Vector2& a,
                                     const Acts::Vector2& b) const {
  return compare<2>(a, b) < 0;
}
}  // namespace Acts::detail
