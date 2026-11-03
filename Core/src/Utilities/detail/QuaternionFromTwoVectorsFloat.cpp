// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// GCC (>=12) emits a spurious -Wmaybe-uninitialized from inside Eigen's SVD
// code used by setFromTwoVectors; suppress it for this translation unit.
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 12
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#include <Eigen/Dense>

namespace Acts::detail {

Eigen::Quaternion<float> quaternionFromTwoVectors(
    const Eigen::Matrix<float, 3, 1>& a, const Eigen::Matrix<float, 3, 1>& b) {
  return Eigen::Quaternion<float>().setFromTwoVectors(a, b);
}

}  // namespace Acts::detail

#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 12
#pragma GCC diagnostic pop
#endif
