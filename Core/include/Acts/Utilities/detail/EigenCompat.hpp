// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Eigen deprecated `eulerAngles` in favor of `canonicalEulerAngles`

#include <Eigen/Core>

namespace Acts::detail::EigenCompat {

template <typename Index, typename Derived>
inline auto canonicalEulerAngles(const Eigen::MatrixBase<Derived>& matrix,
                                 Index a0, Index a1, Index a2)
    -> Eigen::Matrix<typename Eigen::MatrixBase<Derived>::Scalar, 3, 1> {
#if EIGEN_VERSION_AT_LEAST(5, 0, 0)
  return matrix.canonicalEulerAngles(a0, a1, a2);
#else
  return matrix.eulerAngles(a0, a1, a2);
#endif
}

#if EIGEN_VERSION_AT_LEAST(5, 0, 0)
constexpr auto all = Eigen::placeholders::all;
#else
static const auto all = Eigen::all;
#endif

}  // namespace Acts::detail::EigenCompat
