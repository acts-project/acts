// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <Eigen/Core>

namespace Acts::detail {

/// Check that all coefficients of a matrix are finite.
///
/// Same result as Eigen's `allFinite()`, which evaluates `x - x` coefficient
/// by coefficient with an early exit. For a finite coefficient `x - x` is
/// exactly zero, while it is NaN for an infinite or NaN coefficient, so the
/// sum of all differences is zero if and only if every coefficient is finite.
/// The sum is a vectorized, branch-free reduction, which is considerably
/// faster for the fixed-size track parameter matrices.
///
/// @tparam Derived Eigen matrix, map or expression type
/// @param matrix Matrix to check
/// @return Whether all coefficients are finite
template <typename Derived>
bool allFinite(const Eigen::MatrixBase<Derived>& matrix) {
  return (matrix.derived() - matrix.derived()).sum() == 0;
}

}  // namespace Acts::detail
