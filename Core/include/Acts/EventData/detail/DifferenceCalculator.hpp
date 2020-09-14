// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/detail/ParameterTraits.hpp"

#include <utility>

#include <Eigen/Core>

namespace Acts {
namespace detail {

/// Compute the minimal difference between parameters within the nominal range.
///
/// @tparam indices_t Parameter indices enum type
/// @tparam kParameters Parameter indices pack that defines the used parameter
template <typename indices_t, indices_t... kParameters>
struct DifferenceCalculator {
  using StorageSequence = std::make_index_sequence<sizeof...(kParameters)>;

  /// Run the difference calculation.
  ///
  /// @param rhs Right-hand-side input vector
  /// @param lhs Left-hand-side input vector
  /// @return minimal, corrected (lhs - rhs) difference within the nominal range
  template <typename lhs_t, typename rhs_t>
  static auto run(const Eigen::MatrixBase<lhs_t>& lhs,
                  const Eigen::MatrixBase<rhs_t>& rhs) {
    return runImpl(lhs, rhs, StorageSequence{});
  }

  template <typename lhs_t, typename rhs_t, std::size_t... kStorage>
  static auto runImpl(const Eigen::MatrixBase<lhs_t>& lhs,
                      const Eigen::MatrixBase<rhs_t>& rhs,
                      std::index_sequence<kStorage...>)
      -> std::decay_t<decltype((lhs - rhs).eval())> {
    EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(lhs_t, rhs_t);
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(lhs_t, sizeof...(kParameters));
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(rhs_t, sizeof...(kParameters));
    static_assert(sizeof...(kParameters) == sizeof...(kStorage),
                  "Parameters and storage index packs must have the same size");

    std::decay_t<decltype((lhs - rhs).eval())> residuals;
    ((residuals[kStorage] =
          ParameterTraits<indices_t, kParameters>::getDifference(
              lhs[kStorage], rhs[kStorage])),
     ...);
    return residuals;
  }
};

}  // namespace detail
}  // namespace Acts
