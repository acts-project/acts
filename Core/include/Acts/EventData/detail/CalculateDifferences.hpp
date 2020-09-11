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

template <typename lhs_t, typename rhs_t, typename indices_t,
          indices_t... kParameters, std::size_t... kStorage>
constexpr auto calculateDifferencesImpl(
    const Eigen::MatrixBase<lhs_t>& lhs, const Eigen::MatrixBase<rhs_t>& rhs,
    std::integer_sequence<indices_t, kParameters...>,
    std::index_sequence<kStorage...>)
    -> std::decay_t<decltype((lhs - rhs).eval())> {
  using ResidualsVector = std::decay_t<decltype((lhs - rhs).eval())>;

  EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(lhs_t, rhs_t);
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(lhs_t, sizeof...(kParameters));
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(rhs_t, sizeof...(kParameters));
  static_assert(
      sizeof...(kParameters) == sizeof...(kStorage),
      "Parameters and storage index sequences must have the same size");

  ResidualsVector residuals;
  ((residuals[kStorage] =
        ParameterTraits<indices_t, kParameters>::getDifference(lhs[kStorage],
                                                               rhs[kStorage])),
   ...);
  return residuals;
}

/// Compute the minimal difference between parameters within the nominal range.
///
/// @tparam lhs_t Right-hand-side vector type
/// @tparam rhs_t Left-hand-side vector type
/// @tparam indices_t Parameter indices enum type
/// @tparam kIndices Enum values pack; must have same size as the input vectors
///
/// @param rhs Right-hand-side input vector
/// @param lhs Left-hand-side input vector
/// @param paramsSeq Index sequence that identifies the stored parameters
/// @return (lhs - rhs) but corrected to be minimal within the nominal range
template <typename lhs_t, typename rhs_t, typename indices_t,
          indices_t... kParameters>
constexpr auto calculateDifferences(
    const Eigen::MatrixBase<lhs_t>& lhs, const Eigen::MatrixBase<rhs_t>& rhs,
    std::integer_sequence<indices_t, kParameters...> paramsSeq) {
  using StorageSequence = std::make_index_sequence<sizeof...(kParameters)>;
  return calculateDifferencesImpl(lhs, rhs, paramsSeq, StorageSequence{});
}

}  // namespace detail
}  // namespace Acts
