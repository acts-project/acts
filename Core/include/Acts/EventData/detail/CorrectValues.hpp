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

template <typename values_t, typename indices_t, indices_t... kParameters,
          std::size_t... kStorage>
constexpr void correctValuesImpl(
    Eigen::MatrixBase<values_t>& values,
    std::integer_sequence<indices_t, kParameters...>,
    std::index_sequence<kStorage...>) {
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(values_t, sizeof...(kParameters));
  static_assert(
      sizeof...(kParameters) == sizeof...(kStorage),
      "Parameters and storage index sequences must have the same size");

  ((values[kStorage] =
        ParameterTraits<indices_t, kParameters>::getValue(values[kStorage])),
   ...);
}

/// Correct all parameter values to be within their nominal range.
///
/// @tparam values_t Values vector type
/// @tparam indices_t Parameter indices enum type
/// @tparam kParameters Enum values pack; must have same size as values vector
///
/// @param[in,out] values Values vector which will be corrected in-place.
/// @param paramsSeq Index sequence that identifies the stored parameters
///
/// @post All values in `values` are within their nominal range.
template <typename values_t, typename indices_t, indices_t... kParameters>
constexpr void correctValues(
    Eigen::DenseBase<values_t>& values,
    std::integer_sequence<indices_t, kParameters...> paramsSeq) {
  using StorageSequence = std::make_index_sequence<sizeof...(kParameters)>;
  correctValuesImpl(values, paramsSeq, StorageSequence{});
}

}  // namespace detail
}  // namespace Acts
