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

/// Correct all parameter values in a vector to be within their nominal range.
///
/// @tparam indices_t Parameter indices enum type
/// @tparam kParameters Parameter indices pack that defines the used parameter
template <typename indices_t, indices_t... kParameters>
struct ValueCorrector {
  using StorageSequence = std::make_index_sequence<sizeof...(kParameters)>;

  /// Run the value correction.
  ///
  /// @param[in,out] values Values vector which will be corrected in-place.
  /// @post All values in `values` are within their nominal range.
  template <typename values_t>
  static void run(Eigen::MatrixBase<values_t>& values) {
    runImpl(values, StorageSequence{});
  }

  template <typename values_t, std::size_t... kStorage>
  static void runImpl(Eigen::MatrixBase<values_t>& values,
                      std::index_sequence<kStorage...>) {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(values_t, sizeof...(kParameters));
    static_assert(sizeof...(kParameters) == sizeof...(kStorage),
                  "Parameters and storage index packs must have the same size");

    ((values[kStorage] =
          ParameterTraits<indices_t, kParameters>::getValue(values[kStorage])),
     ...);
  }
};

}  // namespace detail
}  // namespace Acts
