// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// System include(s)
#include <ratio>

namespace detray {

/// Helper struct to sum over variadic std::ratio
/// @{
template <typename... ratios>
struct ratio_sum_helper;

/// Base case for a one-element ratio sum.
template <typename R>
struct ratio_sum_helper<R> {
  /// The sole ratio, unchanged.
  using ratio = R;
};

/// Recursive case that adds two ratios before processing any remainder.
template <typename ratio1, typename ratio2, typename... ratios>
struct ratio_sum_helper<ratio1, ratio2, ratios...> {
  /// Sum of the first two ratios before folding in the remainder.
  using first_two_sum = std::ratio_add<ratio1, ratio2>;

  /// Recursive helper used to add any remaining ratios.
  using next_helper = ratio_sum_helper<first_two_sum, ratios...>;

  // recursive summation of first two std::ratio
  /// Result of adding all ratios in the parameter pack.
  using ratio =
      typename std::conditional_t<sizeof...(ratios) == 0, first_two_sum,
                                  typename next_helper::ratio>;
};
/// @}

/// Struct to sum over variadic std::ratio
template <typename... ratios>
struct ratio_sum {
  /// Helper that performs the recursive sum.
  using helper = ratio_sum_helper<ratios...>;
  /// The exact rational sum of all input ratios.
  using ratio = typename helper::ratio;
};

/// Helper trait to check if the ratio is one
template <typename R>
struct is_ratio_one {
  /// Whether the ratio's numerator and denominator are equal.
  static constexpr bool value = (R::num == R::den);
};

/// True when the supplied ratio is exactly one.
template <class R>
inline constexpr bool is_ratio_one_v = is_ratio_one<R>::value;

}  // namespace detray
