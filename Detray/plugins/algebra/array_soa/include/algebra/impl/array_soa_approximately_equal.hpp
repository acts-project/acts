// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "algebra/impl/array_soa_boolean.hpp"
#include "algebra/impl/array_soa_simd.hpp"
#include "detray/algebra/common/math.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <concepts>
#include <cstddef>
#include <limits>

namespace detray::algebra {

/// Elementwise compare two simd types according to a max relative error
/// tolerance
/// @see
/// https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
///
/// @note This is by no means safe for all comparisons. Use with caution!
///
/// @param a first simd type
/// @param b second simd type
/// @param rel_error maximal relative error
///
/// @returns true if the two simd types are elementwise approximately equal
template <concepts::value T, std::size_t W>
DETRAY_HOST_DEVICE constexpr bool approx_equal(
    const array_soa::simd<T, W> &a, const array_soa::simd<T, W> &b,
    const T rel_error = 16.f * std::numeric_limits<T>::epsilon(),
    const T max_error = std::numeric_limits<T>::epsilon()) {
  if constexpr (std::integral<T>) {
    return detray::algebra::boolean::all_of(a == b);
  } else {
    for (std::size_t i = 0u; i < W; ++i) {
      // Calculate the difference.
      const T diff{math::fabs(a[i] - b[i])};
      // If the numbers are close to zero
      if (diff <= max_error) {
        continue;
      }
      // Find the largest entries and scale the epsilon
      const T largest = math::max(math::fabs(a[i]), math::fabs(b[i]));
      if (!(diff <= largest * rel_error)) {
        return false;
      }
    }
    return true;
  }
}

}  // namespace detray::algebra
