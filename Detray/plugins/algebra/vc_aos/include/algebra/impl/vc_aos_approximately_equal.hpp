// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

// System include(s)
#include <concepts>
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
template <typename simd1_t, typename simd2_t>
  requires((Vc::is_simd_vector<simd1_t>::value &&
            std::convertible_to<simd2_t, simd1_t>) ||
           (Vc::is_simd_vector<simd2_t>::value &&
            std::convertible_to<simd1_t, simd2_t>))
DETRAY_HOST_DEVICE constexpr auto approx_equal(
    const simd1_t a, const simd2_t b,
    const typename simd1_t::value_type rel_error =
        16.f * std::numeric_limits<typename simd1_t::value_type>::epsilon(),
    const typename simd1_t::value_type max_error =
        std::numeric_limits<typename simd1_t::value_type>::epsilon()) {
  static_assert(
      std::same_as<typename simd1_t::value_type, typename simd2_t::value_type>);

  if constexpr (std::integral<typename simd1_t::value_type>) {
    return (a == b).isFull();
  } else {
    // Calculate the difference.
    const simd1_t diff{Vc::abs(a - b)};

    // If the numbers are is close to zero
    if ((diff <= simd1_t(max_error)).isFull()) {
      return true;
    }

    // Find the largest entries and scale the epsilon
    const simd1_t largest = Vc::max(Vc::abs(a), Vc::abs(b));

    return (diff <= (largest * rel_error)).isFull();
  }
}

}  // namespace detray::algebra
