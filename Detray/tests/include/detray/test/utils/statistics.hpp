// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// TODO: Remove this when gcc fixes their false positives.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic warning "-Wfree-nonheap-object"
#endif

// Project include(s)
#include "detray/utils/ranges.hpp"

// System include(s).
#include <algorithm>
#include <numeric>
#include <type_traits>

namespace detray::statistics {

/// @returns the sample mean over a range of numbers
template <detray::ranges::range range_t>
inline auto mean(const range_t& r) noexcept(false)
    -> detray::ranges::range_value_t<range_t> {
  using value_t = detray::ranges::range_value_t<range_t>;

  value_t sum = std::accumulate(r.begin(), r.end(), value_t{});
  value_t mean = sum * (1.0f / static_cast<value_t>(r.size()));

  return mean;
}

/// @returns RMS as the variance over a range of numbers with a true mean value
template <detray::ranges::range range_t>
inline auto rms(
    const range_t& r,
    const detray::ranges::range_value_t<range_t> mean) noexcept(false)
    -> detray::ranges::range_value_t<range_t> {
  using value_t = detray::ranges::range_value_t<range_t>;

  std::vector<value_t> diff(r.size());
  std::ranges::transform(r, diff.begin(),
                         [mean](value_t x) { return x - mean; });
  value_t sq_sum =
      std::inner_product(diff.begin(), diff.end(), diff.begin(), value_t{});
  value_t variance = sq_sum * (1.0f / static_cast<value_t>(r.size()));

  return variance;
}

/// @returns the sample variance over a range of numbers with a mean value
/// calculated from the range of numbers
template <detray::ranges::range range_t>
inline auto variance(const range_t& r) noexcept(false)
    -> detray::ranges::range_value_t<range_t> {
  using value_t = detray::ranges::range_value_t<range_t>;

  // Variance is zero for single element
  if (r.size() == 1u) {
    return 0.f;
  }

  value_t mean = statistics::mean(r);

  std::vector<value_t> diff(r.size());
  std::ranges::transform(r, diff.begin(),
                         [mean](value_t x) { return x - mean; });
  value_t sq_sum =
      std::inner_product(diff.begin(), diff.end(), diff.begin(), value_t{});

  // Divide with n-1 to avoid a bias
  value_t variance = sq_sum * (1.0f / static_cast<value_t>(r.size() - 1u));

  return variance;
}

}  // namespace detray::statistics
