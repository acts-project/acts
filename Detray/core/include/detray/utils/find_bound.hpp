// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/ranges.hpp"

// System include(s).
#include <iterator>

namespace detray::detail {

/// Implementation of @c upper_bound
/// @see
/// https://github.com/gcc-mirror/gcc/blob/8f87b3c5ecd47f6ac0d7407ae5d436a12fb169dd/libstdc%2B%2B-v3/include/bits/stl_algo.h
template <std::forward_iterator iterator_t, typename T>
DETRAY_HOST_DEVICE constexpr iterator_t upper_bound(iterator_t first,
                                                    iterator_t last,
                                                    const T& value) {
  using difference_t = std::iter_difference_t<iterator_t>;

  difference_t len{detray::ranges::distance(first, last)};

  // binary search
  while (len > 0) {
    difference_t half{len >> 1};
    iterator_t middle{first};
    detray::ranges::advance(middle, half);
    if (value < *middle) {
      len = half;
    } else {
      first = middle;
      ++first;
      len -= half + 1;
    }
  }
  return first;
}

/// Implementation of @c lower_bound
/// @see
/// https://github.com/gcc-mirror/gcc/blob/8f87b3c5ecd47f6ac0d7407ae5d436a12fb169dd/libstdc%2B%2B-v3/include/bits/stl_algobase.h
template <std::forward_iterator iterator_t, typename T>
DETRAY_HOST_DEVICE constexpr iterator_t lower_bound(iterator_t first,
                                                    iterator_t last,
                                                    const T& value) {
  using difference_t = std::iter_difference_t<iterator_t>;

  difference_t len{detray::ranges::distance(first, last)};

  // binary search
  while (len > 0) {
    difference_t half{len >> 1};
    iterator_t middle{first};
    detray::ranges::advance(middle, half);
    if (*middle < value) {
      first = middle;
      ++first;
      len -= half + 1;
    } else {
      len = half;
    }
  }
  return first;
}

}  // namespace detray::detail
