// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/find_bound.hpp"
#include "detray/utils/ranges/detail/iterator_functions.hpp"

// System include(s).
#include <algorithm>
#include <functional>

namespace detray::detail {

template <std::random_access_iterator RandomIt>
DETRAY_HOST_DEVICE inline void insertion_sort(RandomIt first, RandomIt last) {
  for (auto it = first; it != last; it++) {
    // Searching the upper bound, i.e., first
    // element greater than *it from beginning
    auto const insertion_point = detray::detail::upper_bound(first, it, *it);

    // Shifting the unsorted part
    std::ranges::rotate(insertion_point, it, it + 1);
  }
}

// Function to sort the array
template <template <typename...> class vector_t, typename TYPE>
DETRAY_HOST_DEVICE inline void insertion_sort(vector_t<TYPE> &vec) {
  detray::detail::insertion_sort(vec.begin(), vec.end());
}

template <std::random_access_iterator RandomIt, class Comp = std::less<void>>
DETRAY_HOST_DEVICE inline void selection_sort(RandomIt first, RandomIt last,
                                              Comp &&comp = Comp()) {
  for (RandomIt i = first; i < (last - 1); ++i) {
    RandomIt k = i;

    for (RandomIt j = i + 1; j < last; ++j) {
      if (comp(*j, *k)) {
        k = j;
      }
    }

    if (k != i) {
      auto t = *i;
      *i = *k;
      *k = t;
    }
  }
}

// Function to sort the array
template <template <typename...> class vector_t, typename TYPE>
DETRAY_HOST_DEVICE inline void selection_sort(vector_t<TYPE> &vec) {
  detray::detail::selection_sort(vec.begin(), vec.end());
}

}  // namespace detray::detail
