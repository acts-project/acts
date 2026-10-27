// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/find_bound.hpp"
#include "detray/utils/ranges/detail/iterator_functions.hpp"
#include "detray/utils/sort.hpp"

// System include(s)
#include <algorithm>

namespace detray {

/// Sort a range sequentially on the host or device.
/// @tparam rand_iter_t random-access iterator over the elements to sort
/// @tparam sentinel_t sentinel marking the end of the range
/// @param first iterator to the first element
/// @param last sentinel after the final element
template <std::random_access_iterator rand_iter_t,
          std::sentinel_for<rand_iter_t> sentinel_t>
DETRAY_HOST_DEVICE inline void sequential_sort(rand_iter_t first,
                                               sentinel_t last) {
#if defined(DETRAY_NO_DEVICE)
  std::ranges::sort(first, last);
#else
  detray::detail::selection_sort(first, last);
#endif
}

/// Find the first element satisfying a predicate.
/// @tparam rand_iter_t random-access iterator over the elements
/// @tparam sentinel_t sentinel marking the end of the range
/// @tparam Predicate callable returning true for a matching element
/// @param first iterator to the first element
/// @param last sentinel after the final element
/// @param comp predicate applied to each element
/// @return iterator to the first match, or last when no element matches
template <std::random_access_iterator rand_iter_t,
          std::sentinel_for<rand_iter_t> sentinel_t, class Predicate>
DETRAY_HOST_DEVICE inline auto find_if(rand_iter_t first, sentinel_t last,
                                       Predicate&& comp) {
  for (rand_iter_t i = first; i != last; ++i) {
    if (std::forward<Predicate>(comp)(*i)) {
      return i;
    }
  }

  return last;
}

/// Find the first element that is not less than a value in a sorted range.
/// @tparam forw_iter_t forward iterator over the sorted elements
/// @tparam sentinel_t sentinel marking the end of the range
/// @tparam Value type of the lookup value
/// @param first iterator to the first element
/// @param last sentinel after the final element
/// @param value value to locate
/// @return iterator to the first element not less than value
template <std::forward_iterator forw_iter_t,
          std::sentinel_for<forw_iter_t> sentinel_t, typename Value>
DETRAY_HOST_DEVICE inline auto lower_bound(forw_iter_t first, sentinel_t last,
                                           const Value& value) {
#if defined(DETRAY_NO_DEVICE)
  return std::ranges::lower_bound(first, last, value);
#else
  return detray::detail::lower_bound(first, last, value);
#endif
}

/// Find the first element greater than a value in a sorted range.
/// @tparam forw_iter_t forward iterator over the sorted elements
/// @tparam sentinel_t sentinel marking the end of the range
/// @tparam Value type of the lookup value
/// @param first iterator to the first element
/// @param last sentinel after the final element
/// @param value value to locate
/// @return iterator to the first element greater than value
template <std::forward_iterator forw_iter_t,
          std::sentinel_for<forw_iter_t> sentinel_t, typename Value>
DETRAY_HOST_DEVICE inline auto upper_bound(forw_iter_t first, sentinel_t last,
                                           const Value& value) {
#if defined(DETRAY_NO_DEVICE)
  return std::ranges::upper_bound(first, last, value);
#else
  return detray::detail::upper_bound(first, last, value);
#endif
}

}  // namespace detray
