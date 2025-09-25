// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <tuple>
#include <utility>

namespace Acts {

/// Helper utility to allow indexed enumeration with structured binding
///
/// Usage:
///
/// for (auto [ i, value ] = enumerate(container) ) { ... };
///
/// with 'container' any stl-like container
/// @param iterable Container to enumerate
/// @return Enumerable wrapper with index and value pairs
template <
    typename container_type,
    typename index_type = typename std::decay_t<container_type>::size_type,
    typename container_type_iter =
        decltype(std::begin(std::declval<container_type>())),
    typename = decltype(std::end(std::declval<container_type>()))>
constexpr auto enumerate(container_type &&iterable) {
  struct iterator {
    index_type i;
    container_type_iter iter;

    bool operator!=(const iterator &rhs) const { return iter != rhs.iter; }

    /** Increase index and iterator at once */
    void operator++() {
      ++i;
      ++iter;
    }

    /** Tie them together for returning */
    auto operator*() const { return std::tie(i, *iter); }
  };
  struct iterable_wrapper {
    container_type iterable;
    auto begin() { return iterator{0, std::begin(iterable)}; }
    auto end() { return iterator{0, std::end(iterable)}; }
  };
  return iterable_wrapper{std::forward<container_type>(iterable)};
}

}  // namespace Acts
