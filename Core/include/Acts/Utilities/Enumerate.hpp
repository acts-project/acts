// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <iterator>
#include <ranges>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Acts {

namespace detail {

/// Deduce the index type from the container's `size()` (via `std::size`),
/// falling back to `std::size_t` for ranges without a size (e.g. input
/// ranges).
template <typename T, typename = void>
struct EnumerateIndexType {
  using type = std::size_t;
};
template <typename T>
struct EnumerateIndexType<T,
                          std::void_t<decltype(std::size(std::declval<T>()))>> {
  using type = decltype(std::size(std::declval<T>()));
};

/// Iterator pairing an index with the element of an underlying iterator.
///
/// This is a LegacyInputIterator: the element reference is a proxy (a tuple),
/// which is why it cannot model the C++20 `std::forward_iterator` concept (in
/// C++20 `std::tuple` lacks the `common_reference` support that C++23 adds).
template <typename index_type, typename iterator_t>
struct EnumerateIterator {
  using iterator_category = std::input_iterator_tag;
  using value_type = std::tuple<index_type, std::iter_value_t<iterator_t>>;
  using reference = std::tuple<index_type, std::iter_reference_t<iterator_t>>;
  using pointer = void;
  using difference_type = std::iter_difference_t<iterator_t>;

  index_type i;
  iterator_t iter;

  bool operator==(const EnumerateIterator &rhs) const {
    return iter == rhs.iter;
  }

  /// Pre-increment: advance index and underlying iterator together.
  EnumerateIterator &operator++() {
    ++i;
    ++iter;
    return *this;
  }

  /// Post-increment.
  EnumerateIterator operator++(int) {
    EnumerateIterator tmp = *this;
    ++*this;
    return tmp;
  }

  reference operator*() const { return reference{i, *iter}; }
};

/// Range wrapper that yields `EnumerateIterator`s over a (possibly owned)
/// container.
template <std::ranges::range container_type, typename index_type,
          typename iterator_t>
struct EnumerateWrapper {
  container_type iterable;
  auto begin() {
    return EnumerateIterator<index_type, iterator_t>{0, std::begin(iterable)};
  }
  auto end() {
    return EnumerateIterator<index_type, iterator_t>{0, std::end(iterable)};
  }
};

}  // namespace detail

/// Helper utility to allow indexed enumeration with structured binding
///
/// Usage:
///
/// for (auto [ i, value ] : enumerate(container) ) { ... };
///
/// with 'container' any stl-like container
/// @param iterable Container to enumerate
/// @return Enumerable wrapper with index and value pairs
template <std::ranges::range container_type,
          typename index_type =
              typename detail::EnumerateIndexType<container_type>::type,
          typename container_type_iter =
              decltype(std::begin(std::declval<container_type>()))>
constexpr auto enumerate(container_type &&iterable) {
  return detail::EnumerateWrapper<container_type, index_type,
                                  container_type_iter>{
      std::forward<container_type>(iterable)};
}

}  // namespace Acts
