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

// System include(s)
#include <cassert>
#include <concepts>
#include <iterator>
#include <type_traits>

/// Reimplement some useful iterator functionality with __host__ __device__
/// quialifiers.
///
/// Adapted from libstdc++ e.g. @see
/// https://github.com/gcc-mirror/gcc/blob/16e2427f50c208dfe07d07f18009969502c25dc8/libstdc%2B%2B-v3/include/bits/stl_iterator_base_funcs.h
namespace detray::ranges {

// Reuse iterator tags for the host
using input_iterator_tag = std::input_iterator_tag;
using output_iterator_tag = std::output_iterator_tag;
using forward_iterator_tag = std::forward_iterator_tag;
using bidirectional_iterator_tag = std::bidirectional_iterator_tag;
using random_access_iterator_tag = std::random_access_iterator_tag;

// Iterator concepts
template <class I>
concept input_iterator =
    std::derived_from<typename std::iterator_traits<I>::iterator_category,
                      detray::ranges::input_iterator_tag> &&
    std::input_iterator<I>;

template <class I, typename T>
concept output_iterator =
    std::derived_from<typename std::iterator_traits<I>::iterator_category,
                      detray::ranges::input_iterator_tag> &&
    std::output_iterator<I, T>;

template <class I>
concept forward_iterator =
    std::derived_from<typename std::iterator_traits<I>::iterator_category,
                      detray::ranges::forward_iterator_tag> &&
    std::forward_iterator<I>;

template <class I>
concept bidirectional_iterator =
    std::derived_from<typename std::iterator_traits<I>::iterator_category,
                      detray::ranges::bidirectional_iterator_tag> &&
    std::bidirectional_iterator<I>;

template <class I>
concept random_access_iterator =
    std::derived_from<typename std::iterator_traits<I>::iterator_category,
                      detray::ranges::random_access_iterator_tag> &&
    std::random_access_iterator<I>;

namespace detail {
/// Simply import the std versions of basic iterator functionality where
/// possible
/// @{
using std::begin;
using std::cbegin;
using std::cend;
using std::crbegin;
using std::crend;
using std::end;
using std::rbegin;
using std::rend;

using std::size;
// TODO: ranges::ssize;
using std::data;
using std::empty;
// TODO: ranges::cdata
/// @}

///  @brief Reimplement std::distance.
/// @{
// used by every iterator up to and including bidirectional iterators
template <std::input_iterator iterator_t>
DETRAY_HOST_DEVICE constexpr std::iter_difference_t<iterator_t> distance_impl(
    iterator_t first, iterator_t last,
    detray::ranges::input_iterator_tag /*tag*/) {
  std::iter_difference_t<iterator_t> d{0};
  // simply count
  while (first != last) {
    ++first;
    ++d;
  }
  return d;
}

// random access iterators specialization
template <std::random_access_iterator iterator_t>
DETRAY_HOST_DEVICE constexpr std::iter_difference_t<iterator_t> distance_impl(
    iterator_t first, iterator_t last,
    detray::ranges::random_access_iterator_tag /*tag*/) {
  // use operator-
  return last - first;
}

template <std::input_iterator iterator_t>
DETRAY_HOST_DEVICE constexpr std::iter_difference_t<iterator_t> distance(
    iterator_t first, iterator_t last) {
  return distance_impl(
      first, last,
      typename std::iterator_traits<iterator_t>::iterator_category{});
}
/// @}

///  @brief Reimplement std::advance.
/// @{
template <std::input_iterator iterator_t, typename dist_t>
DETRAY_HOST_DEVICE constexpr void advance_impl(
    iterator_t& itr, dist_t d, detray::ranges::input_iterator_tag /*tag*/) {
  static_assert(std::is_integral_v<dist_t>);
  assert(d > 0);
  // simply count
  while (d--) {
    ++itr;
  }
}

// bidirectional iterators specialization
template <std::bidirectional_iterator iterator_t, typename dist_t>
DETRAY_HOST_DEVICE constexpr void advance_impl(
    iterator_t& itr, dist_t d,
    detray::ranges::bidirectional_iterator_tag /*tag*/) {
  static_assert(std::is_integral_v<dist_t>);
  if (d > 0) {
    while (d--) {
      ++itr;
    }
  } else {
    while (d++) {
      --itr;
    }
  }
}

// random access iterators specialization
template <std::random_access_iterator iterator_t, typename dist_t>
DETRAY_HOST_DEVICE constexpr void advance_impl(
    iterator_t& itr, dist_t d,
    detray::ranges::random_access_iterator_tag /*tag*/) {
  static_assert(std::is_integral_v<dist_t>);
  itr += d;
}

template <std::input_iterator iterator_t, typename dist_t>
DETRAY_HOST_DEVICE constexpr void advance(iterator_t& itr, dist_t d) {
  return advance_impl(
      itr, d, typename std::iterator_traits<iterator_t>::iterator_category{});
}
/// @}

///  @brief Reimplement std::next and std::prev.
/// @{
template <std::input_iterator iterator_t>
DETRAY_HOST_DEVICE constexpr iterator_t next(
    iterator_t itr, std::iter_difference_t<iterator_t> d = 1) {
  detray::ranges::detail::advance(itr, d);
  return itr;
}

template <std::bidirectional_iterator iterator_t>
DETRAY_HOST_DEVICE constexpr iterator_t prev(
    iterator_t itr, std::iter_difference_t<iterator_t> d = 1) {
  detray::ranges::detail::advance(itr, -d);
  return itr;
}
/// @}

}  // namespace detail

}  // namespace detray::ranges
