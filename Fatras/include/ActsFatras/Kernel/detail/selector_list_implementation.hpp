// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace ActsFatras {
namespace detail {
namespace {

template <typename... selectors>
struct selector_list_impl;

/// Recursive call pattern
/// - make sure it bails out
template <typename first, typename... others>
struct selector_list_impl<first, others...> {
  template <typename T, typename detector_t, typename particle_t>
  static bool select(const T &slector_tuple, const detector_t &detector,
                     const particle_t &particle, bool inclusive) {
    // pick the first select
    const auto &this_selector = std::get<first>(slector_tuple);
    bool selected = this_selector(detector, particle);
    // recursive call on the remaining ones, none of the selectors
    // is allowed to fail - a single failed selector rejects the particle
    // @todo check if rhs is actually evaluated if lhs fails (should not!)
    if (inclusive)
      return (selected || selector_list_impl<others...>::select(
                              slector_tuple, detector, particle, inclusive));
    return (selected && selector_list_impl<others...>::select(
                            slector_tuple, detector, particle, inclusive));
  }
};

/// Final call pattern
template <typename last>
struct selector_list_impl<last> {
  template <typename T, typename detector_t, typename particle_t>
  static bool select(const T &slector_tuple, const detector_t &detector,
                     const particle_t &particle, bool) {
    // this is the last select in the tuple
    const auto &this_selector = std::get<last>(slector_tuple);
    return this_selector(detector, particle);
  }
};

/// Empty call pattern
template <>
struct selector_list_impl<> {
  template <typename T, typename detector_t, typename particle_t>
  static bool select(const T &, const detector_t &, const particle_t &, bool) {
    return true;
  }
};

}  // namespace
}  // namespace detail
}  // namespace ActsFatras
