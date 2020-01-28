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

template <typename... processes>
struct physics_list_impl;

/// Recursive call pattern
/// - make sure it bails out if one process triggers a stop
template <typename first, typename... others>
struct physics_list_impl<first, others...> {
  template <typename T, typename generator_t, typename detector_t,
            typename particle_t>
  static bool process(const T &process_tuple, generator_t &gen,
                      const detector_t &det, particle_t &in,
                      std::vector<particle_t> &out) {
    // pick the first process
    const auto &this_process = std::get<first>(process_tuple);
    bool this_process_kills = this_process(gen, det, in, out);
    // recursive call on the remaining ones
    return (this_process_kills || physics_list_impl<others...>::process(
                                      process_tuple, gen, det, in, out));
  }
};

/// Final call pattern
template <typename last>
struct physics_list_impl<last> {
  template <typename T, typename generator_t, typename detector_t,
            typename particle_t>
  static bool process(const T &process_tuple, generator_t &gen,
                      const detector_t &det, particle_t &in,
                      std::vector<particle_t> &out) {
    // this is the last process in the tuple
    const auto &this_process = std::get<last>(process_tuple);
    return this_process(gen, det, in, out);
  }
};

/// Empty call pattern
template <>
struct physics_list_impl<> {
  template <typename T, typename generator_t, typename detector_t,
            typename particle_t>

  static bool process(const T &, generator_t &, const detector_t &,
                      const particle_t &, std::vector<particle_t> &) {
    return false;
  }
};

}  // namespace
}  // namespace detail
}  // namespace ActsFatras
