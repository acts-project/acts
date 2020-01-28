// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <type_traits>

#include "Acts/Utilities/detail/MPL/type_collector.hpp"

namespace ActsFatras {

/// The following operator has to be inplemented in order to satisfy
/// as an sampler for fast simulation
///
/// @code
///  bool
///  operator()(generator_t& generator,
///             const detector_t& detector,
///             const particle_t& in,
///             std::vector<particle_t>& out) const { return false; }
///
/// @endcode
namespace detail {
namespace {

template <typename T, typename generator_t, typename detector_t,
          typename particle_t,
          typename = decltype(std::declval<T>().operator()(
              std::declval<generator_t &>(), std::declval<const detector_t &>(),
              std::declval<particle_t &>(),
              std::declval<std::vector<particle_t> &>()))>

std::true_type test_physics_list(int);

template <typename, typename, typename, typename>
std::false_type test_physics_list(...);

template <typename T, typename generator_t, typename detector_t,
          typename particle_t>
struct process_signature_check
    : decltype(test_physics_list<T, generator_t, detector_t, particle_t>(0)) {};

}  // namespace

template <typename T, typename generator_t, typename detector_t,
          typename particle_t>
constexpr bool process_signature_check_v =
    process_signature_check<T, generator_t, detector_t, particle_t>::value;

}  // namespace detail
}  // namespace ActsFatras
