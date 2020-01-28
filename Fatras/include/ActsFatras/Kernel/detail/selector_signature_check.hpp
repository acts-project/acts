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
/// as a process for fast simulation. The selector can access both,
/// current particle information, but also current detector information,
/// e.g. for deciding if an interaction or process has to take place
///
/// @code
///  bool
///  operator()(const detector_t& detector,
///             const particle_t& particle) const { return true; }
///
/// @endcode
namespace detail {
namespace {

template <typename T, typename detector_t, typename particle_t,
          typename = decltype(
              std::declval<T>().operator()(std::declval<const detector_t &>(),
                                           std::declval<const particle_t &>()))>

std::true_type test_selector_list(int);

template <typename, typename, typename>
std::false_type test_selector_list(...);

template <typename T, typename detector_t, typename particle_t>
struct selector_list_signature_check
    : decltype(test_selector_list<T, detector_t, particle_t>(0)) {};

}  // namespace

template <typename T, typename detector_t, typename particle_t>
constexpr bool selector_list_signature_check_v =
    selector_list_signature_check<T, detector_t, particle_t>::value;

}  // namespace detail
}  // namespace ActsFatras
