// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/concepts.hpp"

namespace detray {

namespace algebra::constants {

/// Utilities to generate values (single value)
/// @{
template <detray::concepts::value T>
consteval T zero() {
  return 0.f;
}

template <detray::concepts::value T>
consteval T one() {
  return 1.f;
}

template <detray::concepts::value T>
consteval T iota() {
  return zero<T>();
}

template <detray::concepts::value T>
consteval std::size_t size() {
  return 1u;
}

template <detray::concepts::value T>
constexpr std::size_t size(T /*unused*/) {
  return 1u;
}
/// @}

}  // namespace algebra::constants

namespace detail {

// Pull the constants directly into the detray namespace
using namespace ::detray::algebra::constants;

}  // namespace detail

}  // namespace detray
