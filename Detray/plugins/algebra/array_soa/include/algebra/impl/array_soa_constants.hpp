// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/array_soa_concepts.hpp"
#include "detray/algebra/common/constants.hpp"

namespace detray::algebra::constants {

/// Pulling in global constants for overload resolution
/// @{
using algebra::constants::iota;
using algebra::constants::one;
using algebra::constants::size;
using algebra::constants::zero;
/// @}

/// Utilities to generate values
/// @{
template <detray::concepts::array_soa_simd scalar_t>
consteval scalar_t zero() {
  return {0.f};
}

template <detray::concepts::array_soa_simd scalar_t>
consteval scalar_t one() {
  return {1.f};
}

template <detray::concepts::array_soa_simd scalar_t>
consteval scalar_t iota() {
  scalar_t ret;
  DETRAY_UNROLL_N(scalar_t::size())
  for (std::size_t i = 0u; i < scalar_t::size(); ++i) {
    ret[i] = static_cast<typename scalar_t::value_type>(i);
  }
  return ret;
}

template <detray::concepts::array_soa_simd scalar_t>
consteval std::size_t size() {
  return scalar_t::size();
}

template <detray::concepts::array_soa_simd scalar_t>
constexpr std::size_t size(const scalar_t& /*unused*/) {
  return scalar_t::size();
}
/// @}

}  // namespace detray::algebra::constants
