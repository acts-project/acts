// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/vc_soa_concepts.hpp"
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
/// @tparam T Vc vector type.
/// @returns zero in all lanes
template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
constexpr T zero() {
  return T::Zero();
}

/// @tparam T Vc vector type.
/// @returns one in all lanes
template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
constexpr T one() {
  return T::One();
}

/// @tparam T Vc vector type.
/// @returns the lane index in every lane
template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
constexpr T iota() {
  return T::IndexesFromZero();
}

/// @tparam T Vc vector type.
/// @returns the number of lanes
template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
constexpr std::size_t size() {
  return T::Size;
}

/// @tparam T Vc vector type.
/// @returns the number of lanes
template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
constexpr std::size_t size(const T& /*unused*/) {
  return T::size();
}
/// @}

}  // namespace detray::algebra::constants
