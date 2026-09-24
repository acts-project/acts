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
#include "detray/algebra/common/boolean.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::algebra::boolean {

/// Boolean utilities on single values
/// @{
using detray::detail::all_of;
using detray::detail::any_of;
using detray::detail::count;
using detray::detail::none_of;
/// @}

/// Lane mask overloads
template <detray::concepts::array_soa_mask M>
DETRAY_HOST_DEVICE constexpr bool none_of(const M& mask) {
  bool ret{true};
  DETRAY_UNROLL_N(M::size())
  for (std::size_t i = 0u; i < M::size(); ++i) {
    ret = ret && !mask[i];
  }
  return ret;
}

template <detray::concepts::array_soa_mask M>
DETRAY_HOST_DEVICE constexpr bool any_of(const M& mask) {
  return !none_of(mask);
}

template <detray::concepts::array_soa_mask M>
DETRAY_HOST_DEVICE constexpr bool all_of(const M& mask) {
  bool ret{true};
  DETRAY_UNROLL_N(M::size())
  for (std::size_t i = 0u; i < M::size(); ++i) {
    ret = ret && mask[i];
  }
  return ret;
}

template <detray::concepts::array_soa_mask M>
DETRAY_HOST_DEVICE constexpr std::size_t count(const M& mask) {
  std::size_t n{0u};
  DETRAY_UNROLL_N(M::size())
  for (std::size_t i = 0u; i < M::size(); ++i) {
    n += mask[i] ? 1u : 0u;
  }
  return n;
}
/// @}

/// Array SoA overloads for masked assignments
/// @{
template <detray::concepts::array_soa_simd S,
          detray::concepts::array_soa_mask M, typename U>
  requires(std::is_scalar_v<U> &&
           std::convertible_to<U, typename S::value_type>)
constexpr void set_if(S& s, const M& mask, const U new_vlaue) {
  DETRAY_UNROLL_N(S::size())
  for (std::size_t i = 0u; i < S::size(); ++i) {
    s[i] = mask[i] ? static_cast<typename S::value_type>(new_vlaue) : s[i];
  }
}

template <detray::concepts::array_soa_simd S,
          detray::concepts::array_soa_mask M, typename U>
  requires(!std::is_scalar_v<U> && std::convertible_to<U, S>)
constexpr void set_if(S& s, const M& mask, const U& new_s) {
  DETRAY_UNROLL_N(S::size())
  for (std::size_t i = 0u; i < S::size(); ++i) {
    s[i] = mask[i] ? static_cast<typename S::value_type>(new_s[i]) : s[i];
  }
}

template <detray::concepts::array_soa_simd S,
          detray::concepts::array_soa_mask M>
constexpr void set_zero(S& s, const M& mask) {
  DETRAY_UNROLL_N(S::size())
  for (std::size_t i = 0u; i < S::size(); ++i) {
    s[i] = mask[i] ? static_cast<typename S::value_type>(0) : s[i];
  }
}

template <detray::concepts::array_soa_simd S,
          detray::concepts::array_soa_mask M>
constexpr void set_zero_inverted(S& s, const M& mask) {
  DETRAY_UNROLL_N(S::size())
  for (std::size_t i = 0u; i < S::size(); ++i) {
    s[i] = mask[i] ? s[i] : static_cast<typename S::value_type>(0);
  }
}
///@}

}  // namespace detray::algebra::boolean
