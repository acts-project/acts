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
#include "algebra/impl/array_soa_simd.hpp"
#include "detray/algebra/common/boolean.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <bit>
#include <concepts>
#include <cstddef>
#include <type_traits>

namespace detray::algebra::array_soa::detail {

/// Select between two values with a lane bit pattern, without a branch
///
/// @tparam T Value type of the lane.
///
/// @param a Value for an unset lane.
/// @param b Value for a set lane.
/// @param k Lane bit pattern: all bits set or all bits zero.
///
/// @returns @p b if @p k is set, otherwise @p a
template <detray::concepts::value T>
DETRAY_HOST_DEVICE constexpr T blend(const T a, const T b,
                                     const lane_int_t<T> k) {
  using uint_t = lane_uint_t<T>;

  const auto bits_k{static_cast<uint_t>(k)};
  const auto bits_a{std::bit_cast<uint_t>(a)};
  const auto bits_b{std::bit_cast<uint_t>(b)};

  return std::bit_cast<T>(static_cast<uint_t>(
      (bits_a & static_cast<uint_t>(~bits_k)) | (bits_b & bits_k)));
}

/// @tparam S Lane bundle type.
/// @tparam M Lane mask type.
///
/// @param mask Lane mask to convert.
///
/// @returns @p mask with the lane integer size of the values in @tparam S
template <detray::concepts::array_soa_simd S,
          detray::concepts::array_soa_mask M>
DETRAY_HOST_DEVICE constexpr typename S::mask_type to_lane_mask(const M& mask) {
  static_assert(S::size() == M::size(), "Number of lanes does not match");

  if constexpr (std::same_as<M, typename S::mask_type>) {
    return mask;
  } else {
    return typename S::mask_type{mask};
  }
}

}  // namespace detray::algebra::array_soa::detail

namespace detray::algebra::boolean {

/// Lane mask overloads
///
/// The reductions combine the lane bits without branches, so that the
/// compiler can vectorize them
/// @{

/// @tparam M Lane mask type.
/// @param mask Lane mask to check.
/// @returns true if at least one lane is set
template <detray::concepts::array_soa_mask M>
DETRAY_HOST_DEVICE constexpr bool any_of(const M& mask) {
  using lane_t = typename M::lane_type;

  lane_t ret{M::s_false};
  DETRAY_UNROLL_N(M::size())
  for (std::size_t i = 0u; i < M::size(); ++i) {
    ret = static_cast<lane_t>(ret | mask.m_data[i]);
  }
  return ret != M::s_false;
}

/// @tparam M Lane mask type.
/// @param mask Lane mask to check.
/// @returns true if no lane is set
template <detray::concepts::array_soa_mask M>
DETRAY_HOST_DEVICE constexpr bool none_of(const M& mask) {
  return !any_of(mask);
}

/// @tparam M Lane mask type.
/// @param mask Lane mask to check.
/// @returns true if all lanes are set
template <detray::concepts::array_soa_mask M>
DETRAY_HOST_DEVICE constexpr bool all_of(const M& mask) {
  return mask.isFull();
}

/// @tparam M Lane mask type.
/// @param mask Lane mask to check.
/// @returns the number of lanes that are set
template <detray::concepts::array_soa_mask M>
DETRAY_HOST_DEVICE constexpr std::size_t count(const M& mask) {
  std::size_t n{0u};
  DETRAY_UNROLL_N(M::size())
  for (std::size_t i = 0u; i < M::size(); ++i) {
    n += static_cast<std::size_t>(mask.m_data[i] & 1);
  }
  return n;
}
/// @}

/// Array SoA overloads for masked assignments
/// @{

/// Assign a single value to the lanes where the mask is set
///
/// @tparam S Lane bundle type.
/// @tparam M Lane mask type.
/// @tparam U Type of the new value, convertible to the lane value type.
///
/// @param s Lane bundle to update.
/// @param mask Lanes to update.
/// @param new_value Value to assign.
template <detray::concepts::array_soa_simd S,
          detray::concepts::array_soa_mask M, typename U>
  requires(std::is_scalar_v<U> &&
           std::convertible_to<U, typename S::value_type>)
DETRAY_HOST_DEVICE constexpr void set_if(S& s, const M& mask,
                                         const U new_value) {
  using value_t = typename S::value_type;

  const auto k{array_soa::detail::to_lane_mask<S>(mask)};
  const auto v{static_cast<value_t>(new_value)};

  DETRAY_UNROLL_N(S::size())
  for (std::size_t i = 0u; i < S::size(); ++i) {
    s[i] = array_soa::detail::blend(s[i], v, k.m_data[i]);
  }
}

/// Assign the lanes of a lane bundle where the mask is set
///
/// @tparam S Lane bundle type.
/// @tparam M Lane mask type.
/// @tparam U Type of the new values, convertible to @tparam S.
///
/// @param s Lane bundle to update.
/// @param mask Lanes to update.
/// @param new_s Lane bundle with the new values.
template <detray::concepts::array_soa_simd S,
          detray::concepts::array_soa_mask M, typename U>
  requires(!std::is_scalar_v<U> && std::convertible_to<U, S>)
DETRAY_HOST_DEVICE constexpr void set_if(S& s, const M& mask, const U& new_s) {
  using value_t = typename S::value_type;

  const auto k{array_soa::detail::to_lane_mask<S>(mask)};

  DETRAY_UNROLL_N(S::size())
  for (std::size_t i = 0u; i < S::size(); ++i) {
    s[i] = array_soa::detail::blend(s[i], static_cast<value_t>(new_s[i]),
                                    k.m_data[i]);
  }
}

/// Set the lanes to zero where the mask is set
///
/// @tparam S Lane bundle type.
/// @tparam M Lane mask type.
///
/// @param s Lane bundle to update.
/// @param mask Lanes to set to zero.
template <detray::concepts::array_soa_simd S,
          detray::concepts::array_soa_mask M>
DETRAY_HOST_DEVICE constexpr void set_zero(S& s, const M& mask) {
  set_if(s, mask, static_cast<typename S::value_type>(0));
}

/// Set the lanes to zero where the mask is not set
///
/// @tparam S Lane bundle type.
/// @tparam M Lane mask type.
///
/// @param s Lane bundle to update.
/// @param mask Lanes to keep.
template <detray::concepts::array_soa_simd S,
          detray::concepts::array_soa_mask M>
DETRAY_HOST_DEVICE constexpr void set_zero_inverted(S& s, const M& mask) {
  set_if(s, !mask, static_cast<typename S::value_type>(0));
}
/// @}

}  // namespace detray::algebra::boolean
