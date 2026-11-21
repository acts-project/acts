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
///
/// @note The functions with a loop take the number of lanes @c W as a template
/// parameter: GCC 14 only accepts a template parameter as the count of
/// @c DETRAY_UNROLL_N, but no other constant expression.
/// @{

/// @tparam W Number of lanes.
/// @tparam I Signed integer type of a lane.
/// @param mask Lane mask to check.
/// @returns true if at least one lane is set
template <std::size_t W, std::signed_integral I>
DETRAY_HOST_DEVICE constexpr bool any_of(const array_soa::mask<W, I>& mask) {
  using mask_t = array_soa::mask<W, I>;

  I ret{mask_t::s_false};
  DETRAY_UNROLL_N(W)
  for (std::size_t i = 0u; i < W; ++i) {
    ret = static_cast<I>(ret | mask.m_data[i]);
  }
  return ret != mask_t::s_false;
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
  return mask.is_full();
}

/// @tparam W Number of lanes.
/// @tparam I Signed integer type of a lane.
/// @param mask Lane mask to check.
/// @returns the number of lanes that are set
template <std::size_t W, std::signed_integral I>
DETRAY_HOST_DEVICE constexpr std::size_t count(
    const array_soa::mask<W, I>& mask) {
  std::size_t n{0u};
  DETRAY_UNROLL_N(W)
  for (std::size_t i = 0u; i < W; ++i) {
    n += static_cast<std::size_t>(mask.m_data[i] & 1);
  }
  return n;
}
/// @}

/// Array SoA overloads for masked assignments
/// @{

/// Assign a single value to the lanes where the mask is set
///
/// @tparam T Value type of a lane.
/// @tparam W Number of lanes.
/// @tparam I Signed integer type of a mask lane.
/// @tparam U Type of the new value, convertible to the lane value type.
///
/// @param s Lane bundle to update.
/// @param mask Lanes to update.
/// @param new_value Value to assign.
template <detray::concepts::value T, std::size_t W, std::signed_integral I,
          typename U>
  requires(std::is_scalar_v<U> && std::convertible_to<U, T>)
DETRAY_HOST_DEVICE constexpr void set_if(array_soa::simd<T, W>& s,
                                         const array_soa::mask<W, I>& mask,
                                         const U new_value) {
  const auto k{array_soa::detail::to_lane_mask<array_soa::simd<T, W>>(mask)};
  const auto v{static_cast<T>(new_value)};

  DETRAY_UNROLL_N(W)
  for (std::size_t i = 0u; i < W; ++i) {
    s[i] = array_soa::detail::blend(s[i], v, k.m_data[i]);
  }
}

/// Assign the lanes of a lane bundle where the mask is set
///
/// @tparam T Value type of a lane.
/// @tparam W Number of lanes.
/// @tparam I Signed integer type of a mask lane.
/// @tparam U Type of the new values, convertible to the lane bundle.
///
/// @param s Lane bundle to update.
/// @param mask Lanes to update.
/// @param new_s Lane bundle with the new values.
template <detray::concepts::value T, std::size_t W, std::signed_integral I,
          typename U>
  requires(!std::is_scalar_v<U> &&
           std::convertible_to<U, array_soa::simd<T, W>>)
DETRAY_HOST_DEVICE constexpr void set_if(array_soa::simd<T, W>& s,
                                         const array_soa::mask<W, I>& mask,
                                         const U& new_s) {
  const auto k{array_soa::detail::to_lane_mask<array_soa::simd<T, W>>(mask)};

  DETRAY_UNROLL_N(W)
  for (std::size_t i = 0u; i < W; ++i) {
    s[i] =
        array_soa::detail::blend(s[i], static_cast<T>(new_s[i]), k.m_data[i]);
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
