// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/common/boolean.hpp"
#include "detray/algebra/concepts.hpp"

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

// System include(s).
#include <concepts>
#include <cstddef>
#include <type_traits>

namespace detray::algebra::boolean {

/// Boolean utilities on single values
/// @{
using algebra::boolean::all_of;
using algebra::boolean::any_of;
using algebra::boolean::count;
using algebra::boolean::none_of;
/// @}

/// Vc overloads of boolean utilities
/// @{
/// @tparam M Vc mask type.
///
/// @param mask Lane mask to check.
/// @returns true if at least one lane is set
template <typename M>
  requires Vc::Traits::is_simd_mask<M>::value
constexpr bool any_of(const M& mask) {
  return Vc::any_of(mask);
}

/// @tparam M Vc mask type.
///
/// @param mask Lane mask to check.
/// @returns true if all lanes are set
template <typename M>
  requires Vc::Traits::is_simd_mask<M>::value
constexpr bool all_of(const M& mask) {
  return Vc::all_of(mask);
}

/// @tparam M Vc mask type.
///
/// @param mask Lane mask to check.
/// @returns true if no lane is set
template <typename M>
  requires Vc::Traits::is_simd_mask<M>::value
constexpr bool none_of(const M& mask) {
  return Vc::none_of(mask);
}

/// @tparam M Vc mask type.
///
/// @param mask Lane mask to check.
/// @returns the number of lanes that are set
template <typename M>
  requires Vc::Traits::is_simd_mask<M>::value
constexpr std::size_t count(const M& mask) {
  return mask.count();
}
/// @}

/// Vc overloads for masked assignments
/// @{
/// Assign a single value to the lanes where the mask is set
///
/// @tparam S Vc vector type.
/// @tparam M Vc mask type.
/// @tparam U Type of the new value, convertible to the vector value type.
///
/// @param s Vector to update.
/// @param mask Lanes to update.
/// @param new_value Value to assign.
template <detray::concepts::scalar S, typename M, typename U>
  requires(Vc::Traits::is_simd_vector<S>::value &&
           Vc::Traits::is_simd_mask<M>::value && std::is_scalar_v<U> &&
           std::convertible_to<U, typename S::value_type>)
constexpr void set_if(S& s, const M& mask, const U new_value) {
  s(mask) = new_value;
}

/// Assign the lanes of a vector where the mask is set
///
/// @tparam S Vc vector type.
/// @tparam M Vc mask type.
///
/// @param s Vector to update.
/// @param mask Lanes to update.
/// @param new_s Vector with the new lane values.
template <detray::concepts::scalar S, typename M>
  requires(Vc::Traits::is_simd_vector<S>::value &&
           Vc::Traits::is_simd_mask<M>::value)
constexpr void set_if(S& s, const M& mask, const S& new_s) {
  s(mask) = new_s;
}

/// Set the lanes to zero where the mask is set
///
/// @tparam S Vc vector type.
/// @tparam M Vc mask type.
///
/// @param s Vector to update.
/// @param mask Lanes to set to zero.
template <detray::concepts::scalar S, typename M>
  requires(Vc::Traits::is_simd_vector<S>::value &&
           Vc::Traits::is_simd_mask<M>::value)
constexpr void set_zero(S& s, const M& mask) {
  s.setZero(mask);
}

/// Set the lanes to zero where the mask is not set
///
/// @tparam S Vc vector type.
/// @tparam M Vc mask type.
///
/// @param s Vector to update.
/// @param mask Lanes to keep.
template <detray::concepts::scalar S, typename M>
  requires(Vc::Traits::is_simd_vector<S>::value &&
           Vc::Traits::is_simd_mask<M>::value)
constexpr void set_zero_inverted(S& s, const M& mask) {
  s.setZeroInverted(mask);
}
///@}

}  // namespace detray::algebra::boolean
