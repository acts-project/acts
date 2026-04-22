// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/algebra/common/boolean.hpp"

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

namespace detray::algebra::boolean {

/// Boolean utilities on single values
/// @{
using algebra::boolean::all_of;
using algebra::boolean::any_of;
using algebra::boolean::none_of;
/// @}

/// Vc overloads of boolean utilities
/// @{
template <typename T>
  requires Vc::Traits::is_simd_mask<T>::value
constexpr bool any_of(T &&mask) {
  return Vc::any_of(std::forward<T>(mask));
}

template <typename T>
  requires Vc::Traits::is_simd_mask<T>::value
constexpr bool all_of(T &&mask) {
  return Vc::all_of(std::forward<T>(mask));
}

template <typename T>
  requires Vc::Traits::is_simd_mask<T>::value
constexpr bool none_of(T &&mask) {
  return Vc::none_of(std::forward<T>(mask));
}
/// @}

}  // namespace detray::algebra::boolean
