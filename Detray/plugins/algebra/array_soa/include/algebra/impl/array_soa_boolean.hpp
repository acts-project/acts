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
using algebra::boolean::all_of;
using algebra::boolean::any_of;
using algebra::boolean::none_of;
/// @}

/// Lane mask overloads
/// @{
template <detray::concepts::array_soa_mask T>
DETRAY_HOST_DEVICE constexpr bool any_of(const T &m) {
  return m.isNotEmpty();
}

template <detray::concepts::array_soa_mask T>
DETRAY_HOST_DEVICE constexpr bool all_of(const T &m) {
  return m.isFull();
}

template <detray::concepts::array_soa_mask T>
DETRAY_HOST_DEVICE constexpr bool none_of(const T &m) {
  return m.isEmpty();
}
/// @}

}  // namespace detray::algebra::boolean
