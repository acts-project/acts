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

namespace detray {

namespace algebra::boolean {

/// Utilities for single booleans: default case
/// @{
constexpr bool any_of(bool b) {
  return b;
}
constexpr bool all_of(bool b) {
  return b;
}
constexpr bool none_of(bool b) {
  return !b;
}
/// @}

}  // namespace algebra::boolean

namespace detail {

// Pull the boolean functions directly into the detray namespace
using namespace ::detray::algebra::boolean;

}  // namespace detail

}  // namespace detray
