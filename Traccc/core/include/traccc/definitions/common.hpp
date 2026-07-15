/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"

// Detray include(s).
#include <detray/definitions/units.hpp>

namespace traccc {

template <typename scalar_t>
using unit = detray::unit<scalar_t>;

template <typename scalar_t>
using constant = detray::constant<scalar_t>;

// epsilon for float variables
constexpr scalar float_epsilon = 1e-5f;

}  // namespace traccc
