/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"

namespace traccc::device {

/// @brief Identity projector functor
///
/// This comes in handy when using an SoA collection with some code that was
/// originally designed for AoS collections.
///
struct identity_projector {

    /// Return the input value unchanged
    template <typename T>
    TRACCC_HOST_DEVICE constexpr T operator()(const T& value) const noexcept {
        return value;
    }

};  // struct identity_projector

}  // namespace traccc::device
