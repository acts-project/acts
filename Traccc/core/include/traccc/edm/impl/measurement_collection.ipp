/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/math.hpp"

namespace traccc::edm {

template <typename BASE>
template <std::integral TYPE>
TRACCC_HOST_DEVICE void measurement<BASE>::set_subspace(
    const std::array<TYPE, 2u>& subs) {

    subspace()[0] = static_cast<std::uint8_t>(subs[0]);
    subspace()[1] = static_cast<std::uint8_t>(subs[1]);
}

template <typename BASE>
template <typename T>
TRACCC_HOST_DEVICE bool measurement<BASE>::operator==(
    const measurement<T>& other) const {

    return ((surface_link() == other.surface_link()) &&
            (math::abs(local_position()[0] - other.local_position()[0]) <
             float_epsilon) &&
            (math::abs(local_position()[1] - other.local_position()[1]) <
             float_epsilon) &&
            (math::abs(local_variance()[0] - other.local_variance()[0]) <
             float_epsilon) &&
            (math::abs(local_variance()[1] - other.local_variance()[1]) <
             float_epsilon));
}

template <typename BASE>
template <typename T>
TRACCC_HOST_DEVICE std::partial_ordering measurement<BASE>::operator<=>(
    const measurement<T>& other) const {

    if (surface_link() != other.surface_link()) {
        return (surface_link() <=> other.surface_link());
    } else if (local_position()[0] != other.local_position()[0]) {
        return (local_position()[0] <=> other.local_position()[0]);
    } else if (local_position()[1] != other.local_position()[1]) {
        return (local_position()[1] <=> other.local_position()[1]);
    } else if (local_variance()[0] != other.local_variance()[0]) {
        return (local_variance()[0] <=> other.local_variance()[0]);
    } else {
        return (local_variance()[1] <=> other.local_variance()[1]);
    }
}

}  // namespace traccc::edm
