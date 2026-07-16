/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/math.hpp"

namespace traccc::edm {

template <typename BASE>
TRACCC_HOST_DEVICE auto& spacepoint<BASE>::x() {

    return global()[0];
}

template <typename BASE>
TRACCC_HOST_DEVICE const auto& spacepoint<BASE>::x() const {

    return global()[0];
}

template <typename BASE>
TRACCC_HOST_DEVICE auto& spacepoint<BASE>::y() {

    return global()[1];
}

template <typename BASE>
TRACCC_HOST_DEVICE const auto& spacepoint<BASE>::y() const {

    return global()[1];
}

template <typename BASE>
TRACCC_HOST_DEVICE auto& spacepoint<BASE>::z() {

    return global()[2];
}

template <typename BASE>
TRACCC_HOST_DEVICE const auto& spacepoint<BASE>::z() const {

    return global()[2];
}

template <typename BASE>
TRACCC_HOST_DEVICE auto spacepoint<BASE>::radius() const {

    const float xx = x();
    const float yy = y();
    return math::sqrt(xx * xx + yy * yy);
}

template <typename BASE>
TRACCC_HOST_DEVICE auto spacepoint<BASE>::phi() const {

    return math::atan2(y(), x());
}

template <typename BASE>
template <typename T>
TRACCC_HOST_DEVICE bool spacepoint<BASE>::operator==(
    const spacepoint<T>& other) const {

    return ((measurement_index_1() == other.measurement_index_1()) &&
            (measurement_index_2() == other.measurement_index_2()) &&
            (math::fabs(x() - other.x()) < 1e-6f) &&
            (math::fabs(y() - other.y()) < 1e-6f) &&
            (math::fabs(z() - other.z()) < 1e-6f) &&
            (z_variance() == other.z_variance()) &&
            (radius_variance() == other.radius_variance()));
}

template <typename BASE>
template <typename T>
TRACCC_HOST_DEVICE std::weak_ordering spacepoint<BASE>::operator<=>(
    const spacepoint<T>& other) const {

    return (radius() <=> other.radius());
}

}  // namespace traccc::edm
