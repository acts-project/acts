/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/math.hpp"

namespace traccc::edm {

template <typename BASE>
TRACCC_HOST_DEVICE silicon_cell<BASE>& silicon_cell<BASE>::operator=(
    const silicon_cell& other) {

    channel0() = other.channel0();
    channel1() = other.channel1();
    activation() = other.activation();
    time() = other.time();
    module_index() = other.module_index();
    return *this;
}

template <typename BASE>
template <typename T, std::enable_if_t<!std::is_same_v<BASE, T>, bool> >
TRACCC_HOST_DEVICE silicon_cell<BASE>& silicon_cell<BASE>::operator=(
    const silicon_cell<T>& other) {

    channel0() = other.channel0();
    channel1() = other.channel1();
    activation() = other.activation();
    time() = other.time();
    module_index() = other.module_index();
    return *this;
}

template <typename BASE>
template <typename T>
TRACCC_HOST_DEVICE bool silicon_cell<BASE>::operator==(
    const silicon_cell<T>& other) const {

    return (channel0() == other.channel0()) &&
           (channel1() == other.channel1()) &&
           (math::fabs(activation() - other.activation()) < 1e-6f) &&
           (math::fabs(time() - other.time()) < 1e-6f) &&
           (module_index() == other.module_index());
}

template <typename BASE>
template <typename T>
TRACCC_HOST_DEVICE std::weak_ordering silicon_cell<BASE>::operator<=>(
    const silicon_cell<T>& other) const {

    if (module_index() != other.module_index()) {
        return (module_index() <=> other.module_index());
    } else if (channel1() != other.channel1()) {
        return (channel1() <=> other.channel1());
    } else {
        return (channel0() <=> other.channel0());
    }
}

}  // namespace traccc::edm
