/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::edm {

template <typename BASE>
TRACCC_HOST_DEVICE bool track_state<BASE>::is_hole() const {

    return (state() & IS_HOLE_MASK);
}

template <typename BASE>
TRACCC_HOST_DEVICE void track_state<BASE>::set_hole(bool value) {

    if (value) {
        state() |= IS_HOLE_MASK;
    } else {
        state() &= ~IS_HOLE_MASK;
    }
}

template <typename BASE>
TRACCC_HOST_DEVICE bool track_state<BASE>::is_smoothed() const {

    return (state() & IS_SMOOTHED_MASK);
}

template <typename BASE>
TRACCC_HOST_DEVICE void track_state<BASE>::set_smoothed(bool value) {

    if (value) {
        state() |= IS_SMOOTHED_MASK;
    } else {
        state() &= ~IS_SMOOTHED_MASK;
    }
}

}  // namespace traccc::edm
