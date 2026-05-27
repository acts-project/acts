/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::edm {

template <typename BASE>
template <typename T>
TRACCC_HOST_DEVICE bool seed<BASE>::operator==(const seed<T>& other) const {

    return ((top_index() == other.top_index()) &&
            (middle_index() == other.middle_index()) &&
            (bottom_index() == other.bottom_index()));
}

template <typename BASE>
template <typename T>
TRACCC_HOST_DEVICE std::strong_ordering seed<BASE>::operator<=>(
    const seed<T>& other) const {

    if (top_index() != other.top_index()) {
        return (top_index() <=> other.top_index());
    } else if (middle_index() != other.middle_index()) {
        return (middle_index() <=> other.middle_index());
    } else {
        return (bottom_index() <=> other.bottom_index());
    }
}

}  // namespace traccc::edm
