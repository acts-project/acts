/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cassert>
#include <cstdint>
#include <tuple>

#include "traccc/definitions/qualifiers.hpp"

namespace traccc::device {
/**
 * @brief Encode the state of our parameter insertion mutex.
 */
TRACCC_HOST_DEVICE inline uint64_t encode_insertion_mutex(const bool locked,
                                                          const uint32_t size,
                                                          const float max) {
    // Assert that the MSB of the size is zero
    assert(size <= 0x7FFFFFFF);

    const uint32_t hi = size | (locked ? 0x80000000 : 0x0);
    const uint32_t lo = std::bit_cast<uint32_t>(max);

    return (static_cast<uint64_t>(hi) << 32) | lo;
}

/**
 * @brief Decode the state of our parameter insertion mutex.
 */
TRACCC_HOST_DEVICE inline std::tuple<bool, uint32_t, float>
decode_insertion_mutex(const uint64_t val) {
    const uint32_t hi = static_cast<uint32_t>(val >> 32);
    const uint32_t lo = val & 0xFFFFFFFF;

    return {static_cast<bool>(hi & 0x80000000), (hi & 0x7FFFFFFF),
            std::bit_cast<float>(lo)};
}
}  // namespace traccc::device
