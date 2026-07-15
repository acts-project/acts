/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <concepts>

namespace traccc::device::concepts {
/**
 * @brief Concept to ensure that a type is barrier-like, i.e. it can
 * synchronize threads in a compute block.
 *
 * @tparam T The barrier-like type.
 */
template <typename T>
concept barrier = requires(T& b) {
    /*
     * Check for the general, nulary barrier function which simply synchronizes
     * threads without return value.
     */
    { b.blockBarrier() } -> std::same_as<void>;

    /*
     * Check for the unary boolean-argument synchronization functions.
     */
    requires requires(bool p) {
        /*
         * `blockOr` should return true iff one or more of the threads in the
         * block issues the call with a truthful argument.
         */
        { b.blockOr(p) } -> std::same_as<bool>;

        /*
         * `blockAnd` should return true iff all of the threads in the block
         * issue the call with a truthful argument.
         */
        { b.blockAnd(p) } -> std::same_as<bool>;

        /*
         * `blockCount` should return the number of threads subject to the
         * barrier that issued the call with a truthful argument.
         */
        { b.blockCount(p) } -> std::integral;
    };
};
}  // namespace traccc::device::concepts
