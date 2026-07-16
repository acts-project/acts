/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cstdint>
#include <vecmem/memory/device_atomic_ref.hpp>

#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/concepts/barrier.hpp"
#include "traccc/device/concepts/thread_id.hpp"

namespace traccc::device {
/**
 * @brief Swap two values of arbitrary type.
 *
 * @tparam T The type of values to swap.
 *
 * @param a The first object in the swap (will take the value of b).
 * @param b The second object in the swap (will take the value of a).
 */
template <std::movable T>
TRACCC_DEVICE void swap(T& a, T& b) {
    T t = std::move(a);
    a = std::move(b);
    b = std::move(t);
}

/**
 * @brief Perform a block-wide odd-even key sorting.
 *
 * This function performs a sorting operation across the entire block, assuming
 * that all the threads in the block are currently active.
 *
 * @warning The behaviour of this function is ill-defined if any of the threads
 * in the block have exited.
 *
 * @warning This method is efficient for sorting small arrays, preferably in
 * shared memory, but given the O(n^2) worst-case performance this should not
 * be used on larger arrays.
 *
 * @tparam T The thread identifier type.
 * @tparam B The barrier type
 * @tparam K The type of keys to sort.
 * @tparam C The type of the comparison function.
 *
 * @param thread_id The thread identifier object.
 * @param barrier The barrier to use for block synchronization.
 * @param keys An array of keys to sort.
 * @param num_keys The number of keys in the array to sort.
 * @param comparison A comparison function.
 */
template <concepts::thread_id1 T, concepts::barrier B, std::movable K,
          std::strict_weak_order<K, K> C>
TRACCC_DEVICE void blockOddEvenSort(const T& thread_id, const B& barrier,
                                    K* keys, uint32_t num_keys,
                                    C&& comparison) {
    bool sorted;

    do {
        sorted = true;

        for (uint32_t j =
                 2 * static_cast<uint32_t>(thread_id.getLocalThreadIdX()) + 1;
             j < num_keys - 1; j += 2 * thread_id.getBlockDimX()) {
            if (comparison(keys[j + 1], keys[j])) {
                swap(keys[j + 1], keys[j]);
                sorted = false;
            }
        }

        barrier.blockBarrier();

        for (uint32_t j =
                 2 * static_cast<uint32_t>(thread_id.getLocalThreadIdX());
             j < num_keys - 1; j += 2 * thread_id.getBlockDimX()) {
            if (comparison(keys[j + 1], keys[j])) {
                swap(keys[j + 1], keys[j]);
                sorted = false;
            }
        }
    } while (barrier.blockOr(!sorted));
}
}  // namespace traccc::device
