/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::cuda {
/**
 * @brief Calculate the index of the current thread in the set of threads in
 * the warp for which a predicate evaluates to true.
 *
 * @param predicate The predicate to evaluate.
 * @param mask The thread mask to sync.
 *
 * @returns The tuple (T, I), where T is the total number of threads in the
 * warp for which the provided predicate was true, and where I is a unique,
 * consecutive index in the set of threads for which this was true.
 *
 * @example If a predicate P evaluates as `[T, F, F, T]` in a 4-lane warp, the
 * call to `warp_indexed_ballot_sync(P)` will evaluate to the values
 * `[(2, 0), (2, ?), (2, ?), (2, 1)]` where `?` indicates an undefined value.
 *
 * @warning The value of T is always well-defined, the value if I is only
 * well-defined if the given predicate was true for the given thread.
 *
 * @note As with all CUDA synchronization barriers, exited threads are treated
 * as having implicitly reached the barrier to avoid deadlock situations.
 *
 * @note This function forces thread synchronization.
 */
__device__ __forceinline__ std::pair<uint32_t, uint32_t>
warp_indexed_ballot_sync(bool predicate, uint32_t mask = 0xFFFFFFFFu) {
    uint32_t vote = __ballot_sync(mask, predicate);

    /*
     * The total number of threads which return true is simply the population
     * count of the voting result, which is to say the number of true bits in
     * its binary expansion.
     */
    uint32_t tot = __popc(vote);

    /*
     * The index is the population count of the vote mask, but only for bits
     * with an index lower than the current thread index. Thus, we shift a bit
     * mask over the voting result, nulling any bits that are _higher_ than
     * the thread index!
     */
    uint32_t idx = __popc(vote & ~(0xFFFFFFFFu << (threadIdx.x % warpSize)));

    return {tot, idx};
}
}  // namespace traccc::cuda
