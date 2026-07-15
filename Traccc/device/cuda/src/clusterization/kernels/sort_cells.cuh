/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cub/cub.cuh>

#include "traccc/device/sort.hpp"

namespace traccc::cuda {
namespace kernels {
/**
 * @brief Sort the cells that will serve as the input to the clustering
 * algorithm.
 *
 * @note This performs a sort per module, and while the sorting might
 * reorder modules within the input, there is no guarantee that the indices
 * of the modules will be sorted in any particular way.
 *
 * @note The permutation map may be filled even if `write_permutation` is
 * false.
 *
 * @tparam SORT_THREADS_PER_BLOCK The number of threads to use per block when
 *     launching this kernel.
 * @tparam PER_THREAD_MAX The maximum number of cells to sort per thread.
 *
 * @param[in] cells_per_thread The target number of cells to allocate per
 *     thread.
 * @param[in] num_cells The number of cells to sort.
 * @param[in] cells The original, potentially unordered cells.
 * @param[out] new_cells The new, sorted vector of cells.
 * @param[out] permutation_map_view The permutation map between new and old
 *     indices.
 * @param[in] write_permutation Boolean flag to enforce writing to the
 *     `permutation_map_view` vector.
 */
template <std::size_t SORT_THREADS_PER_BLOCK, std::size_t PER_THREAD_MAX>
__global__ __launch_bounds__(SORT_THREADS_PER_BLOCK) void sort_cells(
    const unsigned int cells_per_thread, const unsigned int num_cells,
    const edm::silicon_cell_collection::const_view cells,
    edm::silicon_cell_collection::view new_cells,
    vecmem::data::vector_view<unsigned int> permutation_map_view,
    const bool write_permutation) {

    /*
     * We start off defining some types which we will use in sorting. In
     * particular, we will be doing a key-value sort where the key is composed
     * of the channel0 and channel1 coordinates, as well as the module index.
     * Although this can, for many use cases, fit into 32 bits, we use 64 bits
     * to be sure. The value, then, is simply the index of the cell within the
     * block-local frame, which easily fits into 32 bits.
     */
    using key_t = std::uint64_t;
    using value_t = std::uint32_t;
    static_assert(sizeof(key_t) == 8);
    static_assert(sizeof(value_t) == 4);

    /*
     * The meat of the pudding here is CUB's block radix sort algorithm, which
     * is an exact fit for the job we are trying to do here. This struct
     * defines the thread allocation as well as the shared memory sizing.
     */
    using sort_t = cub::BlockRadixSort<key_t, SORT_THREADS_PER_BLOCK,
                                       PER_THREAD_MAX, value_t>;

    /*
     * Compute the partition target size based on the cells per thread
     * parameter and the block dimension. Note that we routinely overshoot
     * or undershoot this value because the partitioning algorithm is greedy.
     */
    unsigned int partition_target_size = cells_per_thread * blockDim.x;

    /*
     * Define a sum type of the storage that we will need in shared memory,
     * which will hold either the CUB sorting scratch space or a small cache
     * of cell coordinates and module indices. Throughout the computation, we
     * use the software cache first and then evict it to repurpose it as the
     * CUB scratchpad.
     *
     * Note that the cache always has size SORT_THREADS_PER_BLOCK times
     * PER_THREAD_MAX, which is also the maximum partition size that can be
     * sorted with the CUB algorithm.
     */
    union smem_union {
        typename sort_t::TempStorage cub_sorting_tmp;
        struct {
            unsigned short channel0[SORT_THREADS_PER_BLOCK * PER_THREAD_MAX];
            unsigned short channel1[SORT_THREADS_PER_BLOCK * PER_THREAD_MAX];
            unsigned int module_index[SORT_THREADS_PER_BLOCK * PER_THREAD_MAX];
        } cell_attribute_cache;
    };

    /*
     * Shared memory allocations. The union type is by far the largest, but we
     * also keep some integers for administrative purposes.
     */
    __shared__ smem_union smem;
    __shared__ unsigned int partition_start, partition_end;
    __shared__ unsigned int min_module_index, max_module_index, max_channel0,
        max_channel1;

    vecmem::device_vector<unsigned int> permutation_map(permutation_map_view);
    const edm::silicon_cell_collection::const_device cells_device(cells);
    edm::silicon_cell_collection::device new_cells_device(new_cells);

    /*
     * Compute the beginning and end of the partition boundary search. Note
     * that this employs a greedy algorithm similar to the CCL algorithm
     * (albeit with different partition points). The actual partition edges
     * will always be equal to or greater than these, never smaller.
     */
    const unsigned int search_start = blockIdx.x * partition_target_size;
    const unsigned int search_end =
        std::min(num_cells, search_start + partition_target_size);

    if (threadIdx.x == 0) {
        /*
         * We'll find the partition point using a parallel algorithm and rely
         * on atomic min operations, so we initialize the partition
         * boundaries by some sentinel values.
         */
        partition_start = std::numeric_limits<unsigned int>::max();
        partition_end = num_cells;

        /*
         * Similarly, the minimum and maximum module indices and channels are
         * initialized to reasonable sentinels.
         */
        min_module_index = std::numeric_limits<unsigned int>::max();
        max_module_index = 0;
        max_channel0 = 0;
        max_channel1 = 0;
    }

    __syncthreads();

    /*
     * First, perform the boundary search for the upper limit of the
     * partition. Simple parallel marching algorithm which modifies the
     * boundary using an atomic minimum operation in shared memory. Computing
     * the end first allows us to find a reasonable upper bound on the start.
     */
    {
        unsigned int i = search_end + threadIdx.x;
        bool valid_split = false;

        do {
            valid_split = (i > 0 && i < num_cells)
                              ? cells_device.module_index().at(i - 1) !=
                                    cells_device.module_index().at(i)
                              : true;

            if (valid_split) {
                atomicMin(&partition_end, i);
            }

            i += blockDim.x;
        } while (!__syncthreads_or(valid_split));
    }

    /*
     * And then the same algorithm to find a reasonable cut-off point for the
     * start of the partition.
     *
     * NOTE: No synchronization is necessary here, the previous do-while loop
     * synchronizes as part of its loop condition.
     */
    {
        unsigned int i = search_start + threadIdx.x;
        bool valid_split = false;

        do {
            valid_split = (i > 0 && i < partition_end)
                              ? (cells_device.module_index().at(i - 1) !=
                                 cells_device.module_index().at(i))
                              : true;

            if (valid_split) {
                atomicMin(&partition_start, i);
            }

            i += blockDim.x;
        } while (!__syncthreads_or(valid_split));
    }

    __syncthreads();

    /*
     * With the partition boundaries found, compute the partition size and
     * exit if it is 0, in which case there is no sorting to be done.
     */
    assert(partition_start <= partition_end);

    const unsigned int partition_size = partition_end - partition_start;

    if (partition_size == 0) {
        return;
    }

    /*
     * Recall that SORT_THREADS_PER_BLOCK * PER_THREAD_MAX is the maximum
     * number of cells that we can sort in one block, so we need a slow path
     * for larger partitions.
     */
    if (partition_size > SORT_THREADS_PER_BLOCK * PER_THREAD_MAX) [[unlikely]] {
        /*
         * We are unlikely to hit the slow path, but if we do we _need_ the
         * permutation map even if we don't explicitly need it as an output,
         * so we populate its values.
         */
        for (unsigned int i = partition_start + threadIdx.x; i < partition_end;
             i += blockDim.x) {
            permutation_map.at(i) = i;
        }

        __syncthreads();

        /*
         * Simple comparison function of cells that will put them in the
         * desired order.
         */
        auto compare_cells = [&cells_device](const unsigned int idx_a,
                                             const unsigned int idx_b) {
            const auto& a = cells_device.at(idx_a);
            const auto& b = cells_device.at(idx_b);
            if (a.module_index() < b.module_index()) {
                return true;
            } else if (a.module_index() == b.module_index()) {
                if (a.channel1() < b.channel1()) {
                    return true;
                } else if (a.channel1() == b.channel1()) {
                    return a.channel0() < b.channel0();
                } else {
                    return false;
                }

            } else {
                return false;
            }
        };

        /*
         * Launch the (slow) odd-even sort to sort the indices.
         */
        traccc::cuda::details::thread_id1 thread_id;
        traccc::cuda::barrier barrier;
        device::blockOddEvenSort(thread_id, barrier,
                                 &permutation_map.at(partition_start),
                                 partition_size, compare_cells);

        /*
         * Using the permutation map which we have now written, put the cells
         * into the output vector.
         */
        for (unsigned int i = partition_start + threadIdx.x; i < partition_end;
             i += blockDim.x) {
            new_cells_device.at(i) = cells_device.at(permutation_map.at(i));
        }
    } else {
        /*
         * Then the fast path. Here, we will compute the maximum channel0 and
         * channel1 values, as well as the range of module indices, to find
         * the number of bits we need to radix sort over.
         *
         * We first compute the limits per thread to reduce atomic contention.
         */
        unsigned int local_max_channel_0 = 0;
        unsigned int local_max_channel_1 = 0;
        unsigned int local_min_module_index =
            std::numeric_limits<unsigned int>::max();
        unsigned int local_max_module_index = 0;

        /*
         * Perform a parallel march over the partition, recording the values
         * into our shared memory cache and then computing the minima and
         * maxima per thread.
         */
        for (unsigned int i = threadIdx.x; i < partition_size;
             i += blockDim.x) {
            const auto& cell = cells_device.at(partition_start + i);
            const auto channel0 = cell.channel0();
            const auto channel1 = cell.channel1();
            const auto module_index = cell.module_index();

            local_max_channel_0 = std::max(local_max_channel_0, channel0);
            local_max_channel_1 = std::max(local_max_channel_1, channel1);
            local_min_module_index =
                std::min(local_min_module_index, module_index);
            local_max_module_index =
                std::max(local_max_module_index, module_index);

            assert(channel0 <= 0xFFFF);
            smem.cell_attribute_cache.channel0[i] = channel0;
            assert(channel1 <= 0xFFFF);
            smem.cell_attribute_cache.channel1[i] = channel1;
            smem.cell_attribute_cache.module_index[i] = module_index;
        }

        /*
         * Aggregate the boundary values to find the partition-wide minima and
         * maxima.
         */
        atomicMax(&max_channel0, local_max_channel_0);
        atomicMax(&max_channel1, local_max_channel_1);
        atomicMin(&min_module_index, local_min_module_index);
        atomicMax(&max_module_index, local_max_module_index);

        __syncthreads();

        /*
         * Use the CLZ (count leading zero) operation to find the number of
         * bits required to store the range of module indices, channel0, and
         * channel1 values. Then add these up to compute the total key size.
         */
        const auto module_index_clz =
            __clz(max_module_index - min_module_index);
        const unsigned int module_index_bits =
            CHAR_BIT * sizeof(std::decay_t<decltype(module_index_clz)>) -
            module_index_clz;

        const auto channel0_clz = __clz(max_channel0);
        const unsigned int channel0_bits =
            CHAR_BIT * sizeof(std::decay_t<decltype(channel0_clz)>) -
            channel0_clz;

        const auto channel1_clz = __clz(max_channel1);
        const unsigned int channel1_bits =
            CHAR_BIT * sizeof(std::decay_t<decltype(channel1_clz)>) -
            channel1_clz;

        const unsigned int total_bits =
            module_index_bits + channel0_bits + channel1_bits;

        /*
         * Accounting for a single sentinel bit, these bits should fit into
         * the chosen key type.
         *
         * Note that the added bit (which any valid cell cannot write to and
         * leaves at 0) serves as the sentinel bit; only invalid cells set
         * this to 1 and are thus automatically sorted past valid cells.
         */
        assert((total_bits + 1) <= (CHAR_BIT * sizeof(key_t)));

        /*
         * Now, extract the cell values from the shared memory cache and use
         * them to populate thread-local arrays.
         */
        key_t keys[PER_THREAD_MAX];
        value_t values[PER_THREAD_MAX];

        for (unsigned int i = 0; i < PER_THREAD_MAX; ++i) {
            unsigned int eff = i * blockDim.x + threadIdx.x;
            if (eff < partition_size) {
                const auto module_index =
                    smem.cell_attribute_cache.module_index[eff] -
                    min_module_index;
                const auto channel0 = smem.cell_attribute_cache.channel0[eff];
                const auto channel1 = smem.cell_attribute_cache.channel1[eff];

                keys[i] = static_cast<key_t>(module_index)
                              << (channel0_bits + channel1_bits) |
                          static_cast<key_t>(channel1) << (channel0_bits) |
                          static_cast<key_t>(channel0);
                values[i] = eff;
            } else {
                keys[i] = std::numeric_limits<key_t>::max();
                values[i] = std::numeric_limits<value_t>::max();
            }
        }

        __syncthreads();

        /*
         * Finally, perform the sorting on the keys and values using an
         * algorithm provided by CUB.
         */
        sort_t(smem.cub_sorting_tmp)
            .SortBlockedToStriped(keys, values, 0, total_bits + 1);

        /*
         * Now, with the values (original indices) sorted, we can populate a
         * new, sorted cell vector.
         */
        for (unsigned int i = 0; i < PER_THREAD_MAX; ++i) {
            unsigned int eff = i * blockDim.x + threadIdx.x;
            if (eff < partition_size) {
                new_cells_device.at(partition_start + eff) =
                    cells_device.at(partition_start + values[i]);
                if (write_permutation) {
                    permutation_map.at(partition_start + eff) =
                        partition_start + values[i];
                }
            }
        }
    }
}
}  // namespace kernels
}  // namespace traccc::cuda
