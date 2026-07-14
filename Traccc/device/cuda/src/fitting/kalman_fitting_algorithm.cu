/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/utils.hpp"
#include "./kernels/fill_fitting_sort_keys.hpp"
#include "traccc/cuda/fitting/kalman_fitting_algorithm.hpp"

// Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/sort.h>

// System include(s).
#include <cassert>
#include <memory_resource>

namespace traccc::cuda {

void kalman_fitting_algorithm::prepare_track_fit_order(
    const edm::track_collection<default_algebra>::const_view& tracks,
    vecmem::data::vector_view<device::sort_key>& track_sort_keys,
    vecmem::data::vector_view<unsigned int>& track_indices) const {

    // Get the number of tracks.
    const unsigned int n_tracks = tracks.capacity();
    assert(n_tracks == copy().get_size(tracks));
    assert(n_tracks == track_indices.capacity());
    assert(track_indices.size_ptr() == nullptr);

    // Launch parameters for the kernel.
    const unsigned int nThreads = warp_size() * 4;
    const unsigned int nBlocks = (n_tracks + nThreads - 1) / nThreads;

    // Fill the keys and indices buffers.
    fill_fitting_sort_keys(nBlocks, nThreads, details::get_stream(stream()),
                           tracks, track_sort_keys, track_indices);

    // Sort the key to get the sorted parameter ids
    vecmem::device_vector<device::sort_key> keys_device(track_sort_keys);
    vecmem::device_vector<unsigned int> track_indices_device(track_indices);
    thrust::sort_by_key(
        thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&mr().main))
            .on(details::get_stream(stream())),
        keys_device.begin(), keys_device.end(), track_indices_device.begin());
}

}  // namespace traccc::cuda
