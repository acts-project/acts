/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/track_container.hpp"
#include "traccc/finding/details/combinatorial_kalman_filter_types.hpp"
#include "traccc/finding/details/expected_layer_pattern_with_extrapolation.hpp"
#include "traccc/finding/finding_config.hpp"

// CUDA include(s).
#include <cuda_runtime.h>

namespace traccc::cuda {
namespace kernels {

template <typename detector_t, typename bfield_t>
__global__ void collect_expected_layer_patterns_kernel(
    const typename detector_t::const_view_type det_data,
    const bfield_t field_data,
    const typename edm::track_container<
        typename detector_t::algebra_type>::const_view tracks_view,
    const finding_config config,
    const vecmem::data::vector_view<const expected_layer_mapping_entry>
        expected_layer_map_view,
    vecmem::data::vector_view<expected_layer_pattern_type>
        output_expected_layer_patterns_view);

}  // namespace kernels

template <typename detector_t, typename bfield_t>
inline void collect_expected_layer_patterns(
    const dim3& grid_size, const dim3& block_size, std::size_t shared_mem_size,
    const cudaStream_t& stream,
    const typename detector_t::const_view_type det_data,
    const bfield_t field_data,
    const typename edm::track_container<
        typename detector_t::algebra_type>::const_view tracks_view,
    const finding_config& config,
    const vecmem::data::vector_view<const expected_layer_mapping_entry>&
        expected_layer_map_view,
    vecmem::data::vector_view<expected_layer_pattern_type>
        output_expected_layer_patterns_view);

}  // namespace traccc::cuda

#include "collect_expected_layer_patterns.cu"
