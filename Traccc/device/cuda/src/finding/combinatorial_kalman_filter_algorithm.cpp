/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"

namespace traccc::cuda {

combinatorial_kalman_filter_algorithm::combinatorial_kalman_filter_algorithm(
    const finding_config& config, const traccc::memory_resource& mr,
    const vecmem::copy& copy, const stream_wrapper& str,
    std::unique_ptr<const Logger> logger,
    std::unique_ptr<traccc::cuda::kalman_fitting_algorithm> kf_fitter)
    : device::combinatorial_kalman_filter_algorithm(
          config, mr, copy, std::move(logger), std::move(kf_fitter)),
      cuda::algorithm_base(str),
      m_expected_layer_config(config) {}

void combinatorial_kalman_filter_algorithm::update_expected_layer_mapping(
    const expected_layer_mapping_entry* map, std::size_t map_size) {
  device::combinatorial_kalman_filter_algorithm::update_expected_layer_mapping(
      map, map_size);
  m_expected_layer_config.expected_layer_map = map;
  m_expected_layer_config.expected_layer_map_size = map_size;
}
}  // namespace traccc::cuda
