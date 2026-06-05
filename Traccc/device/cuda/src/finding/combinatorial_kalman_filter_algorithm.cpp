/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"

#include "../utils/magnetic_field_types.hpp"
#include "../utils/utils.hpp"

namespace traccc::cuda {

combinatorial_kalman_filter_algorithm::combinatorial_kalman_filter_algorithm(
    const config_type& config, const traccc::memory_resource& mr,
    vecmem::copy& copy, stream& str, std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)),
      m_config{config},
      m_mr{mr},
      m_copy{copy},
      m_stream{str},
      m_warp_size(details::get_warp_size(str.device())) {}

}  // namespace traccc::cuda
