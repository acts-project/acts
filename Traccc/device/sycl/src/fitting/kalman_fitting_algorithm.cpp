/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/sycl/fitting/kalman_fitting_algorithm.hpp"

namespace traccc::sycl {

kalman_fitting_algorithm::kalman_fitting_algorithm(
    const config_type& config, const traccc::memory_resource& mr,
    vecmem::copy& copy, queue_wrapper queue,
    std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)),
      m_config{config},
      m_mr{mr},
      m_copy{copy},
      m_queue{queue} {}

}  // namespace traccc::sycl
