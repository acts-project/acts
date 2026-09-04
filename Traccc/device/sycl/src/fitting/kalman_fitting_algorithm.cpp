/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/sycl/fitting/kalman_fitting_algorithm.hpp"

namespace traccc::sycl {

kalman_fitting_algorithm::kalman_fitting_algorithm(
    const config_type& config, const traccc::memory_resource& mr,
    const vecmem::copy& copy, queue_wrapper& q,
    std::unique_ptr<const Logger> logger)
    : device::kalman_fitting_algorithm{config, mr, copy, std::move(logger)},
      sycl::algorithm_base{q} {}

}  // namespace traccc::sycl
