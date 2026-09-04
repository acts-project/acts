/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/sycl/finding/combinatorial_kalman_filter_algorithm.hpp"

namespace traccc::sycl {

combinatorial_kalman_filter_algorithm::combinatorial_kalman_filter_algorithm(
    const config_type& config, const traccc::memory_resource& mr,
    const vecmem::copy& copy, queue_wrapper& queue,
    std::unique_ptr<const Logger> logger,
    std::unique_ptr<traccc::sycl::kalman_fitting_algorithm> kf_fitter)
    : device::combinatorial_kalman_filter_algorithm(
          config, mr, copy, std::move(logger), std::move(kf_fitter)),
      sycl::algorithm_base(queue) {}

}  // namespace traccc::sycl
