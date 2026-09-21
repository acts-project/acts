/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/hip/finding/combinatorial_kalman_filter_algorithm.hpp"

namespace traccc::hip {

combinatorial_kalman_filter_algorithm::combinatorial_kalman_filter_algorithm(
    const finding_config& config, const traccc::memory_resource& mr,
    const vecmem::copy& copy, const stream_wrapper& str,
    std::unique_ptr<const Logger> logger,
    std::unique_ptr<traccc::hip::kalman_fitting_algorithm> kf_fitter)
    : device::combinatorial_kalman_filter_algorithm(
          config, mr, copy, std::move(logger), std::move(kf_fitter)),
      hip::algorithm_base(str) {}

void combinatorial_kalman_filter_algorithm::synchronize() const {
  stream().synchronize();
}

}  // namespace traccc::hip
