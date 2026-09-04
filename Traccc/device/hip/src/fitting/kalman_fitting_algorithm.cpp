/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/hip/fitting/kalman_fitting_algorithm.hpp"

#include "../utils/magnetic_field_types.hpp"
#include "../utils/utils.hpp"

// Project include(s).
#include "traccc/geometry/detector.hpp"

namespace traccc::hip {

kalman_fitting_algorithm::kalman_fitting_algorithm(
    const config_type& config, const traccc::memory_resource& mr,
    const vecmem::copy& copy, const stream_wrapper& str,
    std::unique_ptr<const Logger> logger)
    : device::kalman_fitting_algorithm{config, mr, copy, std::move(logger)},
      hip::algorithm_base{str} {}

auto kalman_fitting_algorithm::prepare_fit_payload(
    const detector_buffer& det, const magnetic_field& field,
    const std::vector<unsigned int>& n_surfaces,
    const device::fit_payload& payload) const -> fit_payload {
  return prepare_fit_payload_helper<detector_type_list,
                                    hip::bfield_type_list<scalar>>(
      det, field, n_surfaces, payload);
}

}  // namespace traccc::hip
