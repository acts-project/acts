/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/alpaka/fitting/kalman_fitting_algorithm.hpp"

#include "../utils/get_queue.hpp"
#include "../utils/magnetic_field_types.hpp"
#include "kalman_fitting.hpp"
#include "traccc/alpaka/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/bfield/magnetic_field_types.hpp"
#include "traccc/utils/detector_buffer_bfield_visitor.hpp"

namespace traccc::alpaka {

kalman_fitting_algorithm::kalman_fitting_algorithm(
    const config_type& config, const traccc::memory_resource& mr,
    vecmem::copy& copy, queue& q, std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)),
      m_config{config},
      m_mr{mr},
      m_copy{copy},
      m_queue{q} {}

kalman_fitting_algorithm::output_type kalman_fitting_algorithm::operator()(
    const detector_buffer& det, const magnetic_field& bfield,
    const edm::track_container<default_algebra>::const_view& track_candidates)
    const {

    // Run the track fitting.
    return detector_buffer_magnetic_field_visitor<
        detector_type_list, alpaka::bfield_type_list<scalar>>(
        det, bfield,
        [&]<typename detector_t, typename bfield_view_t>(
            const typename detector_t::view& detector,
            const bfield_view_t& field) {
            return details::kalman_fitting<typename detector_t::device>(
                detector, field, track_candidates, m_config, m_mr, m_copy.get(),
                details::get_queue(m_queue.get()));
        });
}

}  // namespace traccc::alpaka
