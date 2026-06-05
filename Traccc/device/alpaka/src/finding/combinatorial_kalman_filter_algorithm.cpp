/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/alpaka/finding/combinatorial_kalman_filter_algorithm.hpp"

#include "../utils/get_queue.hpp"
#include "../utils/magnetic_field_types.hpp"
#include "combinatorial_kalman_filter.hpp"
#include "traccc/bfield/magnetic_field_types.hpp"
#include "traccc/utils/detector_buffer_bfield_visitor.hpp"

namespace traccc::alpaka {

combinatorial_kalman_filter_algorithm::combinatorial_kalman_filter_algorithm(
    const config_type& config, const traccc::memory_resource& mr,
    vecmem::copy& copy, queue& q, std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)),
      m_config{config},
      m_mr{mr},
      m_copy{copy},
      m_queue{q} {}

combinatorial_kalman_filter_algorithm::output_type
combinatorial_kalman_filter_algorithm::operator()(
    const detector_buffer& det, const magnetic_field& bfield,
    const edm::measurement_collection::const_view& measurements,
    const bound_track_parameters_collection_types::const_view& seeds) const {

    // Perform the track finding using the templated implementation.
    return detector_buffer_magnetic_field_visitor<
        detector_type_list, alpaka::bfield_type_list<scalar>>(
        det, bfield,
        [&]<typename detector_t, typename bfield_view_t>(
            const typename detector_t::view& detector,
            const bfield_view_t& field) {
            return details::combinatorial_kalman_filter<
                typename detector_t::device>(
                detector, field, measurements, seeds, m_config, m_mr, m_copy,
                logger(), details::get_queue(m_queue.get()));
        });
}

}  // namespace traccc::alpaka
