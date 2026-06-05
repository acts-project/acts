/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"

#include "traccc/bfield/magnetic_field_types.hpp"
#include "traccc/finding/details/combinatorial_kalman_filter.hpp"
#include "traccc/utils/host_detector_bfield_visitor.hpp"

// System include(s).
#include <stdexcept>

namespace traccc::host {

combinatorial_kalman_filter_algorithm::combinatorial_kalman_filter_algorithm(
    const config_type& config, vecmem::memory_resource& mr,
    std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), m_config{config}, m_mr{mr} {

    // Check the configuration.
    if (m_config.min_track_candidates_per_track == 0) {
        throw std::invalid_argument(
            "The minimum number of track candidates per track must be at least "
            "1.");
    }
}

combinatorial_kalman_filter_algorithm::output_type
combinatorial_kalman_filter_algorithm::operator()(
    const host_detector& det, const magnetic_field& bfield,
    const edm::measurement_collection::const_view& measurements,
    const bound_track_parameters_collection_types::const_view& seeds) const {

    // Perform the track finding using the appropriate templated implementation.
    return host_detector_magnetic_field_visitor<detector_type_list,
                                                bfield_type_list<scalar>>(
        det, bfield,
        [&]<typename detector_t, typename bfield_view_t>(
            const typename detector_t::host& detector,
            const bfield_view_t field) {
            return details::combinatorial_kalman_filter(
                detector, field, measurements, seeds, m_config, m_mr.get(),
                logger());
        });
}

}  // namespace traccc::host
