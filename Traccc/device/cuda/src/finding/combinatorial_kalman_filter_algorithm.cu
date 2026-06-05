/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/magnetic_field_types.hpp"
#include "combinatorial_kalman_filter.cuh"
#include "traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/utils/detector_buffer_bfield_visitor.hpp"

// Project include(s).
#include "traccc/bfield/magnetic_field_types.hpp"

// System include(s).
#include <stdexcept>

namespace traccc::cuda {

combinatorial_kalman_filter_algorithm::output_type
combinatorial_kalman_filter_algorithm::operator()(
    const detector_buffer& det, const magnetic_field& field,
    const edm::measurement_collection::const_view& measurements,
    const bound_track_parameters_collection_types::const_view& seeds) const {

    // Perform the track finding using the appropriate templated implementation.
    return detector_buffer_magnetic_field_visitor<
        detector_type_list, cuda::bfield_type_list<scalar>>(
        det, field,
        [&]<typename detector_t, typename bfield_view_t>(
            const typename detector_t::view& detector,
            const bfield_view_t& bfield) {
            return details::combinatorial_kalman_filter<
                typename detector_t::device>(detector, bfield, measurements,
                                             seeds, m_config, m_mr, m_copy,
                                             logger(), m_stream, m_warp_size);
        });
}

}  // namespace traccc::cuda
