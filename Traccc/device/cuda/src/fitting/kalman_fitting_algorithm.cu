/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/magnetic_field_types.hpp"
#include "kalman_fitting.cuh"
#include "traccc/cuda/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/utils/detector_buffer_bfield_visitor.hpp"

namespace traccc::cuda {

kalman_fitting_algorithm::output_type kalman_fitting_algorithm::operator()(
    const detector_buffer& det, const magnetic_field& field,
    const edm::track_container<default_algebra>::const_view& track_candidates)
    const {

    // Run the track fitting.
    return detector_buffer_magnetic_field_visitor<
        detector_type_list, cuda::bfield_type_list<scalar>>(
        det, field,
        [&]<typename detector_t, typename bfield_view_t>(
            const typename detector_t::view& detector,
            const bfield_view_t& bfield) {
            return details::kalman_fitting<typename detector_t::device>(
                detector, bfield, track_candidates, m_config, m_mr,
                m_copy.get(), m_stream, m_warp_size);
        });
}

}  // namespace traccc::cuda
