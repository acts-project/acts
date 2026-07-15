/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/fitting/kalman_fitting_algorithm.hpp"

#include "traccc/bfield/magnetic_field_types.hpp"
#include "traccc/fitting/details/kalman_fitting.hpp"
#include "traccc/fitting/details/kalman_fitting_types.hpp"
#include "traccc/utils/host_detector_bfield_visitor.hpp"

namespace traccc::host {

kalman_fitting_algorithm::kalman_fitting_algorithm(
    const config_type& config, vecmem::memory_resource& mr, vecmem::copy& copy,
    std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), m_config{config}, m_mr{mr}, m_copy(copy) {}

kalman_fitting_algorithm::output_type kalman_fitting_algorithm::operator()(
    const host_detector& det, const magnetic_field& bfield,
    const edm::track_container<default_algebra>::const_view& track_candidates)
    const {

    // Perform the track fitting using the appropriate templated implementation.
    return host_detector_magnetic_field_visitor<detector_type_list,
                                                bfield_type_list<scalar>>(
        det, bfield,
        [&]<typename detector_t, typename bfield_view_t>(
            const typename detector_t::host& detector,
            const bfield_view_t field) {
            traccc::details::kalman_fitter_t<typename detector_t::host,
                                             bfield_view_t>
                fitter{detector, field, m_config};
            return details::kalman_fitting<default_algebra>(
                fitter, track_candidates, m_mr.get(), m_copy.get());
        });
}

}  // namespace traccc::host
