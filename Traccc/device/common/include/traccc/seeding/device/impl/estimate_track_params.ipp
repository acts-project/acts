/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/seeding/detail/track_params_estimation_config.hpp"
#include "traccc/seeding/device/estimate_track_params.hpp"
#include "traccc/seeding/track_params_estimation_helper.hpp"

// System include(s).
#include <cassert>

namespace traccc::device {

template <typename bfield_t>
TRACCC_HOST_DEVICE inline void estimate_track_params(
    const global_index_t globalIndex,
    const track_params_estimation_config& config,
    const edm::measurement_collection::const_view& measurements_view,
    const edm::spacepoint_collection::const_view& spacepoints_view,
    const edm::seed_collection::const_view& seeds_view, const bfield_t& bfield,
    bound_track_parameters_collection_types::view params_view) {

    // Check if anything needs to be done.
    const edm::seed_collection::const_device seeds(seeds_view);
    if (globalIndex >= seeds.size()) {
        return;
    }

    // Create the rest of the device objects.
    const edm::measurement_collection::const_device measurements(
        measurements_view);
    const edm::spacepoint_collection::const_device spacepoints(
        spacepoints_view);
    bound_track_parameters_collection_types::device params(params_view);

    // Figure out the magnetic field at the bottom spacepoint's position.
    const edm::spacepoint bottom_sp =
        spacepoints.at(seeds.bottom_index().at(globalIndex));
    const auto covfie_field_at_sp =
        bfield.at(bottom_sp.x(), bottom_sp.y(), bottom_sp.z());
    const vector3 vector_field_at_sp{
        covfie_field_at_sp[0], covfie_field_at_sp[1], covfie_field_at_sp[2]};

    // Get bound track parameter
    bound_track_parameters<>& track_params = params.at(globalIndex);
    new (&track_params) bound_track_parameters<>();
    seed_to_bound_param_vector(track_params, measurements, spacepoints,
                               seeds.at(globalIndex), vector_field_at_sp);

    // NOTE: The code below uses the covariance of theta in the calculation of
    // the calculation of q/p. Thus, theta must be computed first.
    static_assert(e_bound_qoverp > e_bound_theta);

    // Set Covariance
    for (std::size_t i = 0; i < e_bound_size; i++) {
        scalar var = config.initial_sigma[i] * config.initial_sigma[i];

        if (i == e_bound_qoverp) {
            const scalar var_theta = getter::element(
                track_params.covariance(), e_bound_theta, e_bound_theta);

            // Contribution from sigma(q/pt)
            const scalar sigma_qopt =
                config.initial_sigma_qopt * math::sin(track_params.theta());
            var += sigma_qopt * sigma_qopt;

            // Contribution from sigma(pt)/pt
            const scalar sigma_pt_rel =
                config.initial_sigma_pt_rel * track_params.qop();
            var += sigma_pt_rel * sigma_pt_rel;

            // Contribution from sigma(theta)
            scalar sigma_theta =
                track_params.qop() / math::tan(track_params.theta());
            var += var_theta * sigma_theta * sigma_theta;
        }

        var *= config.initial_inflation[i];

        getter::element(track_params.covariance(), i, i) = var;
    }
}

}  // namespace traccc::device
