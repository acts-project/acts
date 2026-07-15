/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/seeding/track_params_estimation.hpp"

#include "traccc/seeding/track_params_estimation_helper.hpp"

// System include(s).
#include <cassert>

namespace traccc::host {

track_params_estimation::track_params_estimation(
    const track_params_estimation_config& config, vecmem::memory_resource& mr,
    std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), m_config(config), m_mr(mr) {}

track_params_estimation::output_type track_params_estimation::operator()(
    const edm::measurement_collection::const_view& measurements_view,
    const edm::spacepoint_collection::const_view& spacepoints_view,
    const edm::seed_collection::const_view& seeds_view,
    const vector3& bfield) const {

    // Set up the input / output objects.
    const edm::measurement_collection::const_device measurements(
        measurements_view);
    const edm::spacepoint_collection::const_device spacepoints(
        spacepoints_view);
    const edm::seed_collection::const_device seeds(seeds_view);
    const edm::seed_collection::const_device::size_type num_seeds =
        seeds.size();
    output_type result(num_seeds, &m_mr.get());

    // Create a track parameters for each seed.
    for (edm::seed_collection::const_device::size_type i = 0; i < num_seeds;
         ++i) {

        TRACCC_VERBOSE("Creating track parameters for seed " << i + 1 << " / "
                                                             << num_seeds);
        TRACCC_VERBOSE("  - bottom spacepoint: "
                       << spacepoints.at(seeds.at(i).bottom_index()).global()
                       << ", radius: "
                       << spacepoints.at(seeds.at(i).bottom_index()).radius());
        TRACCC_VERBOSE("  - middle spacepoint: "
                       << spacepoints.at(seeds.at(i).middle_index()).global()
                       << ", radius: "
                       << spacepoints.at(seeds.at(i).middle_index()).radius());
        TRACCC_VERBOSE("  - top spacepoint: "
                       << spacepoints.at(seeds.at(i).top_index()).global()
                       << ", radius: "
                       << spacepoints.at(seeds.at(i).top_index()).radius());

        // Calculate the track parameter vector.
        bound_track_parameters<>& track_params = result.at(i);
        seed_to_bound_param_vector(track_params, measurements, spacepoints,
                                   seeds.at(i), bfield);

        // Set Covariance
        for (std::size_t j = 0u; j < e_bound_size; ++j) {
            scalar var =
                m_config.initial_sigma.at(j) * m_config.initial_sigma.at(j);

            if (j == e_bound_qoverp) {
                scalar var_theta = getter::element(
                    track_params.covariance(), e_bound_theta, e_bound_theta);

                var += math::pow(m_config.initial_sigma_qopt *
                                     math::sin(track_params.theta()),
                                 2.f);
                var += math::pow(
                    m_config.initial_sigma_pt_rel * track_params.qop(), 2.f);
                var += var_theta *
                       math::pow(
                           track_params.qop() / math::tan(track_params.theta()),
                           2.f);
            }

            var *= m_config.initial_inflation.at(j);

            getter::element(track_params.covariance(), j, j) = var;
        }
    }

    // Return the result.
    return result;
}

}  // namespace traccc::host
