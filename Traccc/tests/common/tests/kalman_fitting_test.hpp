/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/bfield/magnetic_field_types.hpp"
#include "traccc/definitions/common.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_collection.hpp"
#include "traccc/edm/track_state_collection.hpp"
#include "traccc/fitting/kalman_filter/kalman_fitter.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/simulation/event_generators.hpp"
#include "traccc/utils/propagation.hpp"

// Covfie include(s).
#include <covfie/core/field.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <string>
#include <string_view>
#include <vector>

namespace traccc {

class KalmanFittingTests : public testing::Test {
    public:
    /// Type declarations
    using detector_traits = traccc::default_detector;
    using host_detector_type = detector_traits::host;
    using device_detector_type = detector_traits::device;

    using scalar_type = device_detector_type::scalar_type;
    using b_field_t =
        covfie::field<traccc::const_bfield_backend_t<scalar_type>>;
    using rk_stepper_type =
        detray::rk_stepper<b_field_t::view_t, traccc::default_algebra,
                           detray::constrained_step<scalar_type>>;
    using host_navigator_type =
        detray::caching_navigator<const host_detector_type>;
    using host_fitter_type =
        kalman_fitter<rk_stepper_type, host_navigator_type>;
    using device_navigator_type =
        detray::caching_navigator<const device_detector_type>;
    using device_fitter_type =
        kalman_fitter<rk_stepper_type, device_navigator_type>;

    // Use deterministic random number generator for testing
    using uniform_gen_t =
        detray::detail::random_numbers<scalar,
                                       std::uniform_real_distribution<scalar>>;

    /// Verify that pull distribtions follow the normal distribution
    ///
    /// @param file_name The name of the file holding the distributions
    /// @param hist_names The names of the histograms to process
    ///
    void pull_value_tests(std::string_view file_name,
                          const std::vector<std::string>& hist_names) const;

    /// Verify that P value distribtions follow the uniform
    ///
    /// @param file_name The name of the file holding the distributions
    ///
    void p_value_tests(std::string_view file_name) const;

    /// Validadte the NDF for track fitting output
    ///
    /// @param track Fitting statistics result of a track
    /// @param track_states All track states in the event
    /// @param measurements All measurements in the event
    ///
    void ndf_tests(
        const edm::track_collection<default_algebra>::host::const_proxy_type&
            track,
        const edm::track_state_collection<default_algebra>::host& track_states,
        const edm::measurement_collection::host& measurements);

    /// Count the number of tracks that were successfully fitted
    ///
    /// @param tracks The track fit collection to count on
    /// @return The number of tracks that were successfully fitted
    ///
    std::size_t count_successfully_fitted_tracks(
        const edm::track_collection<default_algebra>::host& tracks) const;

    // The number of tracks successful with KF
    std::size_t n_success{0u};
};

}  // namespace traccc
