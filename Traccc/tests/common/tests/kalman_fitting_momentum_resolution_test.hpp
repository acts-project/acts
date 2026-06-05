/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "kalman_fitting_test.hpp"
#include "test_detectors.hpp"
#include "traccc/io/detector.hpp"

namespace traccc {

/// Momentum Resolution Test with Telescope Detector
///
/// Test parameters:
/// (1) name
/// (2) origin
/// (3) origin stddev
/// (4) momentum range
/// (5) eta range
/// (6) phi range
/// (7) particle type
/// (8) number of tracks per event
/// (9) number of events
/// (10) random charge
/// (11) offset from origin of the first plane in mm
/// (12) Number of planes
/// (13) Spacing between planes in mm
/// (14) Magnetic field
/// (15) Module material
/// (16) Measurement smearing
class KalmanFittingMomentumResolutionTests
    : public KalmanFittingTests,
      public testing::WithParamInterface<std::tuple<
          std::string, std::array<scalar, 3u>, std::array<scalar, 3u>, scalar,
          scalar, scalar, traccc::pdg_particle<scalar>, unsigned int,
          unsigned int, bool, scalar, unsigned int, scalar, vector3,
          detray::material<scalar>, std::array<scalar, 2u>>> {

    public:
    /// Plane thickness
    static constexpr scalar thickness = 0.5f * traccc::unit<scalar>::mm;

    /// Standard deviations for seed track parameters
    std::array<scalar, e_bound_size> stddevs = {
        5.f * traccc::unit<scalar>::mm,
        5.f * traccc::unit<scalar>::mm,
        0.05f,
        0.05f,
        0.1f / traccc::unit<scalar>::GeV,
        1.f * traccc::unit<scalar>::ns};

    void consistency_tests(
        const edm::track_collection<default_algebra>::host::const_proxy_type&
            track,
        const edm::track_state_collection<default_algebra>::host& track_states)
        const;

    void momentum_resolution_tests(std::string_view file_name) const;

    protected:
    virtual void SetUp() override;
};

}  // namespace traccc
