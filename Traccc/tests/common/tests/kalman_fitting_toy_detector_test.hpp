/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "kalman_fitting_test.hpp"
#include "test_detectors.hpp"
#include "toy_detector_fixture.hpp"
#include "traccc/io/detector.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <array>

namespace traccc {

/// Kalman Fitting Test with Toy Geometry
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
class KalmanFittingToyDetectorTests
    : public KalmanFittingTests,
      public ToyDetectorFixture,
      public testing::WithParamInterface<std::tuple<
          std::string, std::array<scalar, 3u>, std::array<scalar, 3u>,
          std::array<scalar, 2u>, std::array<scalar, 2u>,
          std::array<scalar, 2u>, detray::pdg_particle<scalar>, unsigned int,
          unsigned int, bool>> {

    public:
    static constexpr std::array<scalar, 2u> smearing{
        50.f * traccc::unit<scalar>::um, 50.f * traccc::unit<scalar>::um};

    static constexpr std::array<double, e_bound_size> stddevs = {
        0.01 * traccc::unit<double>::mm,
        0.01 * traccc::unit<double>::mm,
        0.001,
        0.001,
        0.001 / traccc::unit<double>::GeV,
        0.01 * traccc::unit<double>::ns};
};

}  // namespace traccc
