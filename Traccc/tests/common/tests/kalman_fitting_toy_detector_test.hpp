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
      public testing::WithParamInterface<std::tuple<
          std::string, std::array<scalar, 3u>, std::array<scalar, 3u>,
          std::array<scalar, 2u>, std::array<scalar, 2u>,
          std::array<scalar, 2u>, traccc::pdg_particle<scalar>, unsigned int,
          unsigned int, bool>> {

    public:
    /// Number of barrel layers
    static constexpr inline unsigned int n_barrels{4u};

    /// Number of endcap layers
    static constexpr inline unsigned int n_endcaps{7u};

    /// B field value and its type
    static constexpr vector3 B{0, 0, 2 * traccc::unit<scalar>::T};

    /// Step constraint
    static const inline scalar step_constraint = 1.f * traccc::unit<scalar>::mm;

    /// Measurement smearing parameters
    static constexpr std::array<scalar, 2u> smearing{
        50.f * traccc::unit<scalar>::um, 50.f * traccc::unit<scalar>::um};

    // Grid search window
    static const inline std::array<detray::dindex, 2> search_window{3u, 3u};

    /// Standard deviations for seed track parameters
    static constexpr std::array<scalar, e_bound_size> stddevs = {
        0.01f * traccc::unit<scalar>::mm,
        0.01f * traccc::unit<scalar>::mm,
        0.001f,
        0.001f,
        0.001f / traccc::unit<scalar>::GeV,
        0.01f * traccc::unit<scalar>::ns};

    protected:
    virtual void SetUp() override {
        vecmem::host_memory_resource host_mr;

        detray::toy_det_config<scalar> toy_cfg{};
        toy_cfg.n_brl_layers(n_barrels).n_edc_layers(n_endcaps).do_check(false);

        // Create the toy geometry
        auto [det, name_map] =
            detray::build_toy_detector<traccc::default_algebra>(host_mr,
                                                                toy_cfg);

        // Write detector file
        auto writer_cfg = detray::io::detector_writer_config{}
                              .format(detray::io::format::json)
                              .replace_files(true)
                              .write_grids(true)
                              .write_material(true)
                              .path(std::get<0>(GetParam()));
        detray::io::write_detector(det, name_map, writer_cfg);
    }
};

}  // namespace traccc
