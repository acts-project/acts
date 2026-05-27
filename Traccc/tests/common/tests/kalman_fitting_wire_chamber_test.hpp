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

/// Kalman Fitting Test with Wire Chamber
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
class KalmanFittingWireChamberTests
    : public KalmanFittingTests,
      public testing::WithParamInterface<std::tuple<
          std::string, std::array<scalar, 3u>, std::array<scalar, 3u>,
          std::array<scalar, 2u>, std::array<scalar, 2u>,
          std::array<scalar, 2u>, traccc::pdg_particle<scalar>, unsigned int,
          unsigned int, bool>> {

    public:
    /// Number of layers
    static const inline unsigned int n_wire_layers{20u};

    /// Half z of cylinder
    static const inline scalar half_z{2000.f * traccc::unit<scalar>::mm};

    /// B field value and its type
    static constexpr vector3 B{0, 0, 2 * traccc::unit<scalar>::T};

    // Set mask tolerance to a large value not to miss the surface during KF
    static const inline scalar mask_tolerance =
        250.f * traccc::unit<scalar>::um;

    // Grid search window
    static const inline std::array<detray::dindex, 2> search_window{3u, 3u};

    /// Measurement smearing parameters
    static constexpr std::array<scalar, 2u> smearing{
        50.f * traccc::unit<scalar>::um, 50.f * traccc::unit<scalar>::um};

    /// Standard deviations for seed track parameters
    static constexpr std::array<scalar, e_bound_size> stddevs = {
        0.1f * traccc::unit<scalar>::mm,
        0.1f * traccc::unit<scalar>::mm,
        0.017f,
        0.017f,
        0.05f / traccc::unit<scalar>::GeV,
        1.f * traccc::unit<scalar>::ns};

    void consistency_tests(
        const edm::track_collection<default_algebra>::host::const_proxy_type&
            track,
        const edm::track_state_collection<default_algebra>::host&) const {

        // The nubmer of track states is supposed be greater than or
        // equal to the number of layers
        ASSERT_GE(track.constituent_links().size(), n_wire_layers);
    }

    protected:
    virtual void SetUp() override {
        vecmem::host_memory_resource host_mr;

        detray::wire_chamber_config<scalar> wire_chamber_cfg;
        wire_chamber_cfg.n_layers(n_wire_layers);
        wire_chamber_cfg.half_z(half_z);

        // Create telescope detector
        auto [det, name_map] = build_wire_chamber<traccc::default_algebra>(
            host_mr, wire_chamber_cfg);

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
