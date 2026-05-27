/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "test_detectors.hpp"
#include "triplet_fitting_test.hpp"

// Detray include(s).
#include <detray/io/frontend/detector_writer.hpp>

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

namespace traccc {

/// Kalman Fitting Test with Telescope Detector
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
class TripletFittingTelescopeTests
    : public TripletFittingTests,
      public testing::WithParamInterface<std::tuple<
          std::string, std::array<scalar, 3u>, std::array<scalar, 3u>,
          std::array<scalar, 2u>, std::array<scalar, 2u>,
          std::array<scalar, 2u>, detray::pdg_particle<scalar>, unsigned int,
          unsigned int, bool, scalar, unsigned int, scalar, vector3>> {

    public:
    /// Plane alignment direction (aligned to x-axis)
    static const inline detray::detail::ray<traccc::default_algebra> traj{
        {0, 0, 0}, 0, {1, 0, 0}, -1};

    /// Plane material and thickness
    static const inline detray::silicon_tml<scalar> mat = {};
    static constexpr scalar thickness = 0.5f * traccc::unit<scalar>::mm;

    // Rectangle mask for the telescope geometry
    static constexpr detray::mask<detray::rectangle2D, traccc::default_algebra>
        rectangle{0u, 100000.f, 100000.f};

    /// Measurement smearing parameters
    static constexpr std::array<scalar, 2u> smearing{
        50 * traccc::unit<scalar>::um, 50 * traccc::unit<scalar>::um};

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

        // The nubmer of track states is supposed be equal to the number
        // of planes
        ASSERT_EQ(track.constituent_links().size(), std::get<11>(GetParam()));
    }

    protected:
    virtual void SetUp() override {

        vecmem::host_memory_resource host_mr;

        const scalar offset = std::get<10>(GetParam());
        const unsigned int n_planes = std::get<11>(GetParam());
        const scalar spacing = std::get<12>(GetParam());

        std::vector<scalar> plane_positions;
        for (unsigned int i = 0; i < n_planes; i++) {
            plane_positions.push_back(offset * unit<scalar>::mm +
                                      static_cast<scalar>(i) * spacing *
                                          unit<scalar>::mm);
        }

        detray::tel_det_config tel_cfg{rectangle};
        tel_cfg.positions(plane_positions);
        tel_cfg.module_material(mat);
        tel_cfg.mat_thickness(thickness);
        tel_cfg.pilot_track(traj);
        tel_cfg.envelope(offset * 2.f);

        // Create telescope detector
        auto [det, name_map] = build_telescope_detector(host_mr, tel_cfg);

        // Write detector file
        auto writer_cfg = detray::io::detector_writer_config{}
                              .format(detray::io::format::json)
                              .replace_files(true)
                              .write_material(true)
                              .path(std::get<0>(GetParam()));
        detray::io::write_detector(det, name_map, writer_cfg);
    }
};

}  // namespace traccc
