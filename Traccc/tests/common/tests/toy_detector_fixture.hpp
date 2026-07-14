/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "test_detectors.hpp"
#include "traccc/io/detector.hpp"

// vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <string>

namespace traccc {

/// Test with Toy Geometry
class ToyDetectorFixture : public testing::Test {

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
        10.f * traccc::unit<scalar>::um, 25.f * traccc::unit<scalar>::um};

    // Grid search window
    static const inline std::array<detray::dindex, 2> search_window{3u, 3u};

    /// Standard deviations for seed track parameters
    static constexpr std::array<double, e_bound_size> stddevs = {
        smearing[0],
        smearing[1],
        0.5 * traccc::unit<double>::degree,
        0.5 * traccc::unit<double>::degree,
        0.01 / traccc::unit<double>::GeV,
        1000. * traccc::unit<double>::ns};

    void WriteDetector(const bool toggle_material_map,
                       const std::string& path = "./") {
        vecmem::host_memory_resource host_mr;

        detray::toy_det_config<scalar> toy_cfg{};
        toy_cfg.n_brl_layers(n_barrels)
            .n_edc_layers(n_endcaps)
            .envelope(2.f * traccc::unit<scalar>::mm)
            .use_material_maps(toggle_material_map)
            .do_check(false);

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
                              .path(path);
        detray::io::write_detector(det, name_map, writer_cfg);
    }
};

}  // namespace traccc
