/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/geometry/detector.hpp"
#include "traccc/io/read_cells.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/read_digitization_config.hpp"
#include "traccc/io/read_measurements.hpp"
#include "traccc/io/read_particles.hpp"
#include "traccc/io/read_spacepoints.hpp"
#include "traccc/io/write.hpp"

// Test include(s).
#include "tests/data_test.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <filesystem>

class io : public traccc::tests::data_test {};

// This defines the local frame test suite
TEST_F(io, csv_read_single_module) {

    vecmem::host_memory_resource resource;

    traccc::edm::silicon_cell_collection::host cells{resource};
    traccc::io::read_cells(cells, get_datafile("single_module/cells.csv"),
                           traccc::getDummyLogger().clone());

    ASSERT_EQ(cells.size(), 6u);

    EXPECT_EQ(cells.channel0().at(0), 123u);
    EXPECT_EQ(cells.channel1().at(0), 32u);
    EXPECT_FLOAT_EQ(static_cast<float>(cells.activation().at(0)), 0.002f);
    EXPECT_EQ(cells.module_index().at(0), 0u);

    EXPECT_EQ(cells.channel0().at(5), 174u);
    EXPECT_EQ(cells.channel1().at(5), 880u);
    EXPECT_FLOAT_EQ(static_cast<float>(cells.activation().at(5)), 0.23f);
    EXPECT_EQ(cells.module_index().at(5), 0u);
}

// This defines the local frame test suite
TEST_F(io, csv_read_two_modules) {

    vecmem::host_memory_resource resource;

    traccc::edm::silicon_cell_collection::host cells{resource};
    traccc::io::read_cells(cells, get_datafile("two_modules/cells.csv"),
                           traccc::getDummyLogger().clone());

    ASSERT_EQ(cells.size(), 14u);

    // Check cells in first module
    EXPECT_EQ(cells.channel0().at(0), 123u);
    EXPECT_EQ(cells.channel1().at(0), 32u);
    EXPECT_EQ(cells.channel0().at(5), 174u);
    EXPECT_EQ(cells.channel1().at(5), 880u);

    // Check cells in second module
    EXPECT_EQ(cells.channel0().at(6), 0u);
    EXPECT_EQ(cells.channel1().at(6), 4u);
    EXPECT_EQ(cells.channel0().at(13), 5u);
    EXPECT_EQ(cells.channel1().at(13), 98u);
}

// This reads in the tml pixel barrel first event
TEST_F(io, csv_read_tml_pixelbarrel) {

    vecmem::host_memory_resource resource;

    traccc::edm::silicon_cell_collection::host cells{resource};
    traccc::io::read_cells(
        cells, get_datafile("tml_pixel_barrel/event000000000-cells.csv"),
        traccc::getDummyLogger().clone(), nullptr, traccc::data_format::csv,
        false);

    EXPECT_EQ(cells.size(), 179961u);
}

/// Tests with ODD "single" muon events.
TEST_F(io, csv_read_odd_single_muon) {

    // Memory resource used by the test.
    vecmem::host_memory_resource mr;

    traccc::host_detector detector;
    traccc::io::read_detector(detector, mr,
                              "geometries/odd/odd-detray_geometry_detray.json");

    // Read the truth particles for the first event.
    traccc::particle_container_types::host particles{&mr};
    traccc::edm::measurement_collection::host measurements{mr};
    traccc::io::read_particles(particles, measurements, 0u,
                               "odd/geant4_1muon_1GeV/", &detector,
                               traccc::data_format::csv);

    // Look at the read container.
    ASSERT_EQ(particles.size(), 265u);
    std::size_t n_muons = 0u;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        // The muon(s) must have measurements associated to it/them.
        if (std::abs(particles.at(i).header.particle_type) == 13) {
            ++n_muons;
            EXPECT_GT(particles.at(i).items.size(), 0u);
        }
    }
    EXPECT_EQ(n_muons, 4u);
}

TEST_F(io, csv_write_odd_single_muon_cells) {

    // Memory resource used by the test.
    vecmem::host_memory_resource mr;

    // Read the ODD detector description.
    traccc::detector_design_description::host det_desc{mr};
    traccc::detector_conditions_description::host det_cond{mr};
    traccc::io::read_detector_description(
        det_desc, det_cond, "geometries/odd/odd-detray_geometry_detray.json",
        "geometries/odd/odd-digi-geometric-config.json",
        "geometries/odd/odd-digi-geometric-config.json");

    // Lambda comparing two cell collections.
    auto compare_cells =
        [](const traccc::edm::silicon_cell_collection::host& a,
           const traccc::edm::silicon_cell_collection::host& b) -> void {
        ASSERT_EQ(a.size(), b.size());
        for (traccc::edm::silicon_cell_collection::host::size_type i = 0;
             i < a.size(); ++i) {
            EXPECT_EQ(a.at(i), b.at(i));
        }
    };

    // Cell collections to use in the test.
    traccc::edm::silicon_cell_collection::host orig{mr}, copy{mr};

    // Lambda performing the test with either using "Acts geometry IDs" or
    // "Detray ones".
    auto perform_test = [&](bool use_acts_geometry_id) {
        // Test the I/O for 10 events.
        for (std::size_t event = 0; event < 10; ++event) {

            // Read the cells for the current event.
            traccc::io::read_cells(orig, event, "odd/geant4_1muon_1GeV/",
                                   traccc::getDummyLogger().clone(), &det_cond);
            // Write the cells into a temporary file.
            traccc::io::write(event,
                              std::filesystem::temp_directory_path().native(),
                              traccc::data_format::csv, vecmem::get_data(orig),
                              vecmem::get_data(det_desc),
                              vecmem::get_data(det_cond), use_acts_geometry_id);

            // Read the cells back in.
            traccc::io::read_cells(
                copy, event, std::filesystem::temp_directory_path().native(),
                traccc::getDummyLogger().clone(), &det_cond,
                traccc::data_format::csv, false, use_acts_geometry_id);

            // Compare the two cell collections.
            compare_cells(orig, copy);
        }
    };

    // Perform the test with both "Acts geometry IDs" and "Detray ones".
    perform_test(true);
    perform_test(false);
}
