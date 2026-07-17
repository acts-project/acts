/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/clusterization/measurement_creation_algorithm.hpp"
#include "traccc/clusterization/sparse_ccl_algorithm.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/edm/silicon_cluster_collection.hpp"
#include "traccc/geometry/detector_design_description.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

TEST(algorithms, seq_single_module) {

    // Memory resource used in the test.
    vecmem::host_memory_resource resource;

    traccc::host::sparse_ccl_algorithm cc(resource);
    traccc::host::measurement_creation_algorithm mc(resource);

    /// Following [DOI: 10.1109/DASIP48288.2019.9049184]
    traccc::edm::silicon_cell_collection::host cells{resource};
    cells.push_back({1, 0, 1.f, 0.f, 0});
    cells.push_back({8, 4, 2.f, 0.f, 0});
    cells.push_back({10, 4, 3.f, 0.f, 0});
    cells.push_back({9, 5, 4.f, 0.f, 0});
    cells.push_back({10, 5, 5.f, 0.f, 0});
    cells.push_back({12, 12, 6.f, 0.f, 0});
    cells.push_back({3, 13, 7.f, 0.f, 0});
    cells.push_back({11, 13, 8.f, 0.f, 0});
    cells.push_back({4, 14, 9.f, 0.f, 0});

    // Create a dummy detector description. With a description of enough
    // detector modules for all the input files that the test uses.
    static constexpr std::size_t NMODULES = 1;
    traccc::detector_design_description::host det_desc{resource};
    traccc::detector_conditions_description::host det_cond{resource};
    det_desc.resize(NMODULES);
    det_cond.resize(NMODULES);
    for (std::size_t i = 0; i < NMODULES; ++i) {
        det_cond.module_to_design_id()[i] = static_cast<unsigned int>(i);
        det_cond.geometry_id()[i] = detray::geometry::identifier{i};
        det_cond.acts_geometry_id()[i] = i;
        det_cond.measurement_translation()[i] = {0.f, 0.f};

        std::vector<traccc::scalar,
                    std::pmr::polymorphic_allocator<traccc::scalar>>
            bin_edges_x(10001), bin_edges_y(10001);
        std::iota(bin_edges_x.begin(), bin_edges_x.end(), -0.5f);
        std::iota(bin_edges_y.begin(), bin_edges_y.end(), -0.5f);
        det_desc.bin_edges_x().back().assign(bin_edges_x.begin(),
                                             bin_edges_x.end());

        det_desc.bin_edges_y().back().assign(bin_edges_y.begin(),
                                             bin_edges_y.end());
        det_desc.dimensions()[i] = 2;
        det_desc.subspace()[i] = {0, 1};
        det_desc.design_id()[i] = static_cast<int>(i);
    }

    auto cells_data = vecmem::get_data(cells);
    auto clusters = cc(cells_data, vecmem::get_data(det_cond));
    EXPECT_EQ(clusters.size(), 4u);

    auto clusters_data = vecmem::get_data(clusters);
    auto det_desc_data = vecmem::get_data(det_desc);
    auto det_cond_data = vecmem::get_data(det_cond);
    auto measurements =
        mc(cells_data, clusters_data, det_desc_data, det_cond_data);

    EXPECT_EQ(measurements.size(), 4u);
}
