/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"

// Test include(s).
#include "tests/cca_test.hpp"

// VecMem include(s).
#include <optional>
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <functional>

namespace {
vecmem::host_memory_resource resource;
traccc::host::sparse_ccl_algorithm cc(resource);
traccc::host::measurement_creation_algorithm mc(resource);

cca_function_t f =
    [](const traccc::edm::silicon_cell_collection::host& cells,
       const traccc::detector_design_description::host& det_desc,
       const traccc::detector_conditions_description::host& det_cond)
    -> std::pair<std::map<traccc::geometry_id,
                          traccc::edm::measurement_collection::host>,
                 std::optional<traccc::edm::silicon_cluster_collection::host>> {
    std::map<traccc::geometry_id, traccc::edm::measurement_collection::host>
        result;

    const traccc::edm::silicon_cell_collection::const_data cells_data =
        vecmem::get_data(cells);
    const traccc::detector_design_description::const_data det_desc_data =
        vecmem::get_data(det_desc);
    const traccc::detector_conditions_description::const_data det_cond_data =
        vecmem::get_data(det_cond);
    const auto clusters = cc(cells_data, det_cond_data);
    const auto clusters_data = vecmem::get_data(clusters);
    auto measurements =
        mc(cells_data, clusters_data, det_desc_data, det_cond_data);

    for (std::size_t i = 0; i < measurements.size(); i++) {
        if (result.contains(measurements.at(i).surface_link().value()) ==
            false) {
            result.insert(
                {measurements.at(i).surface_link().value(),
                 traccc::edm::measurement_collection::host{resource}});
        }
        result.at(measurements.at(i).surface_link().value())
            .push_back(measurements.at(i));
    }

    return {std::move(result), std::nullopt};
};
}  // namespace

TEST_P(ConnectedComponentAnalysisTests, Run) {
    test_connected_component_analysis(GetParam());
}

INSTANTIATE_TEST_SUITE_P(
    SparseCclAlgorithm, ConnectedComponentAnalysisTests,
    ::testing::Combine(
        ::testing::Values(f),
        ::testing::ValuesIn(ConnectedComponentAnalysisTests::get_test_files())),
    ConnectedComponentAnalysisTests::get_test_name);
