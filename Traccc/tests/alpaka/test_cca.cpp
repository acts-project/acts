/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "tests/cca_test.hpp"

// Project include(s).
#include "traccc/alpaka/clusterization/clusterization_algorithm.hpp"
#include "traccc/alpaka/utils/queue.hpp"
#include "traccc/alpaka/utils/vecmem_objects.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GoogleTest include(s).
#include <gtest/gtest.h>

// Standard include(s).
#include <functional>

namespace {

/// Long-lived host memory resource for the tests.
static vecmem::host_memory_resource host_mr;

// template <TAccTag>
cca_function_t get_f_with(traccc::clustering_config cfg) {
    return
        [cfg](const traccc::edm::silicon_cell_collection::host& cells,
              const traccc::detector_design_description::host& det_desc,
              const traccc::detector_conditions_description::host& det_cond)

            -> std::pair<
                std::map<traccc::geometry_id,
                         traccc::edm::measurement_collection::host>,
                std::optional<traccc::edm::silicon_cluster_collection::host>> {
            std::map<traccc::geometry_id,
                     traccc::edm::measurement_collection::host>
                result;

            traccc::alpaka::queue queue;
            traccc::alpaka::vecmem_objects vo(queue);

            vecmem::memory_resource& pinned_host_mr = vo.host_mr();
            vecmem::memory_resource& device_mr = vo.device_mr();
            vecmem::copy& copy = vo.async_copy();

            traccc::alpaka::clusterization_algorithm cc(
                {device_mr, &pinned_host_mr}, copy, queue, cfg);

            traccc::detector_design_description::buffer det_desc_buffer{
                [&]() {
                    std::vector<unsigned int> sizes(det_desc.size());
                    for (std::size_t i = 0; i < det_desc.size(); ++i) {
                        auto this_design = det_desc.at(i);

                        sizes[i] =
                            std::max(static_cast<unsigned int>(
                                         ((this_design.bin_edges_x()).size())),
                                     static_cast<unsigned int>(
                                         ((this_design.bin_edges_y()).size())));
                    }
                    return sizes;
                }(),
                device_mr, &host_mr, vecmem::data::buffer_type::resizable};
            copy.setup(det_desc_buffer)->wait();
            copy(vecmem::get_data(det_desc), det_desc_buffer)->wait();

            traccc::detector_conditions_description::buffer det_cond_buffer{
                static_cast<
                    traccc::detector_conditions_description::buffer::size_type>(
                    det_cond.size()),
                device_mr};
            copy.setup(det_cond_buffer)->ignore();
            copy(vecmem::get_data(det_cond), det_cond_buffer,
                 vecmem::copy::type::host_to_device)
                ->ignore();
            traccc::edm::silicon_cell_collection::buffer cells_buffer{
                static_cast<
                    traccc::edm::silicon_cell_collection::buffer::size_type>(
                    cells.size()),
                device_mr};
            copy.setup(cells_buffer)->wait();
            copy(vecmem::get_data(cells), cells_buffer)->wait();

            auto [measurements_buffer, cluster_buffer] =
                cc(cells_buffer, det_desc_buffer, det_cond_buffer,
                   traccc::device::clustering_keep_disjoint_set{});
            queue.synchronize();
            traccc::edm::measurement_collection::host measurements{
                pinned_host_mr};
            copy(measurements_buffer, measurements)->wait();

            traccc::edm::silicon_cluster_collection::host clusters{host_mr};
            copy(cluster_buffer, clusters)->wait();

            for (std::size_t i = 0; i < measurements.size(); i++) {
                if (result.contains(
                        measurements.at(i).surface_link().value()) == false) {
                    result.insert(
                        {measurements.at(i).surface_link().value(),
                         traccc::edm::measurement_collection::host{host_mr}});
                }
                result.at(measurements.at(i).surface_link().value())
                    .push_back(measurements.at(i));
            }

            return {result, clusters};
        };
}
}  // namespace

TEST_P(ConnectedComponentAnalysisTests, Run) {
    test_connected_component_analysis(GetParam());
}

INSTANTIATE_TEST_SUITE_P(
    AlpakaFastSvAlgorithm, ConnectedComponentAnalysisTests,
    ::testing::Combine(
        ::testing::Values(get_f_with(default_ccl_test_config())),
        ::testing::ValuesIn(ConnectedComponentAnalysisTests::get_test_files())),
    ConnectedComponentAnalysisTests::get_test_name);

INSTANTIATE_TEST_SUITE_P(
    AlpakaFastSvAlgorithmWithScratch, ConnectedComponentAnalysisTests,
    ::testing::Combine(
        ::testing::Values(get_f_with(tiny_ccl_test_config())),
        ::testing::ValuesIn(
            ConnectedComponentAnalysisTests::get_test_files_short())),
    ConnectedComponentAnalysisTests::get_test_name);
