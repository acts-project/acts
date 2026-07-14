/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/clusterization/clustering_config.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/edm/silicon_cluster_collection.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/io/csv/dfe.hpp"
#include "traccc/io/csv/make_cell_reader.hpp"
#include "traccc/io/read_cells.hpp"

// Test include(s).
#include "tests/data_test.hpp"

// VecMem include(s).
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <functional>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

using cca_function_t = std::function<std::pair<
    std::map<traccc::geometry_id, traccc::edm::measurement_collection::host>,
    std::optional<traccc::edm::silicon_cluster_collection::host>>(
    const traccc::edm::silicon_cell_collection::host &,
    const traccc::detector_design_description::host &,
    const traccc::detector_conditions_description::host &)>;

inline traccc::clustering_config default_ccl_test_config() {
    traccc::clustering_config rv;

    rv.threads_per_partition = 128;
    rv.max_cells_per_thread = 16;
    rv.target_cells_per_thread = 8;
    rv.backup_size_multiplier = 256;

    return rv;
}

inline traccc::clustering_config tiny_ccl_test_config() {
    traccc::clustering_config rv;

    rv.threads_per_partition = 128;
    rv.max_cells_per_thread = 1;
    rv.target_cells_per_thread = 1;
    rv.backup_size_multiplier = 16384;

    return rv;
}

class ConnectedComponentAnalysisTests
    : public traccc::tests::data_test,
      public testing::WithParamInterface<
          std::tuple<cca_function_t, std::string>> {
    public:
    struct cca_truth_hit {
        uint64_t geometry_id = 0;
        traccc::measurement_id_type measurement_id = 0;
        uint64_t num_cells = 0;
        traccc::scalar channel0 = 0;
        traccc::scalar channel1 = 0;
        traccc::scalar variance0 = 0.f;
        traccc::scalar variance1 = 0.f;

        DFE_NAMEDTUPLE(cca_truth_hit, geometry_id, measurement_id, num_cells,
                       channel0, channel1, variance0, variance1);
    };

    using cca_truth_hit_reader =
        ::traccc::io::dfe::NamedTupleCsvReader<cca_truth_hit>;

    inline static std::string get_test_name(
        const testing::TestParamInfo<ParamType> &info) {
        return std::get<1>(info.param);
    }

    inline static std::vector<std::string> get_test_files(void) {
        const std::vector<std::pair<std::string, std::size_t>> cases = {
            {"dense", 100},
            {"multiple_module_single_hit", 100},
            {"single_module_multiple_hit_single", 100},
            {"single_module_multiple_hit_single_sparse", 100},
            {"single_module_single_hit", 100},
            {"very_dense", 100},
            {"trackml_like", 30},
        };
        std::vector<std::string> out;

        for (const std::pair<std::string, std::size_t> &c : cases) {
            for (std::size_t i = 0; i < c.second; ++i) {
                std::ostringstream ss;
                ss << c.first << "_" << std::setfill('0') << std::setw(10) << i;
                out.push_back(ss.str());
            }
        }

        return out;
    }

    inline static std::vector<std::string> get_test_files_short(void) {
        const std::vector<std::pair<std::string, std::size_t>> cases = {
            {"trackml_like", 10},
        };
        std::vector<std::string> out;

        for (const std::pair<std::string, std::size_t> &c : cases) {
            for (std::size_t i = 0; i < c.second; ++i) {
                std::ostringstream ss;
                ss << c.first << "_" << std::setfill('0') << std::setw(10) << i;
                out.push_back(ss.str());
            }
        }

        return out;
    }

    inline void test_connected_component_analysis(ParamType p) {
        cca_function_t f = std::get<0>(p);
        std::string file_prefix = std::get<1>(p);

        std::string file_hits =
            get_datafile("cca_test/" + file_prefix + "_hits.csv");
        std::string file_truth =
            get_datafile("cca_test/" + file_prefix + "_truth.csv");

        // Host memory resource for the test.
        vecmem::host_memory_resource mr;

        // Create a dummy detector description. With a description of enough
        // detector modules for all the input files that the test uses.
        static constexpr std::size_t NMODULES = 2500;
        static constexpr traccc::scalar pitch = 1.f;
        traccc::detector_design_description::host det_desc{mr};
        traccc::detector_conditions_description::host det_cond{mr};
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
            std::iota(bin_edges_x.begin(), bin_edges_x.end(), -0.5);
            std::iota(bin_edges_y.begin(), bin_edges_y.end(), -0.5);
            (det_desc.bin_edges_x()[i].assign(bin_edges_x.begin(),
                                              bin_edges_x.end()));
            (det_desc.bin_edges_y()[i])
                .assign(bin_edges_y.begin(), bin_edges_y.end());
            det_desc.dimensions()[i] = 2;
            det_desc.subspace()[i] = {0, 1};
            det_desc.design_id()[i] = static_cast<int>(i);
        }

        ASSERT_EQ(det_desc.size(), det_cond.size()) << "not of same size";

        traccc::edm::silicon_cell_collection::host cells{mr};
        traccc::io::read_cells(cells, file_hits,
                               traccc::getDummyLogger().clone(), &det_cond);

        auto [result, cluster_data] = f(cells, det_desc, det_cond);

        std::size_t total_truth = 0, total_found = 0;

        for (const auto &i : result) {
            total_found += i.second.size();
        }

        cca_truth_hit_reader truth_reader(file_truth);

        traccc::scalar var_adjustment = (pitch * pitch) / 12.f;

        cca_truth_hit io_truth;
        while (truth_reader.read(io_truth)) {
            ASSERT_TRUE(result.find(io_truth.geometry_id) != result.end());

            const traccc::edm::measurement_collection::host &meas =
                result.at(io_truth.geometry_id);

            const traccc::scalar tol = 0.0001f;

            std::size_t meas_idx = static_cast<std::size_t>(-1);
            for (std::size_t i = 0; i < meas.size(); ++i) {

                if ((std::abs(meas.at(i).local_position()[0] -
                              io_truth.channel0) < tol) &&
                    (std::abs(meas.at(i).local_position()[1] -
                              io_truth.channel1) < tol)) {
                    meas_idx = i;
                    break;
                }
            }
            ASSERT_TRUE(meas_idx < meas.size())
                << "measurement not found " << meas.size() << " " << meas_idx;

            const auto match = meas.at(meas_idx);
            EXPECT_NEAR(match.local_position()[0], io_truth.channel0, tol);
            EXPECT_NEAR(match.local_position()[1], io_truth.channel1, tol);
            EXPECT_NEAR(match.local_variance()[0],
                        io_truth.variance0 + var_adjustment, tol);
            EXPECT_NEAR(match.local_variance()[1],
                        io_truth.variance1 + var_adjustment, tol);

            ++total_truth;
        }

        EXPECT_EQ(total_truth, total_found);

        if (cluster_data.has_value()) {
            ASSERT_EQ(cluster_data->size(), total_found);

            std::size_t total_cluster_size = 0;

            for (std::size_t i = 0; i < cluster_data->size(); ++i) {
                total_cluster_size += cluster_data->cell_indices().at(i).size();
            }

            ASSERT_EQ(cells.size(), total_cluster_size);
        }
    }
};
