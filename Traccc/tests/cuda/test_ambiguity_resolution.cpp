/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/ambiguity_resolution/ambiguity_resolution_config.hpp"
#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/cuda/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/device/container_d2h_copy_alg.hpp"
#include "traccc/device/container_h2d_copy_alg.hpp"
#include "traccc/utils/memory_resource.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <chrono>
#include <random>
#include <thread>

using namespace traccc;

void fill_measurements(edm::measurement_collection::host& measurements,
                       const measurement_id_type max_meas_id) {

    measurements.reserve(max_meas_id + 1);
    for (measurement_id_type i = 0; i <= max_meas_id; i++) {
        measurements.push_back({});
        measurements.at(measurements.size() - 1).identifier() = i;
    }
}

void fill_pattern(edm::track_container<default_algebra>::host& track_candidates,
                  const traccc::scalar pval,
                  const std::vector<measurement_id_type>& pattern) {

    track_candidates.tracks.resize(track_candidates.tracks.size() + 1u);
    track_candidates.tracks.pval().back() = pval;

    edm::measurement_collection::const_device measurements{
        track_candidates.measurements};

    for (const auto& meas_id : pattern) {
        const auto meas_iter =
            std::lower_bound(measurements.identifier().begin(),
                             measurements.identifier().end(), meas_id);

        const auto meas_idx =
            std::distance(measurements.identifier().begin(), meas_iter);
        track_candidates.tracks.constituent_links().back().push_back(
            {edm::track_constituent_link::measurement,
             static_cast<measurement_id_type>(meas_idx)});
    }
}

bool find_pattern(
    const edm::track_container<default_algebra>::const_device& tracks,
    const std::vector<measurement_id_type>& pattern) {

    const auto n_tracks = tracks.tracks.size();
    for (unsigned int i = 0; i < n_tracks; i++) {
        std::vector<measurement_id_type> ids;
        for (const auto& [type, meas_idx] :
             tracks.tracks.constituent_links().at(i)) {
            assert(type == edm::track_constituent_link::measurement);
            ids.push_back(tracks.measurements.at(meas_idx).identifier());
        }
        if (pattern == ids) {
            return true;
        }
    }
    return false;
}

std::vector<measurement_id_type> get_pattern(
    const edm::track_container<default_algebra>::host& track_candidates,
    const std::size_t idx) {

    edm::measurement_collection::const_device measurements{
        track_candidates.measurements};
    std::vector<measurement_id_type> ret;
    // A const reference would be fine here. But GCC fears that that would lead
    // to a dangling reference...
    const auto meas_links = track_candidates.tracks.at(idx).constituent_links();
    for (const auto& [type, meas_idx] : meas_links) {
        assert(type == edm::track_constituent_link::measurement);
        ret.push_back(measurements.at(meas_idx).identifier());
    }

    return ret;
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest0) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.23f, {5, 1, 11, 3});
    fill_pattern(trk_cands, 0.85f, {12, 10, 9, 8, 7, 6});
    fill_pattern(trk_cands, 0.42f, {4, 2, 13});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;

    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);
    {
        resolution_alg_cuda.get_config().min_meas_per_track = 3;
        auto res_trk_cands_buffer = resolution_alg_cuda(
            edm::track_container<default_algebra>::const_data(trk_cands));
        stream.synchronize();
        edm::track_container<default_algebra>::const_device res_trk_cands(
            res_trk_cands_buffer);
        // All tracks are accepted as they have more than three measurements
        EXPECT_EQ(res_trk_cands.tracks.size(), 3u);
        ASSERT_TRUE(find_pattern(res_trk_cands, {5, 1, 11, 3}));
        ASSERT_TRUE(find_pattern(res_trk_cands, {12, 10, 9, 8, 7, 6}));
        ASSERT_TRUE(find_pattern(res_trk_cands, {4, 2, 13}));
    }

    {
        resolution_alg_cuda.get_config().min_meas_per_track = 5;
        auto res_trk_cands_buffer = resolution_alg_cuda(
            edm::track_container<default_algebra>::const_data(trk_cands));
        stream.synchronize();
        edm::track_container<default_algebra>::const_device res_trk_cands(
            res_trk_cands_buffer);
        // Only the second track with six measurements is accepted
        ASSERT_EQ(res_trk_cands.tracks.size(), 1u);
        ASSERT_TRUE(find_pattern(res_trk_cands, {12, 10, 9, 8, 7, 6}));
    }
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest1) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.12f, {5, 14, 1, 11, 18, 16, 3});
    fill_pattern(trk_cands, 0.53f, {3, 6, 5, 13});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;

    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    resolution_alg_cuda.get_config().min_meas_per_track = 3;
    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    // All tracks are accepted as they have more than three measurements
    ASSERT_EQ(res_trk_cands.tracks.size(), 1u);

    // The first track is selected over the second one as its relative
    // shared measurement (2/7) is lower than the one of the second track
    // (2/4)
    ASSERT_TRUE(find_pattern(res_trk_cands, {5, 14, 1, 11, 18, 16, 3}));
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest2) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.8f, {1, 3, 5, 11});
    fill_pattern(trk_cands, 0.9f, {3, 5, 6, 13});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 1u);

    // The second track is selected over the first one as their relative
    // shared measurement (2/4) is the same but its p-value is higher
    ASSERT_TRUE(find_pattern(res_trk_cands, {3, 5, 6, 13}));
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest3) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.2f, {5, 1, 11, 3});
    fill_pattern(trk_cands, 0.5f, {6, 2});
    fill_pattern(trk_cands, 0.4f, {3, 21, 12, 6, 19, 14});
    fill_pattern(trk_cands, 0.1f, {13, 16, 2, 7, 11});
    fill_pattern(trk_cands, 0.3f, {1, 7, 8});
    fill_pattern(trk_cands, 0.6f, {1, 3, 11, 22});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 2u);

    ASSERT_TRUE(find_pattern(res_trk_cands, {3, 21, 12, 6, 19, 14}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {13, 16, 2, 7, 11}));
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest5) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.2f, {1, 2, 1, 1});
    fill_pattern(trk_cands, 0.5f, {3, 2, 1});
    fill_pattern(trk_cands, 0.4f, {2, 4, 5, 7, 2});
    fill_pattern(trk_cands, 0.1f, {6, 6, 6, 6});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 2u);

    ASSERT_TRUE(find_pattern(res_trk_cands, {3, 2, 1}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {6, 6, 6, 6}));
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest6) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.2f, {7, 3, 5, 7, 7, 7, 2});
    fill_pattern(trk_cands, 0.5f, {2});
    fill_pattern(trk_cands, 0.4f, {8, 9, 7, 2, 3, 4, 3, 7});
    fill_pattern(trk_cands, 0.1f, {8, 9, 0, 8, 1, 4, 6});
    fill_pattern(trk_cands, 0.9f, {10, 3, 2});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 2u);

    ASSERT_TRUE(find_pattern(res_trk_cands, {7, 3, 5, 7, 7, 7, 2}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {8, 9, 0, 8, 1, 4, 6}));
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest7) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.173853f, {10, 3, 6, 8});
    fill_pattern(trk_cands, 0.548019f, {3, 3, 1});
    fill_pattern(trk_cands, 0.276757f, {2, 8, 5, 4});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 1u);

    ASSERT_TRUE(find_pattern(res_trk_cands, {2, 8, 5, 4}));
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest8) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.0623132f, {10, 4});
    fill_pattern(trk_cands, 0.207417f, {6, 7, 5});
    fill_pattern(trk_cands, 0.325736f, {8, 2, 2});
    fill_pattern(trk_cands, 0.581643f, {5, 7, 9, 7});
    fill_pattern(trk_cands, 0.389551f, {1, 9, 3, 0});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 3u);

    ASSERT_TRUE(find_pattern(res_trk_cands, {6, 7, 5}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {8, 2, 2}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {1, 9, 3, 0}));
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest9) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.542984f, {0, 4, 8, 1, 1});
    fill_pattern(trk_cands, 0.583695f, {10, 6, 8, 7});
    fill_pattern(trk_cands, 0.280232f, {4, 1, 8, 10});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 1u);

    ASSERT_TRUE(find_pattern(res_trk_cands, {0, 4, 8, 1, 1}));
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest10) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.399106f, {14, 51});
    fill_pattern(trk_cands, 0.43899f, {80, 35, 41, 55});
    fill_pattern(trk_cands, 0.0954247f, {73, 63, 49, 89});
    fill_pattern(trk_cands, 0.158046f, {81, 22, 58, 54, 91});
    fill_pattern(trk_cands, 0.349878f, {97, 89, 80});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 3u);

    ASSERT_TRUE(find_pattern(res_trk_cands, {80, 35, 41, 55}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {73, 63, 49, 89}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {81, 22, 58, 54, 91}));
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest11) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.95f, {56, 87});
    fill_pattern(trk_cands, 0.894f, {64, 63});
    fill_pattern(trk_cands, 0.824f, {70, 17});
    fill_pattern(trk_cands, 0.862f, {27, 0});
    fill_pattern(trk_cands, 0.871f, {27, 19});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 0u);
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest12) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.948f, {17, 6, 1, 69, 78});  // 69
    fill_pattern(trk_cands, 0.609f, {17});
    fill_pattern(trk_cands, 0.453f, {84, 45, 81, 69});      // 84, 69
    fill_pattern(trk_cands, 0.910f, {54, 64, 49, 96, 40});  // 64
    fill_pattern(trk_cands, 0.153f, {59, 57, 84, 27, 64});  // 84, 64

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 2u);

    ASSERT_TRUE(find_pattern(res_trk_cands, {17, 6, 1, 69, 78}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {54, 64, 49, 96, 40}));
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest13) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.211f, {46, 92, 74, 58});
    fill_pattern(trk_cands, 0.694f, {15, 78, 9});
    fill_pattern(trk_cands, 0.432f, {15, 4, 58, 68});
    fill_pattern(trk_cands, 0.958f, {38, 93, 68});
    fill_pattern(trk_cands, 0.203f, {57, 64, 57, 36});
    fill_pattern(trk_cands, 0.118f, {4, 85, 65, 14});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 5u);

    ASSERT_TRUE(find_pattern(res_trk_cands, {46, 92, 74, 58}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {15, 78, 9}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {38, 93, 68}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {57, 64, 57, 36}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {4, 85, 65, 14}));
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest14) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.932f, {8, 4, 3});
    fill_pattern(trk_cands, 0.263f, {1, 1, 9, 3});
    fill_pattern(trk_cands, 0.876f, {1, 2, 5});
    fill_pattern(trk_cands, 0.058f, {2, 0, 4, 7});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 2u);

    ASSERT_TRUE(find_pattern(res_trk_cands, {8, 4, 3}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {1, 2, 5}));
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest15) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.293f, {2, 0, 4});
    fill_pattern(trk_cands, 0.362f, {8, 4, 9, 3});
    fill_pattern(trk_cands, 0.011f, {9, 4, 8, 4});
    fill_pattern(trk_cands, 0.843f, {8, 7, 1});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 2u);

    ASSERT_TRUE(find_pattern(res_trk_cands, {2, 0, 4}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {8, 7, 1}));
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest16) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.622598f, {95, 24, 62, 83, 67});
    fill_pattern(trk_cands, 0.541774f, {6, 52, 57, 87, 75});
    fill_pattern(trk_cands, 0.361033f, {14, 52, 29, 79, 89});
    fill_pattern(trk_cands, 0.622598f, {57, 85, 63, 90});
    fill_pattern(trk_cands, 0.481157f, {80, 45, 94});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 4u);

    ASSERT_TRUE(find_pattern(res_trk_cands, {95, 24, 62, 83, 67}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {14, 52, 29, 79, 89}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {57, 85, 63, 90}));
    ASSERT_TRUE(find_pattern(res_trk_cands, {80, 45, 94}));
}

TEST(CUDAAmbiguitySolverTests, GreedyResolverTest17) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    edm::measurement_collection::host measurements{mng_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        mng_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.17975f, {7, 4, 10, 3, 0});
    fill_pattern(trk_cands, 0.924326f, {0, 0, 9});
    fill_pattern(trk_cands, 0.0832954f, {0, 2, 0});
    fill_pattern(trk_cands, 0.303148f, {0, 6, 4, 5, 5});

    traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
        resolution_config, {mng_mr}, copy, stream);

    auto res_trk_cands_buffer = resolution_alg_cuda(
        edm::track_container<default_algebra>::const_data(trk_cands));
    stream.synchronize();
    edm::track_container<default_algebra>::const_device res_trk_cands(
        res_trk_cands_buffer);
    ASSERT_EQ(res_trk_cands.tracks.size(), 1u);

    ASSERT_TRUE(find_pattern(res_trk_cands, {0, 6, 4, 5, 5}));
}

// Test class for the ambiguity resolution comparison with CPU implementation
// Input tuple: < n_event, n_tracks, track_length_range , max_meas_id,
// allow_duplicate >
class CUDAGreedyResolutionCompareToCPU
    : public ::testing::TestWithParam<
          std::tuple<std::size_t, std::size_t, std::array<std::size_t, 2u>,
                     measurement_id_type, bool>> {};

TEST_P(CUDAGreedyResolutionCompareToCPU, Comparison) {

    const std::size_t n_events = std::get<0>(GetParam());
    const std::size_t n_tracks = std::get<1>(GetParam());
    const std::array<std::size_t, 2u> trk_length_range =
        std::get<2>(GetParam());
    const measurement_id_type max_meas_id = std::get<3>(GetParam());
    const bool allow_duplicate = std::get<4>(GetParam());

    // Memory resource used by the EDM.
    vecmem::cuda::device_memory_resource device_mr;
    vecmem::host_memory_resource host_mr;
    traccc::memory_resource mr{device_mr, &host_mr};

    // Cuda stream
    traccc::cuda::stream stream;

    // Cuda copy objects
    vecmem::cuda::async_copy copy{stream.cudaStream()};

    for (std::size_t i_evt = 0u; i_evt < n_events; i_evt++) {

        std::size_t sd = 42u + i_evt;
        std::mt19937 gen(sd);
        std::cout << "Event: " << i_evt << " Seed: " << sd << std::endl;

        edm::measurement_collection::host measurements{host_mr};
        fill_measurements(measurements, max_meas_id);
        edm::track_container<default_algebra>::host trk_cands{
            host_mr, vecmem::get_data(measurements)};

        for (std::size_t i = 0; i < n_tracks; i++) {

            std::uniform_int_distribution<std::size_t> track_length_dist(
                trk_length_range[0], trk_length_range[1]);
            std::uniform_int_distribution<measurement_id_type> meas_id_dist(
                0, max_meas_id);
            std::uniform_real_distribution<traccc::scalar> pval_dist(0.0f,
                                                                     1.0f);

            const std::size_t track_length = track_length_dist(gen);
            const traccc::scalar pval = pval_dist(gen);
            std::vector<measurement_id_type> pattern;
            // std::cout << pval << std::endl;
            while (pattern.size() < track_length) {

                auto mid = meas_id_dist(gen);
                if (!allow_duplicate) {
                    while (std::find(pattern.begin(), pattern.end(), mid) !=
                           pattern.end()) {

                        mid = meas_id_dist(gen);
                    }
                }
                // std::cout << mid << ", ";
                pattern.push_back(mid);
            }
            // std::cout << std::endl;

            // Make sure that partern size is eqaul to the track length
            ASSERT_EQ(pattern.size(), track_length);

            // Fill the pattern
            fill_pattern(trk_cands, pval, pattern);
        }

        // CPU algorithm
        traccc::host::greedy_ambiguity_resolution_algorithm::config_type
            resolution_config;
        traccc::host::greedy_ambiguity_resolution_algorithm resolution_alg_cpu(
            resolution_config, host_mr);

        auto start_cpu = std::chrono::high_resolution_clock::now();

        auto res_trk_cands_cpu = resolution_alg_cpu(
            edm::track_container<default_algebra>::const_data(trk_cands));

        auto end_cpu = std::chrono::high_resolution_clock::now();
        auto duration_cpu =
            std::chrono::duration_cast<std::chrono::milliseconds>(end_cpu -
                                                                  start_cpu);
        std::cout << " Time for the cpu method " << duration_cpu.count()
                  << " ms" << std::endl;

        // CUDA algorithm
        traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
            resolution_config, mr, copy, stream);

        // H2D transfer
        edm::measurement_collection::buffer measurements_buffer =
            copy.to(vecmem::get_data(measurements), device_mr, &host_mr,
                    vecmem::copy::type::host_to_device);
        traccc::edm::track_container<default_algebra>::buffer trk_cands_buffer{
            copy.to(vecmem::get_data(trk_cands.tracks), device_mr, &host_mr,
                    vecmem::copy::type::host_to_device),
            {},
            measurements_buffer};

        auto start_cuda = std::chrono::high_resolution_clock::now();

        // Instantiate output cuda containers/collections
        auto res_trk_cands_buffer = resolution_alg_cuda(trk_cands_buffer);
        stream.synchronize();

        auto end_cuda = std::chrono::high_resolution_clock::now();
        auto duration_cuda =
            std::chrono::duration_cast<std::chrono::milliseconds>(end_cuda -
                                                                  start_cuda);
        std::cout << " Time for the cuda method " << duration_cuda.count()
                  << " ms" << std::endl;

        traccc::edm::track_container<default_algebra>::buffer
            res_trk_cands_cuda{
                copy.to(res_trk_cands_buffer.tracks, host_mr, nullptr,
                        vecmem::copy::type::device_to_host),
                {},
                vecmem::get_data(measurements)};

        const auto n_tracks_cpu = res_trk_cands_cpu.tracks.size();
        ASSERT_EQ(n_tracks_cpu, res_trk_cands_cuda.tracks.capacity());

        // Make sure that CPU and CUDA track candidates have same
        // patterns
        edm::track_container<default_algebra>::const_device
            res_trk_cands_cuda_device{res_trk_cands_cuda};
        for (unsigned int i = 0; i < n_tracks_cpu; i++) {
            ASSERT_TRUE(find_pattern(res_trk_cands_cuda_device,
                                     get_pattern(res_trk_cands_cpu, i)));
        }
    }
};

INSTANTIATE_TEST_SUITE_P(
    CUDAStandard, CUDAGreedyResolutionCompareToCPU,
    ::testing::Values(std::make_tuple(5u, 50000u,
                                      std::array<std::size_t, 2u>{1u, 10u},
                                      20000u, true),
                      std::make_tuple(5u, 50000u,
                                      std::array<std::size_t, 2u>{1u, 10u},
                                      20000u, false)));

INSTANTIATE_TEST_SUITE_P(
    CUDASparse, CUDAGreedyResolutionCompareToCPU,
    ::testing::Values(std::make_tuple(3u, 5000u,
                                      std::array<std::size_t, 2u>{3u, 10u},
                                      1000000u, true),
                      std::make_tuple(3u, 5000u,
                                      std::array<std::size_t, 2u>{3u, 10u},
                                      1000000u, false)));

INSTANTIATE_TEST_SUITE_P(
    CUDADense, CUDAGreedyResolutionCompareToCPU,
    ::testing::Values(std::make_tuple(3u, 5000u,
                                      std::array<std::size_t, 2u>{3u, 10u},
                                      100u, true),
                      std::make_tuple(3u, 5000u,
                                      std::array<std::size_t, 2u>{3u, 10u},
                                      100u, false)));

INSTANTIATE_TEST_SUITE_P(
    CUDALong, CUDAGreedyResolutionCompareToCPU,
    ::testing::Values(std::make_tuple(3u, 10000u,
                                      std::array<std::size_t, 2u>{3u, 500u},
                                      10000u, true),
                      std::make_tuple(3u, 10000u,
                                      std::array<std::size_t, 2u>{3u, 500u},
                                      10000u, false)));

INSTANTIATE_TEST_SUITE_P(
    CUDASimple, CUDAGreedyResolutionCompareToCPU,
    ::testing::Values(
        std::make_tuple(3u, 5u, std::array<std::size_t, 2u>{3u, 5u}, 10u, true),
        std::make_tuple(3u, 5u, std::array<std::size_t, 2u>{3u, 5u}, 10u,
                        false)));
