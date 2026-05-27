/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/ambiguity_resolution/legacy/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/utils/memory_resource.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

#include <chrono>
#include <random>

using namespace traccc;

namespace {
vecmem::host_memory_resource host_mr;
}  // namespace

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

std::vector<std::size_t> get_pattern(
    const edm::track_container<default_algebra>::host& track_candidates,
    const std::size_t idx) {

    edm::measurement_collection::const_device measurements{
        track_candidates.measurements};
    std::vector<std::size_t> ret;
    // A const reference would be fine here. But GCC fears that that would lead
    // to a dangling reference...
    const auto meas_links = track_candidates.tracks.at(idx).constituent_links();
    for (const auto& [type, meas_idx] : meas_links) {
        assert(type == edm::track_constituent_link::measurement);
        ret.push_back(measurements.at(meas_idx).identifier());
    }

    return ret;
}

TEST(AmbiguitySolverTests, GreedyResolverTest0) {

    edm::measurement_collection::host measurements{host_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        host_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.23f, {5, 1, 11, 3});
    fill_pattern(trk_cands, 0.85f, {12, 10, 9, 8, 7, 6});
    fill_pattern(trk_cands, 0.42f, {4, 2, 13});

    traccc::host::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;

    traccc::host::greedy_ambiguity_resolution_algorithm resolution_alg(
        resolution_config, host_mr);
    {
        resolution_alg.get_config().min_meas_per_track = 3;
        auto res_trk_cands = resolution_alg(
            traccc::edm::track_container<traccc::default_algebra>::const_data(
                trk_cands));
        // All tracks are accepted as they have more than three measurements
        ASSERT_EQ(res_trk_cands.tracks.size(), 3u);
        ASSERT_EQ(get_pattern(res_trk_cands, 0),
                  std::vector<std::size_t>({5, 1, 11, 3}));
        ASSERT_EQ(get_pattern(res_trk_cands, 1),
                  std::vector<std::size_t>({12, 10, 9, 8, 7, 6}));
        ASSERT_EQ(get_pattern(res_trk_cands, 2),
                  std::vector<std::size_t>({4, 2, 13}));
    }

    {
        resolution_alg.get_config().min_meas_per_track = 5;
        auto res_trk_cands = resolution_alg(
            traccc::edm::track_container<traccc::default_algebra>::const_data(
                trk_cands));
        // Only the second track with six measurements is accepted
        ASSERT_EQ(res_trk_cands.tracks.size(), 1u);
        ASSERT_EQ(get_pattern(res_trk_cands, 0),
                  std::vector<std::size_t>({12, 10, 9, 8, 7, 6}));
    }

    /*******************
     * Legacy algorithm
     * *****************/

    {
        traccc::legacy::greedy_ambiguity_resolution_algorithm::config_t
            legacy_cfg;
        traccc::legacy::greedy_ambiguity_resolution_algorithm
            legacy_resolution_alg(legacy_cfg, host_mr);

        legacy_resolution_alg.get_config().n_measurements_min = 3;
        auto res_trk_cands = legacy_resolution_alg(trk_cands);
        ASSERT_EQ(res_trk_cands.tracks.size(), 3u);
        ASSERT_EQ(get_pattern(res_trk_cands, 0),
                  std::vector<std::size_t>({5, 1, 11, 3}));
        ASSERT_EQ(get_pattern(res_trk_cands, 1),
                  std::vector<std::size_t>({12, 10, 9, 8, 7, 6}));
        ASSERT_EQ(get_pattern(res_trk_cands, 2),
                  std::vector<std::size_t>({4, 2, 13}));
    }

    {
        traccc::legacy::greedy_ambiguity_resolution_algorithm::config_t
            legacy_cfg;
        traccc::legacy::greedy_ambiguity_resolution_algorithm
            legacy_resolution_alg(legacy_cfg, host_mr);

        legacy_resolution_alg.get_config().n_measurements_min = 5;
        auto res_trk_cands = legacy_resolution_alg(trk_cands);
        ASSERT_EQ(res_trk_cands.tracks.size(), 1u);
        ASSERT_EQ(get_pattern(res_trk_cands, 0),
                  std::vector<std::size_t>({12, 10, 9, 8, 7, 6}));
    }
}

TEST(AmbiguitySolverTests, GreedyResolverTest1) {

    edm::measurement_collection::host measurements{host_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        host_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.12f, {5, 14, 1, 11, 18, 16, 3});
    fill_pattern(trk_cands, 0.53f, {3, 6, 5, 13});

    traccc::host::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::host::greedy_ambiguity_resolution_algorithm resolution_alg(
        resolution_config, host_mr);
    {
        auto res_trk_cands = resolution_alg(
            traccc::edm::track_container<traccc::default_algebra>::const_data(
                trk_cands));
        ASSERT_EQ(res_trk_cands.tracks.size(), 1u);

        // The first track is selected over the second one as its relative
        // shared measurement (2/7) is lower than the one of the second track
        // (2/4)
        ASSERT_EQ(get_pattern(res_trk_cands, 0),
                  std::vector<std::size_t>({5, 14, 1, 11, 18, 16, 3}));
    }

    /*******************
     * Legacy algorithm
     * *****************/

    {
        traccc::legacy::greedy_ambiguity_resolution_algorithm::config_t
            legacy_cfg;
        traccc::legacy::greedy_ambiguity_resolution_algorithm
            legacy_resolution_alg(legacy_cfg, host_mr);

        auto res_trk_cands = legacy_resolution_alg(trk_cands);
        ASSERT_EQ(res_trk_cands.tracks.size(), 1u);

        // The first track is selected over the second one as its relative
        // shared measurement (2/7) is lower than the one of the second track
        // (2/4)
        ASSERT_EQ(get_pattern(res_trk_cands, 0),
                  std::vector<std::size_t>({5, 14, 1, 11, 18, 16, 3}));
    }
}

TEST(AmbiguitySolverTests, GreedyResolverTest2) {

    edm::measurement_collection::host measurements{host_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        host_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.8f, {1, 3, 5, 11});
    fill_pattern(trk_cands, 0.9f, {3, 5, 6, 13});

    traccc::host::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::host::greedy_ambiguity_resolution_algorithm resolution_alg(
        resolution_config, host_mr);
    {
        auto res_trk_cands = resolution_alg(
            traccc::edm::track_container<traccc::default_algebra>::const_data(
                trk_cands));
        ASSERT_EQ(res_trk_cands.tracks.size(), 1u);

        // The second track is selected over the first one as their relative
        // shared measurement (2/4) is the same but its p-value is higher
        ASSERT_EQ(get_pattern(res_trk_cands, 0),
                  std::vector<std::size_t>({3, 5, 6, 13}));
    }
}

TEST(AmbiguitySolverTests, GreedyResolverTest3) {

    edm::measurement_collection::host measurements{host_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        host_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.2f, {5, 1, 11, 3});
    fill_pattern(trk_cands, 0.5f, {6, 2});
    fill_pattern(trk_cands, 0.4f, {3, 21, 12, 6, 19, 14});
    fill_pattern(trk_cands, 0.1f, {13, 16, 2, 7, 11});
    fill_pattern(trk_cands, 0.3f, {1, 7, 8});
    fill_pattern(trk_cands, 0.6f, {1, 3, 11, 22});

    traccc::host::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::host::greedy_ambiguity_resolution_algorithm resolution_alg(
        resolution_config, host_mr);

    {
        auto res_trk_cands = resolution_alg(
            traccc::edm::track_container<traccc::default_algebra>::const_data(
                trk_cands));
        ASSERT_EQ(res_trk_cands.tracks.size(), 2u);

        ASSERT_EQ(get_pattern(res_trk_cands, 0),
                  std::vector<std::size_t>({3, 21, 12, 6, 19, 14}));
        ASSERT_EQ(get_pattern(res_trk_cands, 1),
                  std::vector<std::size_t>({13, 16, 2, 7, 11}));
    }

    // Legacy algorithm
    traccc::legacy::greedy_ambiguity_resolution_algorithm::config_t legacy_cfg;
    traccc::legacy::greedy_ambiguity_resolution_algorithm legacy_resolution_alg(
        legacy_cfg, host_mr);

    {
        auto res_trk_cands = legacy_resolution_alg(trk_cands);
        ASSERT_EQ(res_trk_cands.tracks.size(), 2u);

        ASSERT_EQ(get_pattern(res_trk_cands, 0),
                  std::vector<std::size_t>({3, 21, 12, 6, 19, 14}));
        ASSERT_EQ(get_pattern(res_trk_cands, 1),
                  std::vector<std::size_t>({13, 16, 2, 7, 11}));
    }
}

// Comparison to the legacy algorithm.
TEST(AmbiguitySolverTests, GreedyResolverTest4) {

    edm::measurement_collection::host measurements{host_mr};
    const measurement_id_type max_meas_id = 10000;
    fill_measurements(measurements, max_meas_id);

    edm::track_container<default_algebra>::host trk_cands{
        host_mr, vecmem::get_data(measurements)};

    std::mt19937 gen(42);

    for (std::size_t i = 0; i < 10000u; i++) {

        std::uniform_int_distribution<std::size_t> track_length_dist(1, 20);
        std::uniform_int_distribution<measurement_id_type> meas_id_dist(
            0, max_meas_id);
        std::uniform_real_distribution<traccc::scalar> pval_dist(0.0f, 1.0f);

        const std::size_t track_length = track_length_dist(gen);
        const traccc::scalar pval = pval_dist(gen);
        std::vector<measurement_id_type> pattern;
        while (pattern.size() < track_length) {

            const measurement_id_type meas_id = meas_id_dist(gen);
            if (std::find(pattern.begin(), pattern.end(), meas_id) ==
                pattern.end()) {
                pattern.push_back(meas_id);
            }
        }

        std::sort(pattern.begin(), pattern.end());
        auto last = std::unique(pattern.begin(), pattern.end());

        // There should not be duplicate
        ASSERT_EQ(last, pattern.end());
        pattern.erase(last, pattern.end());

        // Make sure that partern size is eqaul to the track length
        ASSERT_EQ(pattern.size(), track_length);

        // Fill the pattern
        fill_pattern(trk_cands, pval, pattern);
    }

    traccc::host::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::host::greedy_ambiguity_resolution_algorithm resolution_alg(
        resolution_config, host_mr);

    auto start_new = std::chrono::high_resolution_clock::now();

    auto res_trk_cands = resolution_alg(
        traccc::edm::track_container<traccc::default_algebra>::const_data(
            trk_cands));

    auto end_new = std::chrono::high_resolution_clock::now();
    auto duration_new = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_new - start_new);
    std::cout << " Time for the new method " << duration_new.count() << " ms"
              << std::endl;

    // Legacy algorithm
    traccc::legacy::greedy_ambiguity_resolution_algorithm::config_t legacy_cfg;
    traccc::legacy::greedy_ambiguity_resolution_algorithm legacy_resolution_alg(
        legacy_cfg, host_mr);

    auto start_legacy = std::chrono::high_resolution_clock::now();

    auto legacy_res_trk_cands = legacy_resolution_alg(trk_cands);

    auto end_legacy = std::chrono::high_resolution_clock::now();
    auto duration_legacy =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_legacy -
                                                              start_legacy);
    std::cout << " Time for the legacy method " << duration_legacy.count()
              << " ms" << std::endl;

    std::size_t n_res_tracks = res_trk_cands.tracks.size();
    ASSERT_EQ(n_res_tracks, legacy_res_trk_cands.tracks.size());
    for (std::size_t i = 0; i < n_res_tracks; i++) {
        ASSERT_EQ(res_trk_cands.tracks.at(i),
                  legacy_res_trk_cands.tracks.at(i));
    }
}

TEST(AmbiguitySolverTests, GreedyResolverTest5) {

    edm::measurement_collection::host measurements{host_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        host_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.2f, {1, 2, 1, 1});
    fill_pattern(trk_cands, 0.5f, {3, 2, 1});
    fill_pattern(trk_cands, 0.4f, {2, 4, 5, 7, 2});
    fill_pattern(trk_cands, 0.1f, {6, 6, 6, 6});

    traccc::host::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::host::greedy_ambiguity_resolution_algorithm resolution_alg(
        resolution_config, host_mr);

    {
        auto res_trk_cands = resolution_alg(
            traccc::edm::track_container<traccc::default_algebra>::const_data(
                trk_cands));
        ASSERT_EQ(res_trk_cands.tracks.size(), 2u);

        ASSERT_EQ(get_pattern(res_trk_cands, 0),
                  std::vector<std::size_t>({3, 2, 1}));
        ASSERT_EQ(get_pattern(res_trk_cands, 1),
                  std::vector<std::size_t>({6, 6, 6, 6}));
    }
}

TEST(AmbiguitySolverTests, GreedyResolverTest6) {

    edm::measurement_collection::host measurements{host_mr};
    fill_measurements(measurements, 100);

    edm::track_container<default_algebra>::host trk_cands{
        host_mr, vecmem::get_data(measurements)};
    fill_pattern(trk_cands, 0.2f, {7, 3, 5, 7, 7, 7, 2});
    fill_pattern(trk_cands, 0.5f, {2});
    fill_pattern(trk_cands, 0.4f, {8, 9, 7, 2, 3, 4, 3, 7});
    fill_pattern(trk_cands, 0.1f, {8, 9, 0, 8, 1, 4, 6});
    fill_pattern(trk_cands, 0.9f, {10, 3, 2});

    traccc::host::greedy_ambiguity_resolution_algorithm::config_type
        resolution_config;
    traccc::host::greedy_ambiguity_resolution_algorithm resolution_alg(
        resolution_config, host_mr);

    {
        auto res_trk_cands = resolution_alg(
            traccc::edm::track_container<traccc::default_algebra>::const_data(
                trk_cands));
        ASSERT_EQ(res_trk_cands.tracks.size(), 2u);

        ASSERT_EQ(get_pattern(res_trk_cands, 0),
                  std::vector<std::size_t>({7, 3, 5, 7, 7, 7, 2}));
        ASSERT_EQ(get_pattern(res_trk_cands, 1),
                  std::vector<std::size_t>({8, 9, 0, 8, 1, 4, 6}));
    }
}
