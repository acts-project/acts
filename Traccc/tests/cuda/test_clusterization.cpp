/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "tests/cca_test.hpp"
#include "traccc/clusterization/clustering_config.hpp"
#include "traccc/cuda/clusterization/clusterization_algorithm.hpp"
#include "traccc/definitions/common.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/performance/collection_comparator.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/stream_wrapper.hpp>

// GTest include(s).
#include <gtest/gtest.h>

using namespace traccc;

namespace {
void run_clustering_test(
    const traccc::edm::silicon_cell_collection::host& cells,
    const edm::measurement_collection::host& references,
    const std::vector<std::vector<unsigned int>>& reference_disjoint_set,
    auto& mng_mr, const auto& cfg) {
  traccc::memory_resource mr{mng_mr};

  // Cuda stream
  vecmem::cuda::stream_wrapper vecmem_stream;
  traccc::cuda::stream_wrapper stream{vecmem_stream.stream()};

  // Cuda copy objects
  vecmem::cuda::async_copy copy{stream.cudaStream()};

  // Create a dummy detector description.
  traccc::detector_design_description::host det_desc{mng_mr};
  traccc::detector_conditions_description::host det_cond{mng_mr};
  det_desc.resize(1u);
  det_cond.resize(1u);
  det_desc.bin_edges_x()[0] = {0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f};
  det_desc.bin_edges_y()[0] = {0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f};
  det_desc.dimensions()[0] = 2;
  det_cond.geometry_id()[0] = detray::geometry::identifier{0u};
  det_cond.measurement_translation()[0] = {0.f, 0.f};

  // Run Clusterization
  traccc::cuda::clusterization_algorithm ca_cuda(mr, copy, stream, cfg);

  auto measurements_buffer =
      ca_cuda(vecmem::get_data(cells), vecmem::get_data(det_desc),
              vecmem::get_data(det_cond));
  auto [measurements_buffer_wdjs, disjoint_set] = ca_cuda(
      vecmem::get_data(cells), vecmem::get_data(det_desc),
      vecmem::get_data(det_cond), device::clustering_keep_disjoint_set{});
  edm::measurement_collection::const_device measurements(measurements_buffer);
  edm::measurement_collection::const_device measurements_wdjs(
      measurements_buffer_wdjs);

  // Check the results
  ASSERT_EQ(copy.get_size(measurements_buffer), references.size());
  ASSERT_EQ(copy.get_size(measurements_buffer_wdjs), references.size());

  auto check_all_matched = [&references](const auto& meas) {
    for (unsigned int i = 0; i < meas.size(); ++i) {
      const auto test = meas.at(i);
      // 0.01 % uncertainty
      auto iso = traccc::details::is_same_object<
          edm::measurement_collection::const_device::object_type>(test,
                                                                  0.0001f);
      bool matched = false;

      for (std::size_t j = 0; j < references.size(); ++j) {
        const auto ref = references.at(j);
        if (iso(ref)) {
          matched = true;
          break;
        }
      }

      ASSERT_TRUE(matched);
    }
  };

  check_all_matched(measurements);
  check_all_matched(measurements_wdjs);

  edm::silicon_cluster_collection::const_view disjoint_set_view{disjoint_set};
  edm::silicon_cluster_collection::const_device disjoint_set_device{
      disjoint_set_view};

  ASSERT_EQ(reference_disjoint_set.size(), disjoint_set_device.size());

  auto disjoint_set_matched = [&disjoint_set_device](const auto& ref) {
    std::vector<unsigned int> reference_set = ref;
    std::sort(reference_set.begin(), reference_set.end());

    for (unsigned int i = 0; i < disjoint_set_device.size(); ++i) {
      std::vector<unsigned int> found_set;
      for (unsigned int j = 0;
           j < disjoint_set_device.cell_indices().at(i).size(); ++j) {
        found_set.push_back(disjoint_set_device.cell_indices().at(i).at(j));
      }
      std::sort(found_set.begin(), found_set.end());

      if (reference_set.size() != found_set.size()) {
        continue;
      }

      bool equal = true;
      for (unsigned int j = 0; j < reference_set.size(); ++j) {
        if (reference_set.at(j) != found_set.at(j)) {
          equal = false;
          break;
        }
      }

      if (equal) {
        return true;
      }
    }

    return false;
  };

  for (const auto& i : reference_disjoint_set) {
    ASSERT_TRUE(disjoint_set_matched(i));
  }
}
}  // namespace

TEST(CUDAClustering, SingleModule) {
  // Memory resource used by the EDM.
  vecmem::cuda::managed_memory_resource mng_mr;

  // Create cell collection
  traccc::edm::silicon_cell_collection::host cells{mng_mr};
  cells.reserve(8u);
  cells.push_back({1u, 2u, 1.f, 0.f, 0u});
  cells.push_back({2u, 2u, 1.f, 0.f, 0u});
  cells.push_back({3u, 2u, 1.f, 0.f, 0u});
  cells.push_back({6u, 4u, 1.f, 0.f, 0u});
  cells.push_back({5u, 5u, 1.f, 0.f, 0u});
  cells.push_back({6u, 5u, 1.f, 0.f, 0u});
  cells.push_back({7u, 5u, 1.f, 0.f, 0u});
  cells.push_back({6u, 6u, 1.f, 0.f, 0u});

  edm::measurement_collection::host references{mng_mr};
  references.push_back({{2.5f, 2.5f},
                        {0.75f, 0.0833333f},
                        2u,
                        0.f,
                        0.f,
                        0u,
                        detray::geometry::identifier{0u},
                        {1u, 1u},
                        0u});
  references.push_back({{6.5f, 5.5f},
                        {0.483333f, 0.483333f},
                        2u,
                        0.f,
                        0.f,
                        0u,
                        detray::geometry::identifier{0u},
                        {1u, 1u},
                        1u});

  std::vector<std::vector<unsigned int>> reference_disjoint_set{
      {0, 1, 2}, {3, 4, 5, 6, 7}};

  auto cfg = default_ccl_test_config();

  run_clustering_test(cells, references, reference_disjoint_set, mng_mr, cfg);
}

TEST(CUDAClustering, SingleModuleUnsorted) {
  // Memory resource used by the EDM.
  vecmem::cuda::managed_memory_resource mng_mr;

  // Create cell collection
  traccc::edm::silicon_cell_collection::host cells{mng_mr};
  cells.reserve(8u);
  cells.push_back({7u, 5u, 1.f, 0.f, 0u});
  cells.push_back({3u, 2u, 1.f, 0.f, 0u});
  cells.push_back({5u, 5u, 1.f, 0.f, 0u});
  cells.push_back({6u, 5u, 1.f, 0.f, 0u});
  cells.push_back({2u, 2u, 1.f, 0.f, 0u});
  cells.push_back({1u, 2u, 1.f, 0.f, 0u});
  cells.push_back({6u, 6u, 1.f, 0.f, 0u});
  cells.push_back({6u, 4u, 1.f, 0.f, 0u});

  edm::measurement_collection::host references{mng_mr};
  references.push_back({{2.5f, 2.5f},
                        {0.75f, 0.0833333f},
                        2u,
                        0.f,
                        0.f,
                        0u,
                        detray::geometry::identifier{0u},
                        {1u, 1u},
                        0u});
  references.push_back({{6.5f, 5.5f},
                        {0.483333f, 0.483333f},
                        2u,
                        0.f,
                        0.f,
                        0u,
                        detray::geometry::identifier{0u},
                        {1u, 1u},
                        1u});

  std::vector<std::vector<unsigned int>> reference_disjoint_set{
      {5, 4, 1}, {7, 2, 3, 0, 6}};

  auto cfg = default_ccl_test_config();
  cfg.sort_cells = true;

  run_clustering_test(cells, references, reference_disjoint_set, mng_mr, cfg);
}
