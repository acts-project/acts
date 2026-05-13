// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsPlugins/Gnn/Stages.hpp"

using namespace ActsPlugins;

namespace ActsTests {

namespace {

PipelineTensors makeTensors(std::size_t nNodes, std::size_t nFeatures,
                            std::vector<std::int64_t> edges) {
  ExecutionContext ctx{Device::Cpu(), {}};
  const auto nEdges = edges.size() / 2;

  auto nodeFeatures = Tensor<float>::Create({nNodes, nFeatures}, ctx);
  for (std::size_t i = 0; i < nNodes * nFeatures; ++i) {
    nodeFeatures.data()[i] = static_cast<float>(i);
  }

  auto edgeIndex = Tensor<std::int64_t>::Create({2, nEdges}, ctx);
  std::copy(edges.begin(), edges.end(), edgeIndex.data());

  return {std::move(nodeFeatures), std::move(edgeIndex), std::nullopt,
          std::nullopt};
}

}  // namespace

BOOST_AUTO_TEST_SUITE(GnnRemoveUnusedNodesSuite)

BOOST_AUTO_TEST_CASE(test_all_nodes_used) {
  // 3 nodes, all referenced: edges 0→1, 1→2
  auto tensors = makeTensors(3, 4, {0, 1, 1, 2});
  const auto originalNodeData = std::vector<float>(
      tensors.nodeFeatures.data(),
      tensors.nodeFeatures.data() + tensors.nodeFeatures.size());

  std::vector<int> spIds = {10, 11, 12};
  ExecutionContext ctx{Device::Cpu(), {}};
  auto result = removeUnusedNodes(std::move(tensors), spIds, ctx);

  // All 3 nodes survive
  BOOST_CHECK_EQUAL(result.nodeFeatures.shape()[0], 3u);
  BOOST_CHECK_EQUAL(result.nodeFeatures.shape()[1], 4u);
  BOOST_CHECK_EQUAL(result.edgeIndex.shape()[1], 2u);

  // spacePointIds unchanged
  BOOST_REQUIRE_EQUAL(spIds.size(), 3u);
  BOOST_CHECK_EQUAL(spIds[0], 10);
  BOOST_CHECK_EQUAL(spIds[1], 11);
  BOOST_CHECK_EQUAL(spIds[2], 12);

  // Node features unchanged
  BOOST_CHECK_EQUAL_COLLECTIONS(
      result.nodeFeatures.data(),
      result.nodeFeatures.data() + result.nodeFeatures.size(),
      originalNodeData.begin(), originalNodeData.end());

  // Edges remain 0→1, 1→2
  BOOST_CHECK_EQUAL(result.edgeIndex.data()[0], 0);
  BOOST_CHECK_EQUAL(result.edgeIndex.data()[1], 1);
  BOOST_CHECK_EQUAL(result.edgeIndex.data()[2], 1);
  BOOST_CHECK_EQUAL(result.edgeIndex.data()[3], 2);
}

BOOST_AUTO_TEST_CASE(test_unused_nodes_removed) {
  // 5 nodes, but nodes 1 and 3 are never referenced
  // edges: 0→2, 2→4  → used nodes: {0, 2, 4} → new indices 0, 1, 2
  auto tensors = makeTensors(5, 2, {0, 2, 2, 4});
  std::vector<int> spIds = {100, 101, 102, 103, 104};

  ExecutionContext ctx{Device::Cpu(), {}};
  auto result = removeUnusedNodes(std::move(tensors), spIds, ctx);

  // Only 3 nodes should remain
  BOOST_CHECK_EQUAL(result.nodeFeatures.shape()[0], 3u);
  BOOST_CHECK_EQUAL(result.nodeFeatures.shape()[1], 2u);
  BOOST_CHECK_EQUAL(result.edgeIndex.shape()[1], 2u);

  // spacePointIds remapped: original 0→100, 2→102, 4→104
  BOOST_REQUIRE_EQUAL(spIds.size(), 3u);
  BOOST_CHECK_EQUAL(spIds[0], 100);
  BOOST_CHECK_EQUAL(spIds[1], 102);
  BOOST_CHECK_EQUAL(spIds[2], 104);

  // Edges renumbered: 0→1, 1→2 (was 0→2, 2→4)
  BOOST_CHECK_EQUAL(result.edgeIndex.data()[0], 0);  // source of edge 0
  BOOST_CHECK_EQUAL(result.edgeIndex.data()[1], 1);  // source of edge 1
  BOOST_CHECK_EQUAL(result.edgeIndex.data()[2], 1);  // dest of edge 0
  BOOST_CHECK_EQUAL(result.edgeIndex.data()[3], 2);  // dest of edge 1

  // Node features of kept nodes (original rows 0, 2, 4)
  // With nFeatures=2: row i has values {2*i, 2*i+1}
  BOOST_CHECK_EQUAL(result.nodeFeatures.data()[0], 0.f);  // old row 0, col 0
  BOOST_CHECK_EQUAL(result.nodeFeatures.data()[1], 1.f);  // old row 0, col 1
  BOOST_CHECK_EQUAL(result.nodeFeatures.data()[2], 4.f);  // old row 2, col 0
  BOOST_CHECK_EQUAL(result.nodeFeatures.data()[3], 5.f);  // old row 2, col 1
  BOOST_CHECK_EQUAL(result.nodeFeatures.data()[4], 8.f);  // old row 4, col 0
  BOOST_CHECK_EQUAL(result.nodeFeatures.data()[5], 9.f);  // old row 4, col 1
}

BOOST_AUTO_TEST_CASE(test_zero_edges_throws) {
  // edgeIndex with 0 columns — must throw NoEdgesError
  ExecutionContext ctx{Device::Cpu(), {}};
  auto nodeFeatures = Tensor<float>::Create({4, 2}, ctx);
  auto edgeIndex = Tensor<std::int64_t>::Create({2, 0}, ctx);
  PipelineTensors tensors{std::move(nodeFeatures), std::move(edgeIndex),
                          std::nullopt, std::nullopt};
  std::vector<int> spIds = {0, 1, 2, 3};

  BOOST_CHECK_THROW(removeUnusedNodes(std::move(tensors), spIds, ctx),
                    NoEdgesError);
}

BOOST_AUTO_TEST_CASE(test_edge_and_score_tensors_preserved) {
  // Ensure edgeFeatures and edgeScores pass through unchanged
  // 4 nodes, edge 1→2 only (nodes 0 and 3 unused)
  ExecutionContext ctx{Device::Cpu(), {}};

  auto nodeFeatures = Tensor<float>::Create({4, 1}, ctx);
  auto edgeIndex = Tensor<std::int64_t>::Create({2, 1}, ctx);
  edgeIndex.data()[0] = 1;  // source
  edgeIndex.data()[1] = 2;  // dest

  auto edgeFeatures = Tensor<float>::Create({1, 3}, ctx);
  edgeFeatures.data()[0] = 7.f;
  edgeFeatures.data()[1] = 8.f;
  edgeFeatures.data()[2] = 9.f;

  auto edgeScores = Tensor<float>::Create({1, 1}, ctx);
  edgeScores.data()[0] = 0.99f;

  PipelineTensors tensors{std::move(nodeFeatures), std::move(edgeIndex),
                          std::move(edgeFeatures), std::move(edgeScores)};

  std::vector<int> spIds = {10, 11, 12, 13};
  auto result = removeUnusedNodes(std::move(tensors), spIds, ctx);

  // Nodes 1 and 2 survive → new indices 0 and 1
  BOOST_CHECK_EQUAL(result.nodeFeatures.shape()[0], 2u);
  BOOST_REQUIRE_EQUAL(spIds.size(), 2u);
  BOOST_CHECK_EQUAL(spIds[0], 11);
  BOOST_CHECK_EQUAL(spIds[1], 12);

  // Edge renumbered: 0→1 (was 1→2)
  BOOST_CHECK_EQUAL(result.edgeIndex.data()[0], 0);
  BOOST_CHECK_EQUAL(result.edgeIndex.data()[1], 1);

  // edgeFeatures unchanged
  BOOST_REQUIRE(result.edgeFeatures.has_value());
  BOOST_CHECK_EQUAL(result.edgeFeatures->data()[0], 7.f);
  BOOST_CHECK_EQUAL(result.edgeFeatures->data()[1], 8.f);
  BOOST_CHECK_EQUAL(result.edgeFeatures->data()[2], 9.f);

  // edgeScores unchanged
  BOOST_REQUIRE(result.edgeScores.has_value());
  BOOST_CHECK_EQUAL(result.edgeScores->data()[0], 0.99f);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
