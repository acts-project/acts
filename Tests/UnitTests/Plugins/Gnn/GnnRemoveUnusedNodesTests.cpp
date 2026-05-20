// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsPlugins/Gnn/Stages.hpp"

#ifdef ACTS_GNN_WITH_CUDA
#include <cuda_runtime_api.h>
#endif

using namespace ActsPlugins;

namespace ActsTests {

namespace {

const ExecutionContext execContextCpu{Device::Cpu(), {}};
#ifdef ACTS_GNN_WITH_CUDA
const ExecutionContext execContextCuda{Device::Cuda(0), cudaStreamLegacy};
#endif

// Build a PipelineTensors on the given device from flat CPU data.
PipelineTensors makeTensors(std::size_t nNodes, std::size_t nFeatures,
                            std::vector<std::int64_t> edges,
                            const ExecutionContext &ctx) {
  const auto nEdges = edges.size() / 2;

  auto nodeFeaturesCpu =
      Tensor<float>::Create({nNodes, nFeatures}, execContextCpu);
  for (std::size_t i = 0; i < nNodes * nFeatures; ++i) {
    nodeFeaturesCpu.data()[i] = static_cast<float>(i);
  }

  auto edgeIndexCpu = Tensor<std::int64_t>::Create({2, nEdges}, execContextCpu);
  std::copy(edges.begin(), edges.end(), edgeIndexCpu.data());

  return {std::move(nodeFeaturesCpu).clone(ctx),
          std::move(edgeIndexCpu).clone(ctx), std::nullopt, std::nullopt};
}

// Read tensor data back to a host vector regardless of device.
template <typename T>
std::vector<T> toHost(const Tensor<T> &tensor, const ExecutionContext &srcCtx) {
  auto host = tensor.clone({Device::Cpu(), srcCtx.stream});
  return {host.data(), host.data() + host.size()};
}

void testAllNodesUsed(const ExecutionContext &ctx) {
  // 3 nodes, all referenced: edges 0→1, 1→2
  auto tensors = makeTensors(3, 4, {0, 1, 1, 2}, ctx);
  const auto originalNodeData = toHost(tensors.nodeFeatures, ctx);

  std::vector<int> spIds = {10, 11, 12};
  auto result = removeUnusedNodes(std::move(tensors), spIds, ctx);

  BOOST_CHECK_EQUAL(result.nodeFeatures.shape()[0], 3u);
  BOOST_CHECK_EQUAL(result.nodeFeatures.shape()[1], 4u);
  BOOST_CHECK_EQUAL(result.edgeIndex.shape()[1], 2u);

  BOOST_REQUIRE_EQUAL(spIds.size(), 3u);
  BOOST_CHECK_EQUAL(spIds[0], 10);
  BOOST_CHECK_EQUAL(spIds[1], 11);
  BOOST_CHECK_EQUAL(spIds[2], 12);

  const auto nodeData = toHost(result.nodeFeatures, ctx);
  BOOST_CHECK_EQUAL_COLLECTIONS(nodeData.begin(), nodeData.end(),
                                originalNodeData.begin(),
                                originalNodeData.end());

  const auto edgeData = toHost(result.edgeIndex, ctx);
  BOOST_CHECK_EQUAL(edgeData[0], 0);
  BOOST_CHECK_EQUAL(edgeData[1], 1);
  BOOST_CHECK_EQUAL(edgeData[2], 1);
  BOOST_CHECK_EQUAL(edgeData[3], 2);
}

void testUnusedNodesRemoved(const ExecutionContext &ctx) {
  // 5 nodes, nodes 1 and 3 unreferenced: edges 0→2, 2→4
  // Used: {0,2,4} → compacted to {0,1,2}
  auto tensors = makeTensors(5, 2, {0, 2, 2, 4}, ctx);
  std::vector<int> spIds = {100, 101, 102, 103, 104};

  auto result = removeUnusedNodes(std::move(tensors), spIds, ctx);

  BOOST_CHECK_EQUAL(result.nodeFeatures.shape()[0], 3u);
  BOOST_CHECK_EQUAL(result.nodeFeatures.shape()[1], 2u);
  BOOST_CHECK_EQUAL(result.edgeIndex.shape()[1], 2u);

  BOOST_REQUIRE_EQUAL(spIds.size(), 3u);
  BOOST_CHECK_EQUAL(spIds[0], 100);
  BOOST_CHECK_EQUAL(spIds[1], 102);
  BOOST_CHECK_EQUAL(spIds[2], 104);

  const auto edgeData = toHost(result.edgeIndex, ctx);
  BOOST_CHECK_EQUAL(edgeData[0], 0);  // source of edge 0
  BOOST_CHECK_EQUAL(edgeData[1], 1);  // source of edge 1
  BOOST_CHECK_EQUAL(edgeData[2], 1);  // dest of edge 0
  BOOST_CHECK_EQUAL(edgeData[3], 2);  // dest of edge 1

  // With nFeatures=2, row i has values {2*i, 2*i+1}; kept rows: 0, 2, 4
  const auto nodeData = toHost(result.nodeFeatures, ctx);
  BOOST_CHECK_EQUAL(nodeData[0], 0.f);  // old row 0, col 0
  BOOST_CHECK_EQUAL(nodeData[1], 1.f);  // old row 0, col 1
  BOOST_CHECK_EQUAL(nodeData[2], 4.f);  // old row 2, col 0
  BOOST_CHECK_EQUAL(nodeData[3], 5.f);  // old row 2, col 1
  BOOST_CHECK_EQUAL(nodeData[4], 8.f);  // old row 4, col 0
  BOOST_CHECK_EQUAL(nodeData[5], 9.f);  // old row 4, col 1
}

void testEdgeAndScoreTensorsPreserved(const ExecutionContext &ctx) {
  // 4 nodes, single edge 1→2; nodes 0 and 3 unused
  auto nodeFeaturesCpu = Tensor<float>::Create({4, 1}, execContextCpu);
  auto edgeIndexCpu = Tensor<std::int64_t>::Create({2, 1}, execContextCpu);
  edgeIndexCpu.data()[0] = 1;
  edgeIndexCpu.data()[1] = 2;
  auto edgeFeaturesCpu = Tensor<float>::Create({1, 3}, execContextCpu);
  edgeFeaturesCpu.data()[0] = 7.f;
  edgeFeaturesCpu.data()[1] = 8.f;
  edgeFeaturesCpu.data()[2] = 9.f;
  auto edgeScoresCpu = Tensor<float>::Create({1, 1}, execContextCpu);
  edgeScoresCpu.data()[0] = 0.99f;

  PipelineTensors tensors{nodeFeaturesCpu.clone(ctx), edgeIndexCpu.clone(ctx),
                          edgeFeaturesCpu.clone(ctx), edgeScoresCpu.clone(ctx)};

  std::vector<int> spIds = {10, 11, 12, 13};
  auto result = removeUnusedNodes(std::move(tensors), spIds, ctx);

  BOOST_CHECK_EQUAL(result.nodeFeatures.shape()[0], 2u);
  BOOST_REQUIRE_EQUAL(spIds.size(), 2u);
  BOOST_CHECK_EQUAL(spIds[0], 11);
  BOOST_CHECK_EQUAL(spIds[1], 12);

  const auto edgeData = toHost(result.edgeIndex, ctx);
  BOOST_CHECK_EQUAL(edgeData[0], 0);
  BOOST_CHECK_EQUAL(edgeData[1], 1);

  BOOST_REQUIRE(result.edgeFeatures.has_value());
  const auto featData = toHost(*result.edgeFeatures, ctx);
  BOOST_CHECK_EQUAL(featData[0], 7.f);
  BOOST_CHECK_EQUAL(featData[1], 8.f);
  BOOST_CHECK_EQUAL(featData[2], 9.f);

  BOOST_REQUIRE(result.edgeScores.has_value());
  const auto scoreData = toHost(*result.edgeScores, ctx);
  BOOST_CHECK_EQUAL(scoreData[0], 0.99f);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(GnnRemoveUnusedNodesSuite)

BOOST_AUTO_TEST_CASE(test_all_nodes_used_cpu) {
  testAllNodesUsed(execContextCpu);
}

BOOST_AUTO_TEST_CASE(test_unused_nodes_removed_cpu) {
  testUnusedNodesRemoved(execContextCpu);
}

BOOST_AUTO_TEST_CASE(test_zero_edges_throws) {
  auto nodeFeatures = Tensor<float>::Create({4, 2}, execContextCpu);
  auto edgeIndex = Tensor<std::int64_t>::Create({2, 0}, execContextCpu);
  PipelineTensors tensors{std::move(nodeFeatures), std::move(edgeIndex),
                          std::nullopt, std::nullopt};
  std::vector<int> spIds = {0, 1, 2, 3};
  BOOST_CHECK_THROW(
      removeUnusedNodes(std::move(tensors), spIds, execContextCpu),
      NoEdgesError);
}

BOOST_AUTO_TEST_CASE(test_edge_and_score_tensors_preserved_cpu) {
  testEdgeAndScoreTensorsPreserved(execContextCpu);
}

#ifdef ACTS_GNN_WITH_CUDA

BOOST_AUTO_TEST_CASE(test_all_nodes_used_cuda) {
  testAllNodesUsed(execContextCuda);
}

BOOST_AUTO_TEST_CASE(test_unused_nodes_removed_cuda) {
  testUnusedNodesRemoved(execContextCuda);
}

BOOST_AUTO_TEST_CASE(test_edge_and_score_tensors_preserved_cuda) {
  testEdgeAndScoreTensorsPreserved(execContextCuda);
}

#endif  // ACTS_GNN_WITH_CUDA

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
