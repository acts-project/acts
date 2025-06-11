// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <Acts/Plugins/ExaTrkX/detail/CudaUtils.hpp>
#include <Acts/Plugins/ExaTrkX/detail/JunctionRemoval.hpp>

#include <algorithm>
#include <numeric>

using Vi = std::vector<std::int64_t>;
using Vf = std::vector<float>;

void testJunctionRemoval(const Vi &srcNodes, const Vi &dstNodes,
                         const Vf &scores, const Vi &expectedSrcNodes,
                         const Vi &expectedDstNodes) {
  std::size_t nEdges = srcNodes.size();
  std::size_t nNodes =
      std::max(*std::max_element(srcNodes.begin(), srcNodes.end()),
               *std::max_element(dstNodes.begin(), dstNodes.end())) +
      1;
  std::size_t nExpectedEdges = expectedSrcNodes.size();

  cudaStream_t stream{};
  ACTS_CUDA_CHECK(cudaStreamCreate(&stream));

  std::int64_t *cudaSrcNodes{}, *cudaDstNodes{};
  float *cudaScores{};
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&cudaSrcNodes, nEdges * sizeof(std::int64_t), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&cudaDstNodes, nEdges * sizeof(std::int64_t), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaScores, nEdges * sizeof(float), stream));

  ACTS_CUDA_CHECK(cudaMemcpyAsync(cudaSrcNodes, srcNodes.data(),
                                  nEdges * sizeof(std::int64_t),
                                  cudaMemcpyHostToDevice, stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(cudaDstNodes, dstNodes.data(),
                                  nEdges * sizeof(std::int64_t),
                                  cudaMemcpyHostToDevice, stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(cudaScores, scores.data(),
                                  nEdges * sizeof(float),
                                  cudaMemcpyHostToDevice, stream));

  auto [cudaSrcNodesOut, nEdgesOut] = Acts::detail::junctionRemovalCuda(
      nEdges, nNodes, cudaScores, cudaSrcNodes, cudaDstNodes, stream);
  auto cudaDstNodesOut = cudaSrcNodesOut + nEdgesOut;

  Vi srcNodesOut(nEdgesOut);
  Vi dstNodesOut(nEdgesOut);
  ACTS_CUDA_CHECK(cudaMemcpyAsync(srcNodesOut.data(), cudaSrcNodesOut,
                                  nEdgesOut * sizeof(std::int64_t),
                                  cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(dstNodesOut.data(), cudaDstNodesOut,
                                  nEdgesOut * sizeof(std::int64_t),
                                  cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  ACTS_CUDA_CHECK(cudaFreeAsync(cudaSrcNodes, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaDstNodes, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaScores, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaSrcNodesOut, stream));

  ACTS_CUDA_CHECK(cudaStreamDestroy(stream));

  // Sort edges by src and dst nodes before comparison
  Vi idxs(nEdgesOut);
  std::iota(idxs.begin(), idxs.end(), 0);
  std::sort(idxs.begin(), idxs.end(), [&](int a, int b) {
    if (srcNodesOut.at(a) != srcNodesOut.at(b)) {
      return srcNodesOut.at(a) < srcNodesOut.at(b);
    }
    return dstNodesOut.at(a) < dstNodesOut.at(b);
  });

  Vi srcNodesOutSorted(nEdgesOut);
  Vi dstNodesOutSorted(nEdgesOut);
  for (std::size_t i = 0; i < nEdgesOut; ++i) {
    srcNodesOutSorted.at(i) = srcNodesOut.at(idxs.at(i));
    dstNodesOutSorted.at(i) = dstNodesOut.at(idxs.at(i));
  }

  BOOST_REQUIRE_EQUAL(nEdgesOut, nExpectedEdges);
  BOOST_CHECK_EQUAL_COLLECTIONS(
      srcNodesOutSorted.begin(), srcNodesOutSorted.end(),
      expectedSrcNodes.begin(), expectedSrcNodes.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(
      dstNodesOutSorted.begin(), dstNodesOutSorted.end(),
      expectedDstNodes.begin(), expectedDstNodes.end());
}

BOOST_AUTO_TEST_CASE(test_no_junction) {
  Vi srcNodes = {0, 1};
  Vi dstNodes = {1, 2};
  Vf scores = {0.3f, 0.9f};
  Vi expectedSrcNodes = {0, 1};
  Vi expectedDstNodes = {1, 2};
  testJunctionRemoval(srcNodes, dstNodes, scores, expectedSrcNodes,
                      expectedDstNodes);
}

BOOST_AUTO_TEST_CASE(test_junction_2_in) {
  Vi srcNodes = {0, 1, 2};
  Vi dstNodes = {2, 2, 3};
  Vf scores = {0.9f, 0.3f, 0.9f};
  Vi expectedSrcNodes = {0, 2};
  Vi expectedDstNodes = {2, 3};
  testJunctionRemoval(srcNodes, dstNodes, scores, expectedSrcNodes,
                      expectedDstNodes);
}

BOOST_AUTO_TEST_CASE(test_junction_2_out) {
  Vi srcNodes = {0, 1, 1};
  Vi dstNodes = {1, 2, 3};
  Vf scores = {0.3f, 0.3f, 0.9f};
  Vi expectedSrcNodes = {0, 1};
  Vi expectedDstNodes = {1, 3};
  testJunctionRemoval(srcNodes, dstNodes, scores, expectedSrcNodes,
                      expectedDstNodes);
}

BOOST_AUTO_TEST_CASE(test_junction_2_in_2_out) {
  Vi srcNodes = {0, 1, 2, 2};
  Vi dstNodes = {2, 2, 3, 4};
  Vf scores = {0.9f, 0.3f, 0.9f, 0.5f};
  Vi expectedSrcNodes = {0, 2};
  Vi expectedDstNodes = {2, 3};
  testJunctionRemoval(srcNodes, dstNodes, scores, expectedSrcNodes,
                      expectedDstNodes);
}

BOOST_AUTO_TEST_CASE(test_junction_3_in_3_out) {
  Vi srcNodes = {0, 1, 2, 3, 3, 3};
  Vi dstNodes = {3, 3, 3, 4, 5, 6};
  Vf scores = {0.2f, 0.3f, 0.9f, 0.5f, 0.9f, 0.1f};
  Vi expectedSrcNodes = {2, 3};
  Vi expectedDstNodes = {3, 5};
  testJunctionRemoval(srcNodes, dstNodes, scores, expectedSrcNodes,
                      expectedDstNodes);
}

BOOST_AUTO_TEST_CASE(test_junction_leftover_edge) {
  int j = 2;
  Vi srcNodes = {0, 1, 3, j};
  Vi dstNodes = {1, j, j, 4};
  Vf scores = {0.9f, 0.1f, 0.9f, 0.9f};
  Vi expectedSrcNodes = {0, j, 3};
  Vi expectedDstNodes = {1, 4, j};
  testJunctionRemoval(srcNodes, dstNodes, scores, expectedSrcNodes,
                      expectedDstNodes);
}

BOOST_AUTO_TEST_CASE(test_no_in_edges) {
  int j = 0;
  Vi srcNodes = {j, j};
  Vi dstNodes = {1, 2};
  Vf scores = {0.9f, 0.1f};
  Vi expectedSrcNodes = {j};
  Vi expectedDstNodes = {1};
  testJunctionRemoval(srcNodes, dstNodes, scores, expectedSrcNodes,
                      expectedDstNodes);
}

BOOST_AUTO_TEST_CASE(test_no_out_edges) {
  int j = 0;
  Vi srcNodes = {1, 2};
  Vi dstNodes = {j, j};
  Vf scores = {0.9f, 0.1f};
  Vi expectedSrcNodes = {1};
  Vi expectedDstNodes = {j};
  testJunctionRemoval(srcNodes, dstNodes, scores, expectedSrcNodes,
                      expectedDstNodes);
}

BOOST_AUTO_TEST_CASE(test_two_junctions) {
  int j1 = 2;
  int j2 = 3;
  Vi srcNodes = {0, 1, j1, j2, j2};
  Vi dstNodes = {j1, j1, j2, 4, 5};
  Vf scores = {0.9f, 0.1f, 0.9f, 0.1f, 0.9f};
  Vi expectedSrcNodes = {0, j1, j2};
  Vi expectedDstNodes = {j1, j2, 5};

  testJunctionRemoval(srcNodes, dstNodes, scores, expectedSrcNodes,
                      expectedDstNodes);
}
