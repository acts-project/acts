// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <Acts/Plugins/ExaTrkX/detail/ConnectedComponents.cuh>

#include <random>
#include <set>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

using namespace Acts::detail;

using BoostGraph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;

void checkLabeling(const std::vector<int> &src, const std::vector<int> &tgt) {
  std::size_t numNodes = std::max(*std::max_element(src.begin(), src.end()),
                                  *std::max_element(tgt.begin(), tgt.end())) +
                         1;

  int *cudaSrc, *cudaTgt;
  cudaMalloc(&cudaSrc, src.size() * sizeof(int));
  cudaMalloc(&cudaTgt, tgt.size() * sizeof(int));
  cudaMemcpy(cudaSrc, src.data(), src.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(cudaTgt, tgt.data(), src.size() * sizeof(int),
             cudaMemcpyHostToDevice);

  int *cudaLabels;
  cudaMalloc(&cudaLabels, numNodes * sizeof(int));
  int *cudaLabelsNext;
  cudaMalloc(&cudaLabelsNext, numNodes * sizeof(int));

  labelConnectedComponents<<<1, 1024>>>(src.size(), cudaSrc, cudaTgt, numNodes,
                                        cudaLabels, cudaLabelsNext);

  std::vector<int> labelsFromCuda(numNodes);
  cudaMemcpy(labelsFromCuda.data(), cudaLabels, numNodes * sizeof(int),
             cudaMemcpyDeviceToHost);

  BoostGraph G(numNodes);

  for (int i = 0; i < src.size(); ++i) {
    boost::add_edge(src[i], tgt[i], G);
  }

  std::vector<std::size_t> cpuLabels(numNodes);
  boost::connected_components(G, &cpuLabels[0]);

  // print
  std::cout << "cpu labels:     ";
  for (int i = 0; i < numNodes; ++i) {
    std::cout << cpuLabels[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "my CUDA labels: ";
  for (int i = 0; i < numNodes; ++i) {
    std::cout << labelsFromCuda[i] << " ";
  }
  std::cout << std::endl;

  // check systematically
  std::map<int, int> boostToCuda;
  for (int i = 0; i < numNodes; ++i) {
    if (boostToCuda.contains(cpuLabels[i])) {
      BOOST_CHECK_EQUAL(labelsFromCuda[i], boostToCuda.at(cpuLabels[i]));
    } else {
      auto [it, success] =
          boostToCuda.insert({cpuLabels[i], labelsFromCuda[i]});
      BOOST_CHECK(success);
    }
  }
}

using Vi = std::vector<int>;

BOOST_AUTO_TEST_CASE(simple_test_1) {
  Vi src{0, 1, 2, 3};
  Vi tgt{1, 2, 3, 4};
  checkLabeling(src, tgt);
}

BOOST_AUTO_TEST_CASE(simple_test_2) {
  Vi src{0, 1, 2, 4, 5, 6};
  Vi tgt{1, 2, 3, 5, 6, 7};
  checkLabeling(src, tgt);
}

BOOST_AUTO_TEST_CASE(simple_test_3) {
  Vi src{4, 3, 2, 1};
  Vi tgt{3, 2, 1, 0};
  checkLabeling(src, tgt);
}

BOOST_AUTO_TEST_CASE(test_label_mask) {
  Vi labels{0, 3, 5, 3, 0, 0};

  dim3 blockDim = 32;
  dim3 gridDim = (labels.size() + blockDim.x - 1) / blockDim.x;

  // Copy labels to device
  int *cudaLabels;
  cudaMalloc(&cudaLabels, labels.size() * sizeof(int));
  cudaMemcpy(cudaLabels, labels.data(), labels.size() * sizeof(int),
             cudaMemcpyHostToDevice);

  // Init label mask
  int *cudaLabelMask;
  cudaMalloc(&cudaLabelMask, labels.size() * sizeof(int));
  cudaMemset(cudaLabelMask, 0, labels.size() * sizeof(int));

  makeLabelMask<<<1, 256>>>(labels.size(), cudaLabels, cudaLabelMask);

  std::vector<int> labelMask(labels.size());
  cudaMemcpy(labelMask.data(), cudaLabelMask, 6 * sizeof(int),
             cudaMemcpyDeviceToHost);

  BOOST_CHECK((labelMask == Vi{1, 0, 0, 1, 0, 1}));

  // Do not test prefix sum here
  Vi prefixSum{0, 1, 1, 1, 2, 2};

  // copy prefix sum to device
  int *cudaPrefixSum;
  cudaMalloc(&cudaPrefixSum, prefixSum.size() * sizeof(int));
  cudaMemcpy(cudaPrefixSum, prefixSum.data(), prefixSum.size() * sizeof(int),
             cudaMemcpyHostToDevice);

  // Relabel
  mapEdgeLabels<<<1, 256>>>(labels.size(), cudaLabels, cudaPrefixSum);

  std::vector<int> labelsFromCuda(labels.size());
  cudaMemcpy(labelsFromCuda.data(), cudaLabels, labels.size() * sizeof(int),
             cudaMemcpyDeviceToHost);

  BOOST_CHECK((labelsFromCuda == Vi{0, 1, 2, 1, 0, 0}));
}

auto makeRandomGraph(std::size_t nodes, std::size_t edges) {
  std::default_random_engine rng(2345);
  std::uniform_int_distribution<> dist(0, nodes);
  std::set<std::pair<int, int>> set;
  Vi src(edges), tgt(edges);
  for (int n = 0; n < edges; ++n) {
    auto a = dist(rng);
    auto b = dist(rng);
    if (a == b) {
      continue;
    }
    auto s = std::min(a, b);
    auto t = std::max(a, b);
    auto [it, success] = set.insert({s, t});
    if (success) {
      src.at(n) = s;
      tgt.at(n) = t;
    }
  }

  return std::make_pair(src, tgt);
}

BOOST_AUTO_TEST_CASE(test_random_graph) {
  auto [src, tgt] = makeRandomGraph(5, 10);
  checkLabeling(src, tgt);
}

void testFullConnectedComponents(std::size_t nNodes, std::size_t nEdges) {
  auto [src, tgt] = makeRandomGraph(nNodes, nEdges);

  // the random graph can contain less edges than specified
  nEdges = src.size();

  // print src and tgt
  /*
    std::cout << "src: ";
    for (int i = 0; i < src.size(); ++i) {
      std::cout << src[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "tgt: ";
    for (int i = 0; i < tgt.size(); ++i) {
      std::cout << tgt[i] << " ";
    }
    std::cout << std::endl;
  */
  cudaStream_t stream;
  BOOST_REQUIRE_EQUAL(cudaStreamCreate(&stream), cudaSuccess);

  // copy src and tgt to device
  int *cudaSrc, *cudaTgt;
  BOOST_REQUIRE_EQUAL(
      cudaMallocAsync(&cudaSrc, src.size() * sizeof(int), stream), cudaSuccess);
  BOOST_REQUIRE_EQUAL(
      cudaMallocAsync(&cudaTgt, tgt.size() * sizeof(int), stream), cudaSuccess);
  BOOST_REQUIRE_EQUAL(
      cudaMemcpyAsync(cudaSrc, src.data(), src.size() * sizeof(int),
                      cudaMemcpyHostToDevice, stream),
      cudaSuccess);
  BOOST_REQUIRE_EQUAL(
      cudaMemcpyAsync(cudaTgt, tgt.data(), src.size() * sizeof(int),
                      cudaMemcpyHostToDevice, stream),
      cudaSuccess);

  // init label array
  int *cudaLabels;
  BOOST_REQUIRE_EQUAL(
      cudaMallocAsync(&cudaLabels, nNodes * sizeof(int), stream), cudaSuccess);

  // run connected components
  int cudaNumLabels = connectedComponentsCuda(src.size(), cudaSrc, cudaTgt,
                                              nNodes, cudaLabels, stream);
  BOOST_REQUIRE_EQUAL(cudaStreamSynchronize(stream), cudaSuccess);

  // print message from last cuda error code
  std::cout << "CUDA Error msg: " << cudaGetErrorString(cudaPeekAtLastError())
            << std::endl;
  BOOST_REQUIRE_EQUAL(cudaGetLastError(), cudaSuccess);

  // copy labels back
  std::vector<int> labelsFromCuda(nNodes);
  BOOST_REQUIRE_EQUAL(
      cudaMemcpyAsync(labelsFromCuda.data(), cudaLabels, nNodes * sizeof(int),
                      cudaMemcpyDeviceToHost, stream),
      cudaSuccess);

  BOOST_REQUIRE_EQUAL(cudaFreeAsync(cudaSrc, stream), cudaSuccess);
  BOOST_REQUIRE_EQUAL(cudaFreeAsync(cudaTgt, stream), cudaSuccess);
  BOOST_REQUIRE_EQUAL(cudaFreeAsync(cudaLabels, stream), cudaSuccess);

  // sync
  BOOST_REQUIRE_EQUAL(cudaStreamSynchronize(stream), cudaSuccess);
  BOOST_REQUIRE_EQUAL(cudaStreamDestroy(stream), cudaSuccess);

  // print labelsFromCuda
  /*
      std::cout << "CUDA labels: ";
      for (int i = 0; i < nNodes; ++i) {
        std::cout << labelsFromCuda[i] << " ";
      }
      std::cout << std::endl;
  */
  // run boost graph for comparison

  BoostGraph G(nNodes);

  for (int i = 0; i < src.size(); ++i) {
    boost::add_edge(src[i], tgt[i], G);
  }

  std::vector<std::size_t> cpuLabels(boost::num_vertices(G));
  int cpuNumLabels = boost::connected_components(G, &cpuLabels[0]);

  // check
  BOOST_CHECK_EQUAL(cudaNumLabels, cpuNumLabels);

  // check label vectors are the same
  for (int i = 0; i < nNodes; ++i) {
    BOOST_CHECK_EQUAL(labelsFromCuda[i], cpuLabels[i]);
  }
}

BOOST_AUTO_TEST_CASE(full_test_tiny_graph) {
  testFullConnectedComponents(5, 10);
}

BOOST_AUTO_TEST_CASE(full_test_small_graph) {
  testFullConnectedComponents(100, 500);
}

BOOST_AUTO_TEST_CASE(full_test_big_graph) {
  testFullConnectedComponents(100'000, 500'000);
}
