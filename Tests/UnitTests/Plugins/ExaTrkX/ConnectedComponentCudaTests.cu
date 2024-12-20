// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <Acts/Plugins/ExaTrkX/detail/connectedComponents.cuh>

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

auto makeRandomGraph(std::size_t N) {
  std::default_random_engine rng(2345);
  std::uniform_int_distribution<> dist(0, 20);
  Vi src, tgt;
  for (int n = 0; n < N; ++n) {
    src.push_back(dist(rng));
    tgt.push_back(dist(rng));
  }

  return std::make_pair(src, tgt);
}

BOOST_AUTO_TEST_CASE(test_random_graph) {
  auto [src, tgt] = makeRandomGraph(10);
  checkLabeling(src, tgt);
}

BOOST_AUTO_TEST_CASE(test_random_graph_with_relabeling) {
  std::size_t nNodes = 100;
  auto [src, tgt] = makeRandomGraph(nNodes);

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

  // run boost graph for comparison
  BoostGraph G(nNodes);

  for (int i = 0; i < src.size(); ++i) {
    boost::add_edge(src[i], tgt[i], G);
  }

  std::vector<std::size_t> cpuLabels(nNodes);
  int cpuNumLabels = boost::connected_components(G, &cpuLabels[0]);

  // check
  BOOST_CHECK_EQUAL(cudaNumLabels, cpuNumLabels);

  // check label vectors are the same
  for (int i = 0; i < nNodes; ++i) {
    BOOST_CHECK_EQUAL(labelsFromCuda[i], cpuLabels[i]);
  }
}
