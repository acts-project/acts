// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <Acts/Plugins/ExaTrkX/detail/connectedComponents.cuh>

std::vector<int> cudaConnectedComponents(const std::vector<int> &src,
                                         const std::vector<int> &tgt,
                                         int numNodes) {
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

  connectedComponents<<<1, 1024>>>(cudaSrc, cudaTgt, cudaLabels, cudaLabelsNext,
                                   src.size(), numNodes);

  std::vector<int> labels(numNodes);
  cudaMemcpy(labels.data(), cudaLabels, numNodes * sizeof(int),
             cudaMemcpyDeviceToHost);

  return labels;
}

void check(const std::vector<int> &src, const std::vector<int> &tgt) {
  auto numNodes = std::max(*std::max_element(src.begin(), src.end()),
                           *std::max_element(tgt.begin(), tgt.end())) +
                  1;
  std::cout << "numNodes: " << numNodes << std::endl;

  auto cudaLabels = cudaConnectedComponents(src, tgt, numNodes);

  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
      Graph;
  Graph G(numNodes);

  for (int i = 0; i < src.size(); ++i) {
    boost::add_edge(src[i], tgt[i], G);
  }

  std::vector<std::size_t> cpuLabels(numNodes);
  boost::connected_components(G, &cpuLabels[0]);

  std::cout << "cpu labels:     ";
  for (int i = 0; i < numNodes; ++i) {
    std::cout << cpuLabels[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "my CUDA labels: ";
  for (int i = 0; i < numNodes; ++i) {
    std::cout << cudaLabels[i] << " ";
  }
  std::cout << std::endl;
}

using Vi = std::vector<int>;

BOOST_AUTO_TEST_CASE(simple_test_1) {
  Vi src{0, 1, 2, 3};
  Vi tgt{1, 2, 3, 4};
  check(src, tgt);
}

BOOST_AUTO_TEST_CASE(simple_test_2) {
  Vi src{0, 1, 2, 4, 5, 6};
  Vi tgt{1, 2, 3, 5, 6, 7};
  check(src, tgt);
}

BOOST_AUTO_TEST_CASE(simple_test_3) {
  Vi src{4, 3, 2, 1};
  Vi tgt{3, 2, 1, 0};
  check(src, tgt);
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
  check(src, tgt);
}

#if 0
  // Prefix
  std::vector<int> array{1,2,3,4,5};
  std::vector<int> cpuResult(array.size());
  std::exclusive_scan(array.begin(), array.end(), cpuResult.begin(), 0);

  std::cout << "cpu result: ";
  for(auto c : cpuResult) {
    std::cout << c << " ";
  }
  std::cout << std::endl;

  int *cudaVector;
  cudaMalloc(&cudaVector, array.size()*sizeof(int));
  cudaMemcpy(cudaVector, array.data(), array.size()*sizeof(int), cudaMemcpyHostToDevice);

  thrust::exclusive_scan(thrust::device, cudaVector, cudaVector+array.size(), cudaVector);

  std::cout << "my cuda result: ";

  return 0;
}

#endif
