// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsPlugins/Gnn/TruthGraphMetricsHook.hpp"

#include <cassert>

using namespace Acts;
using namespace ActsPlugins;

void testTruthTestGraph(std::vector<std::int64_t> &truthGraph,
                        std::vector<std::int64_t> &testGraph,
                        const std::string &resStr) {
  std::stringstream ss;
  auto logger = getDefaultLogger("Test", Logging::INFO, &ss);

  TruthGraphMetricsHook hook(truthGraph, std::move(logger));

  auto numTestEdges = testGraph.size() / 2;
  auto edgeIndexTensor =
      Tensor<std::int64_t>::Create({2, numTestEdges}, {Device::Cpu(), {}});

  // Transpose the input vector into the tensor
  for (auto i = 0ul; i < numTestEdges; ++i) {
    *(edgeIndexTensor.data() + i) = testGraph.at(2 * i);
    *(edgeIndexTensor.data() + numTestEdges + i) = testGraph.at(2 * i + 1);
  }

  PipelineTensors tensors{Tensor<float>::Create({1, 1}, {Device::Cpu(), {}}),
                          std::move(edgeIndexTensor),
                          {},
                          {}};

  hook(tensors, {Device::Cpu(), {}});

  const auto str = ss.str();

  auto begin = str.begin() + str.find("Efficiency");
  BOOST_CHECK_EQUAL(std::string(begin, str.end() - 1), resStr);
}

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(GnnSuite)

BOOST_AUTO_TEST_CASE(same_graph) {
  // clang-format off
  std::vector<std::int64_t> truthGraph = {
    1,2,
    2,3,
    3,4,
    4,5,
  };
  // clang-format on

  // clang-format off
  std::vector<std::int64_t> testGraph = {
    3,4,
    4,5,
    1,2,
    2,3,
  };
  // clang-format on

  testTruthTestGraph(truthGraph, testGraph, "Efficiency=1, purity=1");
}

// Test large numbers because overflows are a danger for cantor pairing
BOOST_AUTO_TEST_CASE(same_graph_large_numbers) {
  // clang-format off
  std::int64_t k = 100'000;

  std::vector<std::int64_t> truthGraph = {
    1,2,
    2,3,
    3,4,
    4,5,
  };
  // clang-format on
  std::transform(truthGraph.begin(), truthGraph.end(), truthGraph.begin(),
                 [&](auto i) { return k + i; });

  // clang-format off
  std::vector<std::int64_t> testGraph = {
    3,4,
    4,5,
    1,2,
    2,3,
  };
  // clang-format on
  std::transform(testGraph.begin(), testGraph.end(), testGraph.begin(),
                 [&](auto i) { return k + i; });

  testTruthTestGraph(truthGraph, testGraph, "Efficiency=1, purity=1");
}

BOOST_AUTO_TEST_CASE(fifty_fifty) {
  // clang-format off
  std::vector<std::int64_t> truthGraph = {
    1,2,
    2,3,
    3,4,
    4,5,
  };
  // clang-format on

  // clang-format off
  std::vector<std::int64_t> testGraph = {
    3,4,
    4,5,
    6,9,
    5,1,
  };
  // clang-format on

  testTruthTestGraph(truthGraph, testGraph, "Efficiency=0.5, purity=0.5");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
