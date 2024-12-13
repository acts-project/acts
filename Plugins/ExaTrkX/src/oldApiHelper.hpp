// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/detail/TensorVectorConversion.hpp"
#include "Acts/Plugins/ExaTrkX/detail/Utils.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <graph_creator>
#include <vector>

#include <torch/torch.h>

inline TTree_hits<float> makeTTreeHits(
    const std::vector<float> &inputValues,
    const std::vector<std::uint64_t> &moduleIds, float rScale, float phiScale,
    float zScale, const Acts::Logger &logger = Acts::getDummyLogger()) {
  const auto numNodes = moduleIds.size();
  const auto numFeatures = inputValues.size() / numNodes;

  hits<float> hitsCollection(false, false);

  ACTS_DEBUG("Start collecting hits...");

  for (auto i = 0ul; i < numNodes; ++i) {
    // TODO Use std::span when we move to C++20
    const float *hitFeatures = inputValues.data() + i * numFeatures;

    int hitId = static_cast<int>(i);

    // Needs to be rescaled because ModuleMapGraph expects unscaled features
    float r = hitFeatures[0] * rScale;
    float phi = hitFeatures[1] * phiScale;
    float z = hitFeatures[2] * zScale;

    float x = r * std::cos(phi);
    float y = r * std::sin(phi);

    std::uint64_t particleId = 0;  // We do not know
    std::uint64_t moduleId = moduleIds[i];
    std::string hardware = "";      // now hardware
    int barrelEndcap = 0;           // unclear, is this a flag???
    std::uint64_t particleID1 = 0;  // unclear
    std::uint64_t particleID2 = 0;  // unclear

    hit<float> hit(hitId, x, y, z, particleId, moduleId, hardware, barrelEndcap,
                   particleID1, particleID2);

    hitsCollection += hit;
  }

  TTree_hits<float> hitsTree = hitsCollection;
  return hitsTree;
}

template <typename Builder>
std::pair<at::Tensor, at::Tensor> oldApiBuild(
    const std::vector<float> &inputValues,
    const std::vector<std::uint64_t> &moduleIds, const Acts::Logger &logger,
    const Builder &builder, float rScale, float phiScale, float zScale) {
  using namespace torch::indexing;

  const auto numNodes = moduleIds.size();
  const auto numFeatures = inputValues.size() / numNodes;

  const auto t0 = std::chrono::high_resolution_clock::now();
  auto hitsTree =
      makeTTreeHits(inputValues, moduleIds, rScale, phiScale, zScale, logger);

  ACTS_DEBUG("Hits tree has " << hitsTree.size()
                              << " hits, now build graph...");
  const auto t1 = std::chrono::high_resolution_clock::now();
  auto graph = builder(hitsTree, logger().level() <= Acts::Logging::DEBUG);
  const auto t2 = std::chrono::high_resolution_clock::now();
  const auto numEdges = boost::num_edges(graph.graph_impl());

  ACTS_DEBUG("Got " << numEdges << " edges, put them in a vector...");
  std::vector<std::int64_t> edgeIndexVector;
  edgeIndexVector.reserve(2 * numEdges);

  constexpr std::size_t numEdgeFeatures = 6;
  std::vector<float> edgeFeatureVector;
  edgeFeatureVector.reserve(numEdgeFeatures * numEdges);

  // TODO I think this is already somewhere in the codebase
  const float pi = static_cast<float>(3.141592654);
  auto resetAngle = [pi](float angle) {
    if (angle > pi) {
      return angle - 2.f * pi;
    }
    if (angle < -pi) {
      return angle + 2.f * pi;
    }
    return angle;
  };

  std::ofstream of("old_api_edges.csv");
  of << "src,tgt\n";

  auto [begin, end] = boost::edges(graph.graph_impl());
  for (auto it = begin; it != end; ++it) {
    const auto &edge = *it;

    auto src = graph.graph_impl()[boost::source(edge, graph.graph_impl())];
    auto dst = graph.graph_impl()[boost::target(edge, graph.graph_impl())];

    // Edge index
    assert(src.hit_id() >= 0 && src.hit_id() < static_cast<int>(numNodes));
    assert(dst.hit_id() >= 0 && dst.hit_id() < static_cast<int>(numNodes));

    edgeIndexVector.push_back(src.hit_id());
    edgeIndexVector.push_back(dst.hit_id());
    of << src.hit_id() << "," << dst.hit_id() << "\n";

    // Edge features
    // See
    // https://gitlab.cern.ch/gnn4itkteam/acorn/-/blob/dev/acorn/utils/loading_utils.py?ref_type=heads#L288
    const float *srcFeatures = inputValues.data() + src.hit_id() * numFeatures;
    const float *dstFeatures = inputValues.data() + dst.hit_id() * numFeatures;

    const float deltaR = dstFeatures[0] - srcFeatures[0];
    const float deltaPhi =
        resetAngle((dstFeatures[1] - srcFeatures[1]) * phiScale) / phiScale;
    const float deltaZ = dstFeatures[2] - srcFeatures[2];
    const float deltaEta = dstFeatures[3] - srcFeatures[3];

    // In acorn, nan gets converted to 0
    float phiSlope = 0.f;
    float rPhiSlope = 0.f;

    if (deltaR != 0) {
      phiSlope = std::clamp(deltaPhi / deltaR, -100.f, 100.f);

      const float avgR = 0.5f * (dstFeatures[0] + srcFeatures[0]);
      rPhiSlope = avgR * phiSlope;
    }

    for (auto f : {deltaR, deltaPhi, deltaZ, deltaEta, phiSlope, rPhiSlope}) {
      edgeFeatureVector.push_back(f);
    }
  }

  const auto t3 = std::chrono::high_resolution_clock::now();
  // Build final tensors
  ACTS_DEBUG("Construct final tensors...");
  assert(edgeIndexVector.empty() || edgeIndexVector.size() % numEdges == 0);
  auto edgeIndex = Acts::detail::vectorToTensor2D(edgeIndexVector, 2)
                       .t()
                       .contiguous()
                       .to(torch::kInt64);

  assert(edgeIndexVector.empty() ||
         edgeFeatureVector.size() % numEdgeFeatures == 0);
  auto edgeFeatures =
      Acts::detail::vectorToTensor2D(edgeFeatureVector, numEdgeFeatures)
          .clone();

  ACTS_DEBUG("edgeIndex: " << Acts::detail::TensorDetails{edgeIndex});
  ACTS_DEBUG("edgeFeatures: " << Acts::detail::TensorDetails{edgeFeatures});

  ACTS_VERBOSE("Edge vector:\n"
               << (Acts::detail::RangePrinter{edgeIndexVector.begin(),
                                              edgeIndexVector.begin() + 10}));
  ACTS_VERBOSE("Edge index slice:\n"
               << edgeIndex.index({Slice(0, 2), Slice(0, 9)}));

  ACTS_VERBOSE("Edge features slice:\n"
               << edgeFeatures.index({Slice(0, 9), Slice()}));

  const auto t4 = std::chrono::high_resolution_clock::now();
  auto count_ms = [](auto ta, auto tb) {
    return std::chrono::duration<double, std::milli>(tb - ta).count();
  };
  ACTS_DEBUG("Preparation: " << count_ms(t0, t1) << " ms");
  ACTS_DEBUG("Build MM: " << count_ms(t1, t2) << " ms");
  ACTS_DEBUG("Edge features: " << count_ms(t2, t3) << " ms");
  ACTS_DEBUG("Tensors: " << count_ms(t3, t4) << " ms");

  return {edgeIndex, edgeFeatures};
}
