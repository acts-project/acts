// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingExaTrkX/GraphConstructionWithTruthClassifier.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/ExaTrkX/TorchGraphStoreHook.hpp"
#include "Acts/Plugins/ExaTrkX/TorchTruthGraphMetricsHook.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <algorithm>
#include <chrono>
#include <numeric>
#include <ranges>

#include <torch/torch.h>

#include "createFeatures.hpp"

using namespace ActsExamples;
using namespace Acts::UnitLiterals;

namespace {

struct LoopHook : public Acts::ExaTrkXHook {
  std::vector<Acts::ExaTrkXHook*> hooks;

  ~LoopHook() {}

  void operator()(const std::any& nodes, const std::any& edges,
                  const std::any& weights) const override {
    for (auto hook : hooks) {
      (*hook)(nodes, edges, weights);
    }
  }
};

}  // namespace

ActsExamples::GraphConstructionWithTruthClassifier::
    GraphConstructionWithTruthClassifier(Config config,
                                         Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("GraphConstructionWithTruthClassifier", level),
      m_cfg(std::move(config)) {
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_inputClusters.maybeInitialize(m_cfg.inputClusters);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
}

/// Allow access to features with nice names

ActsExamples::ProcessCode
ActsExamples::GraphConstructionWithTruthClassifier::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  auto spacepoints = m_inputSpacePoints(ctx);
  auto measParticlesMap = m_inputMeasurementParticlesMap(ctx);
  auto clusters = m_inputClusters(ctx);

  // Build spacepoint particles map
  std::multimap<std::size_t, ActsFatras::Barcode> spacepointParticlesMap;
  for (std::size_t isp = 0; isp < spacepoints.size(); ++isp) {
    const auto& sp = spacepoints.at(isp);
    for (auto sl : sp.sourceLinks()) {
      auto measIdx = sl.template get<IndexSourceLink>().index();

      auto [b, e] = measParticlesMap.equal_range(measIdx);
      for (auto it = b; it != e; ++it) {
        spacepointParticlesMap.emplace(isp, it->second);
      }
    }
  }

  // Convert Input data to a list of size [num_measurements x
  // measurement_features]
  const std::size_t numSpacepoints = spacepoints.size();
  const std::size_t numFeatures = m_cfg.nodeFeatures.size();
  ACTS_DEBUG("Received " << numSpacepoints << " spacepoints");
  ACTS_DEBUG("Construct " << numFeatures << " node features");

  std::vector<int> spacepointIDs;
  std::vector<std::uint64_t> moduleIds;

  spacepointIDs.reserve(spacepoints.size());
  moduleIds.reserve(spacepoints.size());

  for (auto isp = 0ul; isp < numSpacepoints; ++isp) {
    const auto& sp = spacepoints[isp];

    // For now just take the first index since does require one single index
    // per spacepoint
    // TODO does it work for the module map construction to use only the first
    // sp?
    const auto& sl1 = sp.sourceLinks()[0].template get<IndexSourceLink>();

    // TODO this makes it a bit useless, refactor so we do not need to pass this
    // to the pipeline
    spacepointIDs.push_back(isp);

    if (m_cfg.geometryIdMap != nullptr) {
      moduleIds.push_back(m_cfg.geometryIdMap->right.at(sl1.geometryId()));
    } else {
      moduleIds.push_back(sl1.geometryId().value());
    }
  }

  auto features = createFeatures(spacepoints, clusters, m_cfg.nodeFeatures,
                                 m_cfg.featureScales);

  auto [nodeFeaturesAny, edgeIndexAny, edgeFeaturesAny] =
      (*m_cfg.graphConstructor)(features, spacepointIDs.size(), moduleIds,
                                m_cfg.graphConstructor->device());

  auto edgeIndex = std::any_cast<at::Tensor>(edgeIndexAny);

  using namespace torch::indexing;
  edgeIndex = edgeIndex.to(torch::kCPU);

  auto scores = torch::zeros(edgeIndex.size(1));
  auto scoreAccessor = scores.accessor<float, 1>();

  auto srcNodes = edgeIndex.index({0, Slice()});
  auto srcNodeAccessor = srcNodes.accessor<std::int64_t, 1>();

  auto dstNodes = edgeIndex.index({1, Slice()});
  auto dstNodeAccessor = dstNodes.accessor<std::int64_t, 1>();

  std::vector<std::uint64_t> srcParticles, dstParticles, intersection;
  for (auto i = 0; i < edgeIndex.size(1); ++i) {
    srcParticles.clear();
    dstParticles.clear();
    intersection.clear();

    auto src = srcNodeAccessor[i];
    auto dst = dstNodeAccessor[i];

    auto [srcA, srcB] = spacepointParticlesMap.equal_range(src);
    std::transform(srcA, srcB, std::back_inserter(srcParticles),
                   [](auto& a) { return a.second.value(); });
    std::ranges::sort(srcParticles);

    auto [dstA, dstB] = spacepointParticlesMap.equal_range(dst);
    std::transform(dstA, dstB, std::back_inserter(dstParticles),
                   [](auto& a) { return a.second.value(); });
    std::ranges::sort(dstParticles);

    std::ranges::set_intersection(srcParticles, dstParticles,
                                  std::back_inserter(intersection));

    if (!intersection.empty()) {
      scoreAccessor[i] = 1.f;
    }
  }

  torch::Tensor mask = scores > 0.5f;
  std::any edgeIndexAfterCut = edgeIndex.index({Slice(), mask});
  scores = scores.index({mask});

  auto trackCandidates = (*m_cfg.trackBuilder)(
      std::move(nodeFeaturesAny), std::move(edgeIndexAfterCut), scores,
      spacepointIDs, m_cfg.trackBuilder->device());

  ACTS_DEBUG("Done with pipeline, received " << trackCandidates.size()
                                             << " candidates");

  // Make the prototracks
  std::vector<ProtoTrack> protoTracks;
  protoTracks.reserve(trackCandidates.size());

  int nShortTracks = 0;

  /// TODO the whole conversion back to meas idxs should be pulled out of the
  /// track trackBuilder
  for (auto& candidate : trackCandidates) {
    ProtoTrack onetrack;
    onetrack.reserve(candidate.size());

    for (auto i : candidate) {
      for (const auto& sl : spacepoints[i].sourceLinks()) {
        onetrack.push_back(sl.template get<IndexSourceLink>().index());
      }
    }

    if (onetrack.size() < m_cfg.minMeasurementsPerTrack) {
      nShortTracks++;
      continue;
    }

    protoTracks.push_back(std::move(onetrack));
  }

  ACTS_INFO("Removed " << nShortTracks << " with less then "
                       << m_cfg.minMeasurementsPerTrack << " hits");
  ACTS_INFO("Created " << protoTracks.size() << " proto tracks");
  m_outputProtoTracks(ctx, std::move(protoTracks));

  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode GraphConstructionWithTruthClassifier::finalize() {
  return {};
}
