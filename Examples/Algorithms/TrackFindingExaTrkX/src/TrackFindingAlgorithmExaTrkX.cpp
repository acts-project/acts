// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingExaTrkX/TrackFindingAlgorithmExaTrkX.hpp"

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

// TODO do we have these function in the repo somewhere?
float theta(float r, float z) {
  return std::atan2(r, z);
}
float eta(float r, float z) {
  return -std::log(std::tan(0.5 * theta(r, z)));
}

}  // namespace

ActsExamples::TrackFindingAlgorithmExaTrkX::TrackFindingAlgorithmExaTrkX(
    Config config, Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("TrackFindingMLBasedAlgorithm", level),
      m_cfg(std::move(config)),
      m_pipeline(m_cfg.graphConstructor, m_cfg.edgeClassifiers,
                 m_cfg.trackBuilder, logger().clone()) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing spacepoint input collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing protoTrack output collection");
  }

  // Sanitizer run with dummy input to detect configuration issues
  // TODO This would be quite helpful I think, but currently it does not work
  // in general because the stages do not expose the number of node features.
  // However, this must be addressed anyway when we also want to allow to
  // configure this more flexible with e.g. cluster information as input. So
  // for now, we disable this.
#if 0
  if( m_cfg.sanitize ) {
  Eigen::VectorXf dummyInput = Eigen::VectorXf::Random(3 * 15);
  std::vector<float> dummyInputVec(dummyInput.data(),
                                   dummyInput.data() + dummyInput.size());
  std::vector<int> spacepointIDs;
  std::iota(spacepointIDs.begin(), spacepointIDs.end(), 0);

  runPipeline(dummyInputVec, spacepointIDs);
  }
#endif

  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_inputClusters.maybeInitialize(m_cfg.inputClusters);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);

  m_inputTruthGraph.maybeInitialize(m_cfg.inputTruthGraph);
  m_outputGraph.maybeInitialize(m_cfg.outputGraph);

  // reserve space for timing
  m_timing.classifierTimes.resize(
      m_cfg.edgeClassifiers.size(),
      decltype(m_timing.classifierTimes)::value_type{0.f});

  // Check if we want cluster features but do not have them
  const static std::array clFeatures = {
      NodeFeature::eClusterLoc0, NodeFeature::eClusterLoc0,
      NodeFeature::eCellCount,   NodeFeature::eChargeSum,
      NodeFeature::eCluster1R,   NodeFeature::eCluster2R};

  auto wantClFeatures = std::ranges::any_of(
      m_cfg.nodeFeatures,
      [&](const auto& f) { return Acts::rangeContainsValue(clFeatures, f); });

  if (wantClFeatures && !m_inputClusters.isInitialized()) {
    throw std::invalid_argument("Cluster features requested, but not provided");
  }

  if (m_cfg.nodeFeatures.size() != m_cfg.featureScales.size()) {
    throw std::invalid_argument(
        "Number of features mismatches number of scale parameters.");
  }
}

/// Allow access to features with nice names

ActsExamples::ProcessCode ActsExamples::TrackFindingAlgorithmExaTrkX::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  using Clock = std::chrono::high_resolution_clock;
  using Duration = std::chrono::duration<double, std::milli>;
  auto t0 = Clock::now();

  // Setup hooks
  LoopHook hook;

  std::unique_ptr<Acts::TorchTruthGraphMetricsHook> truthGraphHook;
  if (m_inputTruthGraph.isInitialized()) {
    truthGraphHook = std::make_unique<Acts::TorchTruthGraphMetricsHook>(
        m_inputTruthGraph(ctx).edges, this->logger().clone());
    hook.hooks.push_back(&*truthGraphHook);
  }

  std::unique_ptr<Acts::TorchGraphStoreHook> graphStoreHook;
  if (m_outputGraph.isInitialized()) {
    graphStoreHook = std::make_unique<Acts::TorchGraphStoreHook>();
    hook.hooks.push_back(&*graphStoreHook);
  }

  // Read input data
  auto spacepoints = m_inputSpacePoints(ctx);

  std::optional<ClusterContainer> clusters;
  if (m_inputClusters.isInitialized()) {
    clusters = m_inputClusters(ctx);
  }

  // Convert Input data to a list of size [num_measurements x
  // measurement_features]
  const std::size_t numSpacepoints = spacepoints.size();
  const std::size_t numFeatures = m_cfg.nodeFeatures.size();
  ACTS_DEBUG("Received " << numSpacepoints << " spacepoints");
  ACTS_DEBUG("Construct " << numFeatures << " node features");

  std::vector<float> features(numSpacepoints * numFeatures);
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

    // This should be fine, because check in constructor
    Cluster* cl1 = clusters ? &clusters->at(sl1.index()) : nullptr;
    Cluster* cl2 = cl1;

    if (sp.sourceLinks().size() == 2) {
      const auto& sl2 = sp.sourceLinks()[1].template get<IndexSourceLink>();
      cl2 = clusters ? &clusters->at(sl2.index()) : nullptr;
    }

    // I would prefer to use a std::span or boost::span here once available
    float* f = features.data() + isp * numFeatures;

    using NF = NodeFeature;

    auto r = [](const auto& p) {
      return std::hypot(p[Acts::ePos0], p[Acts::ePos1]);
    };
    auto phi = [](const auto& p) {
      return std::atan2(p[Acts::ePos1], p[Acts::ePos0]);
    };

    // clang-format off
#define MAKE_CLUSTER_FEATURES(n) \
    break; case NF::eCluster##n##X:   f[ift] = cl##n->globalPosition[Acts::ePos0]; \
    break; case NF::eCluster##n##Y:   f[ift] = cl##n->globalPosition[Acts::ePos1]; \
    break; case NF::eCluster##n##R:   f[ift] = r(cl##n->globalPosition); \
    break; case NF::eCluster##n##Phi: f[ift] = phi(cl##n->globalPosition); \
    break; case NF::eCluster##n##Z:   f[ift] = cl##n->globalPosition[Acts::ePos2]; \
    break; case NF::eCluster##n##Eta: f[ift] = eta(r(cl##n->globalPosition), cl##n->globalPosition[Acts::ePos2]); \
    break; case NF::eCellCount##n:    f[ift] = cl##n->channels.size(); \
    break; case NF::eChargeSum##n:    f[ift] = cl##n->sumActivations(); \
    break; case NF::eLocDir0##n:      f[ift] = cl##n->localDirection[0]; \
    break; case NF::eLocDir1##n:      f[ift] = cl##n->localDirection[1]; \
    break; case NF::eLocDir2##n:      f[ift] = cl##n->localDirection[2]; \
    break; case NF::eLengthDir0##n:   f[ift] = cl##n->lengthDirection[0]; \
    break; case NF::eLengthDir1##n:   f[ift] = cl##n->lengthDirection[1]; \
    break; case NF::eLengthDir2##n:   f[ift] = cl##n->lengthDirection[2]; \
    break; case NF::eLocEta##n:       f[ift] = cl##n->localEta; \
    break; case NF::eLocPhi##n:       f[ift] = cl##n->localPhi; \
    break; case NF::eGlobEta##n:      f[ift] = cl##n->globalEta; \
    break; case NF::eGlobPhi##n:      f[ift] = cl##n->globalPhi; \
    break; case NF::eEtaAngle##n:     f[ift] = cl##n->etaAngle; \
    break; case NF::ePhiAngle##n:     f[ift] = cl##n->phiAngle;
    // clang-format on

    for (auto ift = 0ul; ift < numFeatures; ++ift) {
      // clang-format off
      switch(m_cfg.nodeFeatures[ift]) {
        // Spacepoint features
        break; case NF::eR:           f[ift] = std::hypot(sp.x(), sp.y());
        break; case NF::ePhi:         f[ift] = std::atan2(sp.y(), sp.x());
        break; case NF::eZ:           f[ift] = sp.z();
        break; case NF::eX:           f[ift] = sp.x();
        break; case NF::eY:           f[ift] = sp.y();
        break; case NF::eEta:         f[ift] = eta(std::hypot(sp.x(), sp.y()), sp.z());
        // Single cluster features
        break; case NF::eClusterLoc0: f[ift] = cl1->sizeLoc0;
        break; case NF::eClusterLoc1: f[ift] = cl1->sizeLoc1;
        break; case NF::eCellCount:   f[ift] = cl1->channels.size();
        break; case NF::eChargeSum:   f[ift] = cl1->sumActivations();
        // Features for split clusters
        MAKE_CLUSTER_FEATURES(1)
        MAKE_CLUSTER_FEATURES(2)
      }
      // clang-format on

      assert(std::isfinite(f[ift]));
      f[ift] /= m_cfg.featureScales[ift];
    }
#undef MAKE_CLUSTER_FEATURES
  }

  auto t1 = Clock::now();

  // Run the pipeline
  const auto trackCandidates = [&]() {
    std::lock_guard<std::mutex> lock(m_mutex);

    Acts::ExaTrkXTiming timing;
    auto res =
        m_pipeline.run(features, moduleIds, spacepointIDs, hook, &timing);

    m_timing.graphBuildingTime(timing.graphBuildingTime.count());

    assert(timing.classifierTimes.size() == m_timing.classifierTimes.size());
    for (auto [aggr, a] :
         Acts::zip(m_timing.classifierTimes, timing.classifierTimes)) {
      aggr(a.count());
    }

    m_timing.trackBuildingTime(timing.trackBuildingTime.count());

    return res;
  }();

  auto t2 = Clock::now();

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

  if (m_outputGraph.isInitialized()) {
    auto graph = graphStoreHook->storedGraph();
    std::transform(
        graph.first.begin(), graph.first.end(), graph.first.begin(),
        [&](const auto& a) -> std::int64_t { return spacepointIDs.at(a); });
    m_outputGraph(ctx, {graph.first, graph.second});
  }

  auto t3 = Clock::now();
  m_timing.preprocessingTime(Duration(t1 - t0).count());
  m_timing.postprocessingTime(Duration(t3 - t2).count());

  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode TrackFindingAlgorithmExaTrkX::finalize() {
  namespace ba = boost::accumulators;

  auto print = [](const auto& t) {
    std::stringstream ss;
    ss << ba::mean(t) << " +- " << std::sqrt(ba::variance(t)) << " ";
    ss << "[" << ba::min(t) << ", " << ba::max(t) << "]";
    return ss.str();
  };

  ACTS_INFO("Exa.TrkX timing info");
  ACTS_INFO("- preprocessing:  " << print(m_timing.preprocessingTime));
  ACTS_INFO("- graph building: " << print(m_timing.graphBuildingTime));
  // clang-format off
  for (const auto& t : m_timing.classifierTimes) {
  ACTS_INFO("- classifier:     " << print(t));
  }
  // clang-format on
  ACTS_INFO("- track building: " << print(m_timing.trackBuildingTime));
  ACTS_INFO("- postprocessing: " << print(m_timing.postprocessingTime));

  return {};
}
