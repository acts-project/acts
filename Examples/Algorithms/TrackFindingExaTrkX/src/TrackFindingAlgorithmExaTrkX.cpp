// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingExaTrkX/TrackFindingAlgorithmExaTrkX.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/ExaTrkX/TorchGraphStoreHook.hpp"
#include "Acts/Plugins/ExaTrkX/TorchTruthGraphMetricsHook.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

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
}

/// Allow access to features with nice names
enum feat : std::size_t {
  eR = 0,
  ePhi,
  eZ,
  eCellCount,
  eCellSum,
  eClusterX,
  eClusterY
};

ActsExamples::ProcessCode ActsExamples::TrackFindingAlgorithmExaTrkX::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Setup hooks
  LoopHook hook;

  std::unique_ptr<Acts::TorchTruthGraphMetricsHook> truthGraphHook;
  if (m_inputTruthGraph.isInitialized()) {
    truthGraphHook = std::make_unique<Acts::TorchTruthGraphMetricsHook>(
        m_inputTruthGraph(ctx), this->logger().clone());
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
  const std::size_t numFeatures = clusters ? 7 : 3;
  ACTS_INFO("Received " << numSpacepoints << " spacepoints");

  std::vector<float> features(numSpacepoints * numFeatures);
  std::vector<int> spacepointIDs;

  spacepointIDs.reserve(spacepoints.size());

  double sumCells = 0.0;
  double sumActivation = 0.0;

  for (auto i = 0ul; i < numSpacepoints; ++i) {
    const auto& sp = spacepoints[i];

    // I would prefer to use a std::span or boost::span here once available
    float* featurePtr = features.data() + i * numFeatures;

    // For now just take the first index since does require one single index
    // per spacepoint
    const auto& sl = sp.sourceLinks()[0].template get<IndexSourceLink>();
    spacepointIDs.push_back(sl.index());

    featurePtr[eR] = std::hypot(sp.x(), sp.y()) / m_cfg.rScale;
    featurePtr[ePhi] = std::atan2(sp.y(), sp.x()) / m_cfg.phiScale;
    featurePtr[eZ] = sp.z() / m_cfg.zScale;

    if (clusters) {
      const auto& cluster = clusters->at(sl.index());
      const auto& chnls = cluster.channels;

      featurePtr[eCellCount] = cluster.channels.size() / m_cfg.cellCountScale;
      featurePtr[eCellSum] =
          std::accumulate(chnls.begin(), chnls.end(), 0.0,
                          [](double s, const Cluster::Cell& c) {
                            return s + c.activation;
                          }) /
          m_cfg.cellSumScale;
      featurePtr[eClusterX] = cluster.sizeLoc0 / m_cfg.clusterXScale;
      featurePtr[eClusterY] = cluster.sizeLoc1 / m_cfg.clusterYScale;

      sumCells += featurePtr[eCellCount];
      sumActivation += featurePtr[eCellSum];
    }
  }

  ACTS_DEBUG("Avg cell count: " << sumCells / spacepoints.size());
  ACTS_DEBUG("Avg activation: " << sumActivation / sumCells);

  // Run the pipeline
  const auto trackCandidates = [&]() {
    std::lock_guard<std::mutex> lock(m_mutex);

    Acts::ExaTrkXTiming timing;
    auto res =
        m_pipeline.run(features, spacepointIDs, hook, &timing);

    m_timing.graphBuildingTime(timing.graphBuildingTime.count());

    assert(timing.classifierTimes.size() == m_timing.classifierTimes.size());
    for (auto [aggr, a] :
         Acts::zip(m_timing.classifierTimes, timing.classifierTimes)) {
      aggr(a.count());
    }

    m_timing.trackBuildingTime(timing.trackBuildingTime.count());

    return res;
  }();

  ACTS_DEBUG("Done with pipeline, received " << trackCandidates.size()
                                             << " candidates");

  // Make the prototracks
  std::vector<ProtoTrack> protoTracks;
  protoTracks.reserve(trackCandidates.size());

  int nShortTracks = 0;

  for (auto& x : trackCandidates) {
    if (m_cfg.filterShortTracks && x.size() < 3) {
      nShortTracks++;
      continue;
    }

    ProtoTrack onetrack;
    onetrack.reserve(x.size());

    std::copy(x.begin(), x.end(), std::back_inserter(onetrack));
    protoTracks.push_back(std::move(onetrack));
  }

  ACTS_INFO("Removed " << nShortTracks << " with less then 3 hits");
  ACTS_INFO("Created " << protoTracks.size() << " proto tracks");
  m_outputProtoTracks(ctx, std::move(protoTracks));

  if (m_outputGraph.isInitialized()) {
    auto graph = graphStoreHook->storedGraph();
    std::transform(
        graph.first.begin(), graph.first.end(), graph.first.begin(),
        [&](const auto& a) -> std::int64_t { return spacepointIDs.at(a); });
    m_outputGraph(ctx, std::move(graph));
  }

  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode TrackFindingAlgorithmExaTrkX::finalize() {
  namespace ba = boost::accumulators;

  ACTS_INFO("Exa.TrkX timing info");
  {
    const auto& t = m_timing.graphBuildingTime;
    ACTS_INFO("- graph building: " << ba::mean(t) << " +- "
                                   << std::sqrt(ba::variance(t)));
  }
  for (const auto& t : m_timing.classifierTimes) {
    ACTS_INFO("- classifier:     " << ba::mean(t) << " +- "
                                   << std::sqrt(ba::variance(t)));
  }
  {
    const auto& t = m_timing.trackBuildingTime;
    ACTS_INFO("- track building: " << ba::mean(t) << " +- "
                                   << std::sqrt(ba::variance(t)));
  }

  return {};
}
