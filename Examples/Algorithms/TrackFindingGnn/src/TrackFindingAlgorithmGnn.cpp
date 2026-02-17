// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingGnn/TrackFindingAlgorithmGnn.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsPlugins/Gnn/GraphStoreHook.hpp"
#include "ActsPlugins/Gnn/TruthGraphMetricsHook.hpp"
#include "ActsPlugins/Gnn/detail/NvtxUtils.hpp"

#include <algorithm>
#include <chrono>
#include <numeric>

#include "createFeatures.hpp"

#ifdef ACTS_GNN_WITH_CUDA
#include <cuda_runtime_api.h>
#endif

using namespace Acts;
using namespace ActsPlugins;
using namespace Acts::UnitLiterals;

namespace ActsExamples {

namespace {

struct LoopHook : public GnnHook {
  std::vector<GnnHook*> hooks;

  void operator()(const PipelineTensors& tensors,
                  const ExecutionContext& ctx) const override {
    for (auto hook : hooks) {
      (*hook)(tensors, ctx);
    }
  }
};

}  // namespace

TrackFindingAlgorithmGnn::TrackFindingAlgorithmGnn(Config config,
                                                   Logging::Level level)
    : IAlgorithm("TrackFindingMLBasedAlgorithm", level),
      m_cfg(std::move(config)),
      m_pipeline(m_cfg.graphConstructor, m_cfg.edgeClassifiers,
                 m_cfg.trackBuilder, logger().clone()) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing spacepoint input collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing protoTrack output collection");
  }

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
      [&](const auto& f) { return rangeContainsValue(clFeatures, f); });

  if (wantClFeatures && !m_inputClusters.isInitialized()) {
    throw std::invalid_argument("Cluster features requested, but not provided");
  }

  if (m_cfg.nodeFeatures.size() != m_cfg.featureScales.size()) {
    throw std::invalid_argument(
        "Number of features mismatches number of scale parameters.");
  }

// Reset error state to prevent failure due to previous runs
#ifdef ACTS_GNN_WITH_CUDA
  cudaGetLastError();
#endif
}

/// Allow access to features with nice names

ProcessCode TrackFindingAlgorithmGnn::execute(
    const AlgorithmContext& ctx) const {
  ACTS_NVTX_START(data_preparation);

  using Clock = std::chrono::high_resolution_clock;
  using Duration = std::chrono::duration<double, std::milli>;
  auto t0 = Clock::now();

  // Setup hooks
  LoopHook hook;

  std::unique_ptr<TruthGraphMetricsHook> truthGraphHook;
  if (m_inputTruthGraph.isInitialized()) {
    truthGraphHook = std::make_unique<TruthGraphMetricsHook>(
        m_inputTruthGraph(ctx).edges, this->logger().clone());
    hook.hooks.push_back(&*truthGraphHook);
  }

  std::unique_ptr<GraphStoreHook> graphStoreHook;
  if (m_outputGraph.isInitialized()) {
    graphStoreHook = std::make_unique<GraphStoreHook>();
    hook.hooks.push_back(&*graphStoreHook);
  }

  // Read input data
  const auto& spacepoints = m_inputSpacePoints(ctx);

  const ClusterContainer* clusters = nullptr;
  if (m_inputClusters.isInitialized()) {
    clusters = &m_inputClusters(ctx);
  }

  // Convert Input data to a list of size [num_measurements x
  // measurement_features]
  const std::size_t numSpacepoints = spacepoints.size();
  const std::size_t numFeatures = m_cfg.nodeFeatures.size();
  ACTS_DEBUG("Received " << numSpacepoints << " spacepoints");
  ACTS_DEBUG("Construct " << numFeatures << " node features");

  auto t01 = Clock::now();

  std::vector<std::uint64_t> moduleIds;
  moduleIds.reserve(spacepoints.size());

  for (auto isp = 0ul; isp < numSpacepoints; ++isp) {
    const auto& sp = spacepoints[isp];

    // For now just take the first index since does require one single index
    // per spacepoint
    // TODO does it work for the module map construction to use only the first
    // sp?
    const auto& sl1 = sp.sourceLinks().at(0).template get<IndexSourceLink>();

    if (m_cfg.geometryIdMap != nullptr) {
      moduleIds.push_back(m_cfg.geometryIdMap->right.at(sl1.geometryId()));
    } else {
      moduleIds.push_back(sl1.geometryId().value());
    }
  }

  auto t02 = Clock::now();

  // Sort the spacepoints by module ide. Required by module map
  std::vector<int> idxs(numSpacepoints);
  std::iota(idxs.begin(), idxs.end(), 0);
  std::ranges::sort(idxs, {}, [&](auto i) { return moduleIds[i]; });

  std::ranges::sort(moduleIds);

  SimSpacePointContainer sortedSpacepoints;
  sortedSpacepoints.reserve(spacepoints.size());
  std::ranges::transform(idxs, std::back_inserter(sortedSpacepoints),
                         [&](auto i) { return spacepoints[i]; });

  auto t03 = Clock::now();

  auto features = createFeatures(sortedSpacepoints, clusters,
                                 m_cfg.nodeFeatures, m_cfg.featureScales);

  auto t1 = Clock::now();

  auto ms = [](auto a, auto b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
  };
  ACTS_DEBUG("Setup time:              " << ms(t0, t01));
  ACTS_DEBUG("ModuleId mapping & copy: " << ms(t01, t02));
  ACTS_DEBUG("Spacepoint sort:         " << ms(t02, t03));
  ACTS_DEBUG("Feature creation:        " << ms(t03, t1));

  // Run the pipeline
  ACTS_NVTX_STOP(data_preparation);
  GnnTiming timing;
#ifdef ACTS_GNN_CPUONLY
  Device device = {Device::Type::eCPU, 0};
#else
  Device device = {Device::Type::eCUDA, 0};
#endif
  auto trackCandidates =
      m_pipeline.run(features, moduleIds, idxs, device, hook, &timing);
  ACTS_NVTX_START(post_processing);

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
      for (const auto& sl : spacepoints.at(i).sourceLinks()) {
        onetrack.push_back(sl.template get<IndexSourceLink>().index());
      }
    }

    if (onetrack.size() < m_cfg.minMeasurementsPerTrack) {
      nShortTracks++;
      continue;
    }

    protoTracks.push_back(std::move(onetrack));
  }

  ACTS_DEBUG("Removed " << nShortTracks << " with less then "
                        << m_cfg.minMeasurementsPerTrack << " hits");
  ACTS_DEBUG("Created " << protoTracks.size() << " proto tracks");

  m_outputProtoTracks(ctx, std::move(protoTracks));

  if (m_outputGraph.isInitialized()) {
    auto graph = graphStoreHook->storedGraph();
    std::transform(graph.first.begin(), graph.first.end(), graph.first.begin(),
                   [&](const auto& a) -> std::int64_t { return idxs.at(a); });
    m_outputGraph(ctx, {graph.first, graph.second});
  }

  auto t3 = Clock::now();

  {
    std::lock_guard<std::mutex> lock(m_mutex);

    m_timing.preprocessingTime(Duration(t1 - t0).count());
    m_timing.graphBuildingTime(timing.graphBuildingTime.count());

    assert(timing.classifierTimes.size() == m_timing.classifierTimes.size());
    for (auto [aggr, a] :
         zip(m_timing.classifierTimes, timing.classifierTimes)) {
      aggr(a.count());
    }

    m_timing.trackBuildingTime(timing.trackBuildingTime.count());
    m_timing.postprocessingTime(Duration(t3 - t2).count());
    m_timing.fullTime(Duration(t3 - t0).count());
  }

  ACTS_NVTX_STOP(post_processing);
  return ProcessCode::SUCCESS;
}

ProcessCode TrackFindingAlgorithmGnn::finalize() {
  namespace ba = boost::accumulators;

  auto print = [](const auto& t) {
    std::stringstream ss;
    ss << ba::mean(t) << " +- " << std::sqrt(ba::variance(t)) << " ";
    ss << "[" << ba::min(t) << ", " << ba::max(t) << "]";
    return ss.str();
  };

  ACTS_INFO("GNN timing info");
  ACTS_INFO("- preprocessing:  " << print(m_timing.preprocessingTime));
  ACTS_INFO("- graph building: " << print(m_timing.graphBuildingTime));
  // clang-format off
  for (const auto& t : m_timing.classifierTimes) {
  ACTS_INFO("- classifier:     " << print(t));
  }
  // clang-format on
  ACTS_INFO("- track building: " << print(m_timing.trackBuildingTime));
  ACTS_INFO("- postprocessing: " << print(m_timing.postprocessingTime));
  ACTS_INFO("- full timing:    " << print(m_timing.fullTime));

  return {};
}

}  // namespace ActsExamples
