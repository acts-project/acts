// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingExaTrkX/TrackFindingAlgorithmExaTrkX.hpp"

#include "Acts/Plugins/ExaTrkX/TorchTruthGraphMetricsHook.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <torch/torch.h>

using namespace ActsExamples;

namespace {

auto makeHook(const SimSpacePointContainer& spacepoints,
              const IndexMultimap<Index>& measHitMap,
              const SimHitContainer& truthHits, const Acts::Logger& logger) {
  std::unordered_map<SimBarcode, std::pair<double, std::vector<std::size_t>>>
      tracks;

  for (auto i = 0ul; i < spacepoints.size(); ++i) {
    const auto measId =
        spacepoints[i].sourceLinks()[0].template get<IndexSourceLink>().index();

    auto [a, b] = measHitMap.equal_range(measId);
    for (auto it = a; it != b; ++it) {
      const auto& hit = *truthHits.nth(it->second);
      const auto pid = hit.particleId();
      tracks[pid].first = std::max(tracks[pid].first, hit.momentum4Before()[3]);
      tracks[pid].second.push_back(i);
    }
  }

  // for(const auto &[pid, track] : tracks) {
  //     const auto &[m, t] = track;
  //     ACTS_VERBOSE(pid << ": " << t);
  // }
  // for(auto it=tracks.begin(); it != tracks.end();) {
  //     if( (it->second.first < 0.5) && (it->second.second.size() < 3) ) {
  //         it = tracks.erase(it);
  //     } else {
  //         ++it;
  //     }
  // }

  std::vector<int64_t> truthGraph;
  for (auto& [_, v] : tracks) {
    auto& [mom, track] = v;
    std::sort(track.begin(), track.end());
    for (auto i = 0ul; i < track.size() - 1; ++i) {
      truthGraph.push_back(track[i]);
      truthGraph.push_back(track[i + 1]);
    }
  }

  return std::make_unique<Acts::TorchTruthGraphMetricsHook>(truthGraph,
                                                            logger.clone());
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
  // However, this must be addressed anyways when we also want to allow to
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

  m_inputSimHits.maybeInitialize(m_cfg.inputSimhits);
  m_inputMeasurementMap.maybeInitialize(m_cfg.inputMeasurementSimhitMap);
}

enum feat : std::size_t {
  eR,
  ePhi,
  eZ,
  eCellCount,
  eCellSum,
  eClusterX,
  eClusterY
};

ActsExamples::ProcessCode ActsExamples::TrackFindingAlgorithmExaTrkX::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  auto spacepoints = m_inputSpacePoints(ctx);
  std::sort(spacepoints.begin(), spacepoints.end(),
            [](const auto& a, const auto& b) {
              return std::hypot(a.x(), a.y()) < std::hypot(b.x(), b.y());
            });

  std::optional<ClusterContainer> clusters;
  if (m_inputClusters.isInitialized()) {
    clusters = m_inputClusters(ctx);
  }

  auto hook = std::make_unique<Acts::PipelineHook>();
  if (m_inputSimHits.isInitialized() && m_inputMeasurementMap.isInitialized()) {
    hook = makeHook(spacepoints, m_inputMeasurementMap(ctx),
                    m_inputSimHits(ctx), logger());
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

  ACTS_INFO("First spacepoint x,y,z=" << spacepoints[0].x() << ", "
                                      << spacepoints[0].y() << ", "
                                      << spacepoints[0].z());

  for (auto i = 0ul; i < numSpacepoints; ++i) {
    const auto& sp = spacepoints[i];

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
  const auto trackCandidates = m_pipeline.run(features, spacepointIDs, *hook);

  ACTS_DEBUG("Done with pipeline");

  // Make the prototracks
  std::vector<ProtoTrack> protoTracks;
  protoTracks.reserve(trackCandidates.size());

  int nShortTracks = 0;

  for (auto& x : trackCandidates) {
    ProtoTrack onetrack;
    std::copy(x.begin(), x.end(), std::back_inserter(onetrack));

    if (onetrack.size() < 3) {
      nShortTracks++;
      continue;
    }
    protoTracks.push_back(std::move(onetrack));
  }

  ACTS_INFO("Removed " << nShortTracks << " with less then 3 hits");
  ACTS_INFO("Created " << protoTracks.size() << " proto tracks");
  m_outputProtoTracks(ctx, std::move(protoTracks));

  return ActsExamples::ProcessCode::SUCCESS;
}
