// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingExaTrkX/TrackFindingAlgorithmExaTrkX.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/ExaTrkX/TorchTruthGraphMetricsHook.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

using namespace ActsExamples;
using namespace Acts::UnitLiterals;

namespace {

class ExamplesEdmHook : public Acts::PipelineHook {
  constexpr static double m_targetPT = 0.5_GeV;
  constexpr static std::size_t m_targetSize = 3;

  std::unique_ptr<const Acts::Logger> m_logger;
  std::unique_ptr<Acts::TorchTruthGraphMetricsHook> m_truthGraphHook;
  std::unique_ptr<Acts::TorchTruthGraphMetricsHook> m_targetGraphHook;

  const Acts::Logger& logger() const { return *m_logger; }

  struct SpacepointInfo {
    std::size_t index;
    const ActsFatras::Hit* hit;
  };

 public:
  ExamplesEdmHook(const SimSpacePointContainer& spacepoints,
                  const IndexMultimap<Index>& measHitMap,
                  const SimHitContainer& simhits, const Acts::Logger& logger)
      : m_logger(logger.clone("MetricsHook")) {
    // Associate tracks to graph, collect momentum
    std::unordered_map<ActsFatras::Barcode, std::vector<SpacepointInfo>> tracks;

    for (auto i = 0ul; i < spacepoints.size(); ++i) {
      const auto measId = spacepoints[i]
                              .sourceLinks()[0]
                              .template get<IndexSourceLink>()
                              .index();

      auto [a, b] = measHitMap.equal_range(measId);
      for (auto it = a; it != b; ++it) {
        const auto& hit = *simhits.nth(it->second);

        tracks[hit.particleId()].push_back({i, &hit});
      }
    }

    // Collect edges for truth graph and target graph
    std::vector<int64_t> truthGraph;
    std::vector<int64_t> targetGraph;

    for (auto& [pid, track] : tracks) {
      // Sort by hit index, so the edges are connected correctly
      std::sort(track.begin(), track.end(), [](const auto& a, const auto& b) {
        return a.hit->index() < b.hit->index();
      });

      // Assume that the truth momentum before the first hit is close to the
      // truth momentum of the particles. Avoids pulling in the whole particles
      // collection as well
      const auto& mom4 = track.front().hit->momentum4Before();
      const auto pT = std::hypot(mom4[Acts::eMom0], mom4[Acts::eMom1]);

      for (auto i = 0ul; i < track.size() - 1; ++i) {
        truthGraph.push_back(track[i].index);
        truthGraph.push_back(track[i + 1].index);

        if (pT > m_targetPT && track.size() >= m_targetSize) {
          targetGraph.push_back(track[i].index);
          targetGraph.push_back(track[i + 1].index);
        }
      }
    }

    m_truthGraphHook = std::make_unique<Acts::TorchTruthGraphMetricsHook>(
        truthGraph, logger.clone());
    m_targetGraphHook = std::make_unique<Acts::TorchTruthGraphMetricsHook>(
        targetGraph, logger.clone());
  }

  ~ExamplesEdmHook(){};

  void operator()(const std::any& nodes, const std::any& edges) const override {
    ACTS_INFO("Metrics for total graph:");
    (*m_truthGraphHook)(nodes, edges);
    ACTS_INFO("Metrics for target graph (pT > "
              << m_targetPT / Acts::UnitConstants::GeV
              << " GeV, nHits >= " << m_targetSize << "):");
    (*m_targetGraphHook)(nodes, edges);
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
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);

  m_inputSimHits.maybeInitialize(m_cfg.inputSimHits);
  m_inputMeasurementMap.maybeInitialize(m_cfg.inputMeasurementSimhitsMap);
}

ActsExamples::ProcessCode ActsExamples::TrackFindingAlgorithmExaTrkX::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  auto spacepoints = m_inputSpacePoints(ctx);

  auto hook = std::make_unique<Acts::PipelineHook>();
  if (m_inputSimHits.isInitialized() && m_inputMeasurementMap.isInitialized()) {
    hook = std::make_unique<ExamplesEdmHook>(
        spacepoints, m_inputMeasurementMap(ctx), m_inputSimHits(ctx), logger());
  }

  // Convert Input data to a list of size [num_measurements x
  // measurement_features]
  size_t num_spacepoints = spacepoints.size();
  ACTS_INFO("Received " << num_spacepoints << " spacepoints");

  std::vector<float> inputValues;
  std::vector<int> spacepointIDs;
  inputValues.reserve(spacepoints.size() * 3);
  spacepointIDs.reserve(spacepoints.size());
  for (const auto& sp : spacepoints) {
    float x = sp.x();
    float y = sp.y();
    float z = sp.z();
    float r = sp.r();
    float phi = std::atan2(y, x);

    inputValues.push_back(r / m_cfg.rScale);
    inputValues.push_back(phi / m_cfg.phiScale);
    inputValues.push_back(z / m_cfg.zScale);

    // For now just take the first index since does require one single index per
    // spacepoint
    const auto& islink = sp.sourceLinks()[0].template get<IndexSourceLink>();
    spacepointIDs.push_back(islink.index());
  }

  // Run the pipeline
  const auto trackCandidates =
      m_pipeline.run(inputValues, spacepointIDs, *hook);

  ACTS_DEBUG("Done with pipeline, received " << trackCandidates.size()
                                             << " candidates");

  // Make the prototracks
  std::vector<ProtoTrack> protoTracks;
  protoTracks.reserve(trackCandidates.size());
  for (auto& x : trackCandidates) {
    ProtoTrack onetrack;
    std::copy(x.begin(), x.end(), std::back_inserter(onetrack));
    protoTracks.push_back(std::move(onetrack));
  }

  ACTS_INFO("Created " << protoTracks.size() << " proto tracks");
  m_outputProtoTracks(ctx, std::move(protoTracks));

  return ActsExamples::ProcessCode::SUCCESS;
}
