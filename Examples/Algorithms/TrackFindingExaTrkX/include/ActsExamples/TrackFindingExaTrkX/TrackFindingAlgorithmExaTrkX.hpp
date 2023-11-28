// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/ExaTrkX/ExaTrkXPipeline.hpp"
#include "Acts/Plugins/ExaTrkX/Stages.hpp"
#include "Acts/Plugins/ExaTrkX/TorchGraphStoreHook.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <mutex>
#include <string>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

namespace ActsExamples {

class TrackFindingAlgorithmExaTrkX final : public IAlgorithm {
 public:
  struct Config {
    /// Input spacepoints collection.
    std::string inputSpacePoints;

    /// Input cluster information (Optional). If given, the following features
    /// are added:
    /// * cell count
    /// * sum cell activations
    /// * cluster size in local x
    /// * cluster size in local y
    std::string inputClusters;

    /// Input simhits (Optional).
    std::string inputSimHits;
    /// Input measurement simhit map (Optional).
    std::string inputParticles;
    /// Input measurement simhit map (Optional).
    std::string inputMeasurementSimhitsMap;

    /// Output protoTracks collection.
    std::string outputProtoTracks;

    /// Output graph (optional)
    std::string outputGraph;

    std::shared_ptr<Acts::GraphConstructionBase> graphConstructor;

    std::vector<std::shared_ptr<Acts::EdgeClassificationBase>> edgeClassifiers;

    std::shared_ptr<Acts::TrackBuildingBase> trackBuilder;

    /// Scaling of the input features
    float rScale = 1.f;
    float phiScale = 1.f;
    float zScale = 1.f;
    float cellCountScale = 1.f;
    float cellSumScale = 1.f;
    float clusterXScale = 1.f;
    float clusterYScale = 1.f;

    /// Remove track candidates with 2 or less hits
    bool filterShortTracks = false;

    /// Target graph properties
    std::size_t targetMinHits = 3;
    double targetMinPT = 500 * Acts::UnitConstants::MeV;
  };

  /// Constructor of the track finding algorithm
  ///
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  TrackFindingAlgorithmExaTrkX(Config cfg, Acts::Logging::Level lvl);

  ~TrackFindingAlgorithmExaTrkX() override = default;

  /// Framework execute method of the track finding algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final;

  /// Finalize and print timing
  ActsExamples::ProcessCode finalize() final;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  Acts::ExaTrkXPipeline m_pipeline;
  mutable std::mutex m_mutex;

  using Accumulator = boost::accumulators::accumulator_set<
      float, boost::accumulators::features<boost::accumulators::tag::mean,
                                           boost::accumulators::tag::variance>>;

  mutable struct {
    Accumulator graphBuildingTime;
    std::vector<Accumulator> classifierTimes;
    Accumulator trackBuildingTime;
  } m_timing;

  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this,
                                                            "InputSpacePoints"};
  ReadDataHandle<ClusterContainer> m_inputClusters{this, "InputClusters"};

  WriteDataHandle<ProtoTrackContainer> m_outputProtoTracks{this,
                                                           "OutputProtoTracks"};
  WriteDataHandle<Acts::TorchGraphStoreHook::Graph> m_outputGraph{
      this, "OutputGraph"};

  // for truth graph
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<IndexMultimap<Index>> m_inputMeasurementMap{
      this, "InputMeasurementMap"};
};

}  // namespace ActsExamples
