// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/ExaTrkX/ExaTrkXPipeline.hpp"
#include "Acts/Plugins/ExaTrkX/Stages.hpp"
#include "Acts/Plugins/ExaTrkX/TorchGraphStoreHook.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Graph.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/TrackFindingExaTrkX/TrackFindingAlgorithmExaTrkX.hpp"
#include "ActsExamples/TrackFindingExaTrkX/TruthGraphBuilder.hpp"

#include <mutex>
#include <string>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

namespace ActsExamples {

class GraphConstructionWithTruthClassifier final : public IAlgorithm {
  using NodeFeature = TrackFindingAlgorithmExaTrkX::NodeFeature;

 public:
  struct Config {
    /// Input spacepoints collection.
    std::string inputSpacePoints;
    std::string inputClusters;

    std::string inputMeasurementParticlesMap;

    /// Output prototracks
    std::string outputProtoTracks;

    /// Graph constructor
    std::shared_ptr<Acts::GraphConstructionBase> graphConstructor;

    /// The track builder
    std::shared_ptr<Acts::TrackBuildingBase> trackBuilder;

    /// Node features
    std::vector<NodeFeature> nodeFeatures = {NodeFeature::eR, NodeFeature::ePhi,
                                             NodeFeature::eZ};

    /// Feature scales
    std::vector<float> featureScales = {1.f, 1.f, 1.f};

    /// Remove track candidates with 2 or less hits
    std::size_t minMeasurementsPerTrack = 3;

    /// Optionally remap the geometry Ids that are put into the chain
    std::shared_ptr<GeometryIdMapActsAthena> geometryIdMap;
  };

  /// Constructor of the track finding algorithm
  ///
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  GraphConstructionWithTruthClassifier(Config cfg, Acts::Logging::Level lvl);

  ~GraphConstructionWithTruthClassifier() override = default;

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

  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this,
                                                            "InputSpacePoints"};
  ReadDataHandle<ClusterContainer> m_inputClusters{this, "InputClusters"};
  ReadDataHandle<IndexMultimap<ActsFatras::Barcode>>
      m_inputMeasurementParticlesMap{this, "InputMeasPartMap"};
  WriteDataHandle<ProtoTrackContainer> m_outputProtoTracks{this,
                                                           "OutputProtoTracks"};
};

}  // namespace ActsExamples
