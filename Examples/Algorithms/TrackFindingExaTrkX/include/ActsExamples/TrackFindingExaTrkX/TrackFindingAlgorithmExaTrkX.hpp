// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/Pipeline.hpp"
#include "Acts/Plugins/ExaTrkX/Stages.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

class TrackFindingAlgorithmExaTrkX final : public IAlgorithm {
 public:
  struct Config {
    /// Input spacepoints collection.
    std::string inputSpacePoints;

    /// Input simhits (Optional).
    std::string inputSimHits;
    /// Input measurement simhit map (Optional).
    std::string inputParticles;
    /// Input measurement simhit map (Optional).
    std::string inputMeasurementSimhitsMap;

    /// Output protoTracks collection.
    std::string outputProtoTracks;

    std::shared_ptr<Acts::GraphConstructionBase> graphConstructor;

    std::vector<std::shared_ptr<Acts::EdgeClassificationBase>> edgeClassifiers;

    std::shared_ptr<Acts::TrackBuildingBase> trackBuilder;

    /// Scaling of the input features
    float rScale = 1.f;
    float phiScale = 1.f;
    float zScale = 1.f;
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

  const Config& config() const { return m_cfg; }

 private:
  // configuration
  Config m_cfg;

  Acts::Pipeline m_pipeline;

  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this,
                                                            "InputSpacePoints"};
  WriteDataHandle<ProtoTrackContainer> m_outputProtoTracks{this,
                                                           "OutputProtoTracks"};

  // for truth graph
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<IndexMultimap<Index>> m_inputMeasurementMap{
      this, "InputMeasurementMap"};
};

}  // namespace ActsExamples
