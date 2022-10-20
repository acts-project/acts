// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFindingBase.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

class TrackFindingAlgorithmExaTrkX final : public BareAlgorithm {
 public:
  struct Config {
    /// Input spacepoints collection.
    std::string inputSpacePoints;

    /// Output protoTracks collection.
    std::string outputProtoTracks;

    /// ML based track finder
    std::shared_ptr<Acts::ExaTrkXTrackFindingBase> trackFinderML;

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

  virtual ~TrackFindingAlgorithmExaTrkX() {}

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
};

}  // namespace ActsExamples
