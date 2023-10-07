// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

namespace ActsExamples {

class SeedsToPrototracks final : public IAlgorithm {
 public:
  struct Config {
    std::string inputSeeds = "seeds";
    std::string outputProtoTracks = "tracks-from-seeds";
  };

  /// Construct the algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SeedsToPrototracks(Config cfg, Acts::Logging::Level lvl);

  /// Run the algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<SimSeedContainer> m_inputSeeds{this, "InputSeeds"};
  WriteDataHandle<ProtoTrackContainer> m_outputProtoTracks{this,
                                                           "OutputProtoTracks"};
};

}  // namespace ActsExamples
