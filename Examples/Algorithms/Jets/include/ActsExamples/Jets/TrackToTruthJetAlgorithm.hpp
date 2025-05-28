// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TrackJet.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

namespace ActsExamples {
struct AlgorithmContext;

/// Print all particles.
class TrackToTruthJetAlgorithm : public IAlgorithm {
 public:
  struct Config {
    /// Input tracks collection.
    std::string inputTracks;
    /// Input jets collection.
    std::string inputJets;
    /// Output track jets collection.
    std::string outputTrackJets;
    /// Maximum delta R for track to jet matching.
    double maxDeltaR = 0.4;
  };

  TrackToTruthJetAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const AlgorithmContext& ctx) const override;
  ProcessCode finalize() const;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "inputTracks"};
  ReadDataHandle<TrackJetContainer> m_inputJets{this,
                                                              "inputJets"};
  WriteDataHandle<TrackJetContainer> m_outputTrackJets{this, "outputTrackJets"};
};

}  // namespace ActsExamples
