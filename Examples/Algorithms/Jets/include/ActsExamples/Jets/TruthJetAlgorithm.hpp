// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsPlugins/FastJet/Jets.hpp"

#include <string>

namespace ActsExamples {

using TruthJetContainer = std::vector<ActsPlugins::FastJet::TruthJet>;

class TruthJetAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input particles collection.
    std::string inputTruthParticles;
    /// Input tracks collection.
    std::string inputTracks;
    /// Output jets collection.
    std::string outputJets;
    /// Minimum jet pT.
    double jetPtMin = 20 * Acts::UnitConstants::GeV;
    /// Jet eta range a pair of doubles defaulted to -inf/+inf
    std::pair<double, double> jetEtaRange = {
        -std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    /// Jet clustering radius
    double jetClusteringRadius = 0.4;
    /// Only cluster HS particles
    bool clusterHSParticlesOnly = true;
    /// Do jet labeling
    bool doJetLabeling = true;
    /// Delta R for labeling
    double jetLabelingDeltaR = 0.4;
    /// Minimum hadron pT for labeling
    double jetLabelingHadronPtMin = 5 * Acts::UnitConstants::GeV;
    /// Only label HS hadrons
    bool jetLabelingHSHadronsOnly = true;
    /// Enable track-jet matching
    bool doTrackJetMatching = false;
  };

  TruthJetAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const AlgorithmContext& ctx) const override;
  ProcessCode finalize() override;

  void trackJetMatching(const ConstTrackContainer& tracks,
                        TruthJetContainer& jets) const;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  ReadDataHandle<SimParticleContainer> m_inputTruthParticles{
      this, "inputTruthParticles"};
  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "inputTracks"};
  WriteDataHandle<TruthJetContainer> m_outputJets{this, "outputJets"};
};

}  // namespace ActsExamples
