// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/FastJet/Jets.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>

namespace fastjet {
class PseudoJet;
}

namespace HepMC3 {
class GenEvent;
}

namespace ActsFastJet {
class TruthJetBuilder;
class JetProperties;
}  // namespace ActsFastJet

namespace ActsExamples {
struct AlgorithmContext;

class TruthJetAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input particles collection.
    std::string inputTruthParticles;
    /// Output jets collection.
    std::string outputJets;
    /// Minimum jet pT.
    double jetPtMin = 20 * Acts::UnitConstants::GeV;
    std::pair<std::optional<double>, std::optional<double>> jetEtaRange = {
        std::nullopt, std::nullopt};
    double jetClusteringR = 0.4;
    bool clusterHardScatterParticlesOnly = true;

    std::optional<std::string> inputHepMC3Event;
    bool doJetLabeling = true;
    double jetLabelingDeltaR = 0.4;
    double jetLabelingHadronPtMin = 5 * Acts::UnitConstants::GeV;
    bool jetLabelingHardScatterHadronsOnly = true;

    // Isolated TRUTH lepton overlap removal
    bool doOverlapRemoval = false;
    double overlapRemovalDeltaR = 0.2;
    double overlapRemovalIsolationDeltaR = 0.2;
    double overlapRemovalIsolation = 0.1;

    bool debugCsvOutput = false;
  };

  TruthJetAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode initialize() override;

  ProcessCode execute(const AlgorithmContext& ctx) const override;

  const Config& config() const { return m_cfg; }

 private:
  void overlapRemoval(const SimParticleContainer& truthParticles,
                      Acts::FastJet::TrackJetContainer& jets) const;
  Config m_cfg;
  ReadDataHandle<SimParticleContainer> m_inputTruthParticles{
      this, "inputTruthParticles"};
  WriteDataHandle<Acts::FastJet::TrackJetContainer> m_outputJets{this,
                                                                 "outputJets"};

  ReadDataHandle<std::shared_ptr<HepMC3::GenEvent>> m_inputHepMC3Event{
      this, "inputHepMC3Event"};

  // Statistics
  mutable std::atomic<std::size_t> m_numJets = 0;
  mutable std::atomic<std::size_t> m_numJetsAfterOverlapRemoval = 0;
  mutable std::atomic<std::size_t> m_numLightJets = 0;
  mutable std::atomic<std::size_t> m_numCJets = 0;
  mutable std::atomic<std::size_t> m_numBJets = 0;
};

}  // namespace ActsExamples
