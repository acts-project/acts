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
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>

namespace fastjet {
class PseudoJet;
}

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
    double jetPtMin;
  };

  TruthJetAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const AlgorithmContext& ctx) const override;
  ProcessCode finalize() override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  ReadDataHandle<SimParticleContainer> m_inputTruthParticles{
      this, "inputTruthParticles"};
  WriteDataHandle<std::vector<fastjet::PseudoJet>> m_outputJets{this,
                                                                "outputJets"};
};

}  // namespace ActsExamples
