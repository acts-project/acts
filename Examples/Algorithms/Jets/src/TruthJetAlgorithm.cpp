// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Jets/TruthJetAlgorithm.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <ostream>
#include <stdexcept>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

namespace ActsExamples {

TruthJetAlgorithm::TruthJetAlgorithm(const Config& cfg,
                                     Acts::Logging::Level lvl)
    : IAlgorithm("TruthJetAlgorithm", lvl), m_cfg(cfg) {
  if (m_cfg.inputTruthParticles.empty()) {
    throw std::invalid_argument("Input particles collection is not configured");
  }
  m_inputTruthParticles.initialize(m_cfg.inputTruthParticles);
  m_outputJets.initialize(m_cfg.outputJets);
}

ProcessCode ActsExamples::TruthJetAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  TrackJetContainer outputJets;

  const auto& truthParticles = m_inputTruthParticles(ctx);

  ACTS_DEBUG("Number of truth particles: " << truthParticles.size());

  const fastjet::JetDefinition defaultJetDefinition =
      fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4);

  // Get the 4-momentum information from the simulated truth particles
  // and create fastjet::PseudoJet objects
  std::vector<fastjet::PseudoJet> inputPseudoJets;

  int particleIndex = 0;
  for (const auto& particle : truthParticles) {
    fastjet::PseudoJet pseudoJet(particle.momentum().x(),
                                 particle.momentum().y(),
                                 particle.momentum().z(), particle.energy());

    pseudoJet.set_user_index(particleIndex);
    inputPseudoJets.push_back(pseudoJet);
    particleIndex++;
  }
  ACTS_DEBUG("Number of input pseudo jets: " << inputPseudoJets.size());

  // Run the jet clustering
  fastjet::ClusterSequence clusterSeq(inputPseudoJets, defaultJetDefinition);

  // Get the jets above a certain pt threshold
  std::vector<fastjet::PseudoJet> jets =
      sorted_by_pt(clusterSeq.inclusive_jets(m_cfg.jetPtMin));
  ACTS_DEBUG("Number of clustered jets: " << jets.size());

  // Prepare jets for the storage - conversion of jets to custom track jet class
  // (and later add here the jet classification)

  for (unsigned int i = 0; i < jets.size(); i++) {
    // Get information on the jet constituents
    std::vector<fastjet::PseudoJet> jetConstituents = jets[i].constituents();
    std::vector<int> constituentIndices;
    constituentIndices.reserve(jetConstituents.size());

    // Get the jet classification label later here! For now, we use "unknown"
    ActsExamples::jetlabel label = ActsExamples::unknown;
    
    Acts::Vector4 jetFourMomentum(jets[i].px(), jets[i].py(), jets[i].pz(),
                                  jets[i].e());

    // Initialize the (track) jet with 4-momentum and jet label
    ActsExamples::TrackJet storedJet(jetFourMomentum, label);

    // Add the jet constituents to the (track)jet
    for (unsigned int j = 0; j < jetConstituents.size(); j++) {
      // Get the index of the constituent in the original input pseudo jets
      constituentIndices.push_back(jetConstituents[j].user_index());
    }

    storedJet.setConstituents(constituentIndices);

    outputJets.push_back(storedJet);
    ACTS_DEBUG("Stored jet " << i << " with 4-momentum: " << jetFourMomentum(0)
                             << ", " << jetFourMomentum(1) << ", "
                             << jetFourMomentum(2) << ", " << jetFourMomentum(3)
                             << " and " << constituentIndices.size()
                             << " constituents.");
  }

  m_outputJets(ctx, std::move(outputJets));
  return ProcessCode::SUCCESS;
}

ProcessCode ActsExamples::TruthJetAlgorithm::finalize() {
  ACTS_INFO("Finalizing truth jet clustering");
  return ProcessCode::SUCCESS;
}

};  // namespace ActsExamples
