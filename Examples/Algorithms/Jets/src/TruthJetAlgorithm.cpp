// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Jets/TruthJetAlgorithm.hpp"

#include "Acts/Definitions/ParticleData.hpp"
#include "ActsExamples/Utilities/ParticleId.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "Acts/Plugins/FastJet/Jets.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Util.hpp"

#include <algorithm>
#include <fstream>
#include <mutex>
#include <ostream>
#include <ranges>
#include <stdexcept>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/Print.h>
#include <boost/container/flat_map.hpp>
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

  if (m_cfg.doJetLabeling) {
    if (!m_cfg.inputHepMC3Event) {
      throw std::invalid_argument(
          "Input HepMC3 event is not configured for jet labeling");
    }
    m_inputHepMC3Event.initialize(m_cfg.inputHepMC3Event.value());
  }
}

namespace {

Acts::FastJet::JetLabel jetLabelFromHadronType(ActsExamples::ParticleId::HadronType hType) {
  using enum ActsExamples::ParticleId::HadronType;
  switch (hType) {
    case BBbarMeson:
    case BottomMeson:
    case BottomBaryon:
      return Acts::FastJet::JetLabel::BJet;
    case CCbarMeson:
    case CharmedMeson:
    case CharmedBaryon:
      return Acts::FastJet::JetLabel::CJet;
    case StrangeMeson:
    case StrangeBaryon:
    case LightMeson:
    case LightBaryon:
      return Acts::FastJet::JetLabel::LightJet;
    default:
      return Acts::FastJet::JetLabel::Unknown;
  }
}
}  // namespace

ProcessCode ActsExamples::TruthJetAlgorithm::initialize() {
  if (m_cfg.debugCsvOutput) {
    std::ofstream outfile;
    outfile.open("particles.csv");
    outfile << "event,pt,eta,phi,pdg,label" << std::endl;

    outfile.flush();
    outfile.close();

    outfile.open("jets.csv");
    outfile << "event,pt,eta,phi,label" << std::endl;

    outfile.flush();
    outfile.close();

    outfile.open("hadrons.csv");
    outfile << "event,pt,eta,phi,pdg,label" << std::endl;

    outfile.flush();
    outfile.close();
  }

  return ProcessCode::SUCCESS;
}


ProcessCode ActsExamples::TruthJetAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  Acts::FastJet::TrackJetContainer outputJets;

  const auto& truthParticles = m_inputTruthParticles(ctx);

  ACTS_DEBUG("Number of truth particles: " << truthParticles.size());

  const fastjet::JetDefinition defaultJetDefinition =
      fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4);

  // Get the 4-momentum information from the simulated truth particles
  // and create fastjet::PseudoJet objects
  std::vector<fastjet::PseudoJet> inputPseudoJets;

  std::vector<SimParticle> inputParticles;
  std::ranges::copy(truthParticles, std::back_inserter(inputParticles));

  for (unsigned int i = 0; i < inputParticles.size(); i++) {
    const auto& particle = inputParticles.at(i);
    fastjet::PseudoJet pseudoJet(particle.momentum().x(),
                                 particle.momentum().y(),
                                 particle.momentum().z(), particle.energy());

    pseudoJet.set_user_index(i);
    inputPseudoJets.push_back(pseudoJet);
  }
  ACTS_DEBUG("Number of input pseudo jets from truth particles: "
             << inputPseudoJets.size());

  // Run the jet clustering
  fastjet::ClusterSequence clusterSeq(inputPseudoJets, defaultJetDefinition);

  // Get the jets above a certain pt threshold
  std::vector<fastjet::PseudoJet> jets =
      sorted_by_pt(clusterSeq.inclusive_jets(m_cfg.jetPtMin));
  ACTS_DEBUG("Number of clustered truth jets: " << jets.size());

  std::vector<std::pair<const HepMC3::GenParticle*, ParticleId::HadronType>>
      hadrons;
  if (m_cfg.doJetLabeling) {
    ACTS_DEBUG("Jet labeling is enabled");
    const auto& genEvent = *m_inputHepMC3Event(ctx);

    for (const auto& particle : genEvent.particles()) {
      if (!ParticleId::isHadron(particle->pdg_id())) {
        continue;
      }
      // ACTS_VERBOSE("Particle " << particle->pdg_id() << " with status "
      // HepMC3::Print::line(particle);

      hadrons.emplace_back(particle.get(),
                           ParticleId::hadronLabel(particle->pdg_id()));
      // auto hadronLabel = ParticleId::hadronLabel(particle->pdg_id());
      // std::cout << "Hadron label: " << hadronLabel << std::endl;
    }
  }

  // Prepare jets for the storage - conversion of jets to custom track jet class
  // (and later add here the jet classification)

  auto deltaR = [](const auto& a, const auto& b) {
    double dphi = abs(a.phi() - b.phi());
    if (dphi > std::numbers::pi) {
      dphi = std::numbers::pi * 2 - dphi;
    }
    double drap = a.rap() - b.rap();
    return std::sqrt(dphi * dphi + drap * drap);
  };

  for (unsigned int i = 0; i < jets.size(); i++) {
    // Get information on the jet constituents
    const auto& jet = jets[i];
    std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();
    std::vector<int> constituentIndices;
    constituentIndices.reserve(jetConstituents.size());

    // Get the jet classification label later here! For now, we use "unknown"
    Acts::FastJet::JetLabel label = Acts::FastJet::JetLabel::Unknown;

    decltype(hadrons) hadronsInJet;

    if (m_cfg.doJetLabeling) {
      for (const auto& hadron : hadrons) {
        auto dR = deltaR(jet, hadron.first->momentum());
        if (dR < m_cfg.jetLabelingDeltaR) {
          hadronsInJet.emplace_back(hadron.first, hadron.second);
        }
      }
    }

    ACTS_DEBUG("Jet " << i << " has " << hadronsInJet.size() << " hadrons");
    if (logger().doPrint(Acts::Logging::VERBOSE)) {
      for (const auto& hadron : hadronsInJet) {
        auto dR = deltaR(jet, hadron.first->momentum());
        ACTS_VERBOSE("  - type: " << hadron.second << " dR: " << dR
                                  << " p=" << hadron.first->momentum().px()
                                  << ", " << hadron.first->momentum().py()
                                  << ", " << hadron.first->momentum().pz()
                                  << ", " << hadron.first->momentum().e());
      }
    }

    Acts::Vector4 jetFourMomentum(jets[i].px(), jets[i].py(), jets[i].pz(),
                                  jets[i].e());

    // Initialize the (track) jet with 4-momentum
    Acts::FastJet::TruthJetBuilder storedJet(jetFourMomentum);

    // Add the jet constituents to the (track)jet
    for (unsigned int j = 0; j < jetConstituents.size(); j++) {
      // Get the index of the constituent in the original input pseudo jets
      constituentIndices.push_back(jetConstituents[j].user_index());
    }

    Acts::FastJet::JetProperties jetProps(storedJet);
    jetProps.setConstituents(constituentIndices);

    outputJets.push_back(storedJet);
    ACTS_DEBUG("Stored jet properties: " << jetProps);
  }

  m_outputJets(ctx, std::move(outputJets));
  return ProcessCode::SUCCESS;
}

};  // namespace ActsExamples
