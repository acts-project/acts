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
#include "ActsExamples/Io/HepMC3/HepMC3Util.hpp"
#include "ActsExamples/Utilities/ParticleId.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <algorithm>
#include <fstream>
#include <mutex>
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

  m_inputTruthParticles.initialize(m_cfg.inputTruthParticles);

  if (m_cfg.doJetLabeling && !m_cfg.inputHepMC3Event.has_value()) {
    throw std::invalid_argument("Input HepMC3 event is not configured ");
  }

  m_inputHepMC3Event.initialize(m_cfg.inputHepMC3Event.value());
}

namespace {

constexpr int kBeamParticleStatus = 4;

JetLabel maxJetLabel(JetLabel a, JetLabel b) {
  return static_cast<JetLabel>(
      std::max(Acts::toUnderlying(a), Acts::toUnderlying(b)));
}

JetLabel jetLabelFromHadronType(ParticleId::HadronType hType) {
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

void findHadrons(
    const std::shared_ptr<const HepMC3::GenParticle>& particle,
    std::vector<std::shared_ptr<const HepMC3::GenParticle>>& hadrons,
    const Acts::Logger& logger, const std::string& prefix = "") {
  // JetLabel thisLabel = JetLabel::Unknown;  // light jet is our baseline, we
  //                                          // assume it IS a jet in any case

  // ACTS_VERBOSE(
  //     prefix << "-> particle="
  //            << Acts::findName(particle->pdg_id()).value_or("UNKNOWN"));

  bool isHadron = ActsExamples::ParticleId::isHadron(particle->pdg_id());
  if (isHadron) {
    auto hadronType = ActsExamples::ParticleId::hadronType(particle->pdg_id());

    auto label = jetLabelFromHadronType(hadronType);
    if (label != JetLabel::Unknown) {
      // ACTS_VERBOSE(prefix << "-> particle is a hadron with label=" << label);
      // ACTS_VERBOSE(prefix << "<== ADDING particle, don't recurse further");
      hadrons.push_back(particle);
      // do not recurse further, we have found a hadron
      return;
    }
  }

  if (const auto& vtx = particle->end_vertex(); vtx != nullptr) {
    for (const auto& child : vtx->particles_out()) {
      findHadrons(child, hadrons, logger, prefix + "  ");
    }
  }
  // ACTS_VERBOSE(prefix << "<-");
}

}  // namespace

ProcessCode ActsExamples::TruthJetAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  Acts::FastJet::TrackJetContainer outputJets;

  const SimParticleContainer& truthParticlesRaw = m_inputTruthParticles(ctx);
  std::vector<const SimParticle*> truthParticles;
  std::ranges::transform(truthParticlesRaw, std::back_inserter(truthParticles),
                         [](const auto& particle) { return &particle; });

  // std::vector<const HepMC3::GenParticle*> truthParticlesHepMC3;
  // findInputParticles(genEvent, truthParticlesHepMC3);

  ACTS_DEBUG("Number of truth particles: " << truthParticles.size());
  // ACTS_DEBUG(
  //     "Number of truth particles (HepMC3): " << truthParticlesHepMC3.size());

  // if (truthParticlesHepMC3.size() != truthParticles.size()) {
  //   ACTS_ERROR("Number of truth particles (HepMC3) and (SimParticle)
  //   differ"); return ProcessCode::ABORT;
  // }

  const fastjet::JetDefinition defaultJetDefinition =
      fastjet::JetDefinition(fastjet::antikt_algorithm, m_cfg.jetClusteringR);

  // Get the 4-momentum information from the simulated truth particles
  // and create fastjet::PseudoJet objects
  std::vector<fastjet::PseudoJet> inputPseudoJets;

  static std::mutex mtxPseudoJets;
  {
    std::lock_guard lock(mtxPseudoJets);
    std::ofstream outfile;
    outfile.open("pseudojets.csv",
                 std::ios_base::app);  // append instead of overwrite

    for (unsigned int i = 0; i < truthParticles.size(); i++) {
      const auto* particle = truthParticles.at(i);
      //   fastjet::PseudoJet pseudoJet(
      //       particle->momentum().px(), particle->momentum().py(),
      // particle->momentum().pz(), particle->momentum().e());

      fastjet::PseudoJet pseudoJet(
          particle->momentum().x(), particle->momentum().y(),
          particle->momentum().z(), particle->energy());

      outfile << ctx.eventNumber << "," << pseudoJet.pt() << ","
              << pseudoJet.eta() << "," << pseudoJet.phi() << ","
              << particle->pdg();
      outfile << std::endl;

    pseudoJet.set_user_index(i);
    inputPseudoJets.push_back(pseudoJet);
  }
  ACTS_DEBUG("Number of input pseudo jets: " << inputPseudoJets.size());

  std::vector<fastjet::PseudoJet> jets;
  fastjet::ClusterSequence clusterSeq;
  {
    Acts::ScopedTimer timer("Jet clustering", logger(), Acts::Logging::DEBUG);
    // Run the jet clustering, only once
    clusterSeq =
        fastjet::ClusterSequence(inputPseudoJets, defaultJetDefinition);

<<<<<<< HEAD
  // Get the jets above a certain pt threshold
  std::vector<fastjet::PseudoJet> jets =
      sorted_by_pt(clusterSeq.inclusive_jets(m_cfg.jetPtMin));
  ACTS_DEBUG("Number of clustered truth jets: " << jets.size());
=======
    // Get the jets above a certain pt threshold
    jets = sorted_by_pt(clusterSeq.inclusive_jets(m_cfg.jetPtMin));
    ACTS_DEBUG("Number of clustered jets: " << jets.size());
  }

  std::vector<std::pair<const HepMC3::GenParticle*, ParticleId::HadronType>>
      hadrons;
  if (m_cfg.doJetLabeling) {
    ACTS_DEBUG("Jet labeling is enabled");
    const auto& genEvent = *m_inputHepMC3Event(ctx);

    // for (const auto& vertex : genEvent.vertices()) {
    //   if (vertex->particles_in().empty() ||
    //       std::ranges::all_of(vertex->particles_in(), [](const auto&
    //       particle) {
    //         return particle->status() == kBeamParticleStatus;
    //       })) {
    //     ACTS_VERBOSE("PRIMARY VERTEX: ");
    //     HepMC3::Print::line(vertex);
    //     ACTS_VERBOSE("-> incoming particles: ");
    //     for (const auto& particle : vertex->particles_in()) {
    //       std::stringstream ss;
    //       HepMC3::Print::line(ss, particle);
    //       ACTS_VERBOSE("  - " << ss.str());
    //     }

    //     std::size_t nHadrons = hadrons.size();
    //     for (const auto& particle : vertex->particles_out()) {
    //       findHadrons(particle, hadrons, logger());
    //       ACTS_VERBOSE("-> out n_hadrons=" << (hadrons.size() - nHadrons));
    //     }
    //   }
    // }

    // ACTS_VERBOSE("-> total n_hadrons=" << hadrons.size());
    // for (const auto& hadron : hadrons) {
    //   ACTS_VERBOSE(
    //       "  - " << hadron->pdg_id() << " "
    //              << Acts::findName(hadron->pdg_id()).value_or("UNKNOWN")
    //              << " label="
    //              <<
    //              jetLabelFromHadronType(ActsExamples::ParticleId::hadronType(
    //                     hadron->pdg_id())));
    //   std::stringstream ss;
    //   HepMC3::Print::line(ss, hadron);
    //   ACTS_VERBOSE("  - " << ss.str());
    // }

    // for (const auto& particle : genEvent.particles()) {
    //   if (!ParticleId::isHadron(particle->pdg_id())) {
    //     continue;
    //   }
    //   // ACTS_VERBOSE("Particle " << particle->pdg_id() << " with status "
    //   // HepMC3::Print::line(particle);

    //   auto hadronType =
    //       ActsExamples::ParticleId::hadronType(particle->pdg_id());
    //   hadrons.emplace_back(hadronType, particle.get());
    //   // auto hadronLabel = ParticleId::hadronLabel(particle->pdg_id());
    //   // std::cout << "Hadron label: " << hadronLabel << std::endl;
    // }

    auto hadronView =
        genEvent.particles() | std::views::filter([](const auto& particle) {
          return ParticleId::isHadron(particle->pdg_id());
        }) |
        std::views::transform([](const auto& particle) {
          auto type = ActsExamples::ParticleId::hadronType(particle->pdg_id());
          auto label = jetLabelFromHadronType(type);
          return std::pair{label, particle};
        }) |
        std::views::filter([](const auto& hadron) {
          return hadron.first > JetLabel::LightJet;
        });

    std::ranges::copy(hadronView, std::back_inserter(hadrons));

    // deduplicate hadrons
    std::ranges::sort(hadrons, [](const auto& a, const auto& b) {
      return a.second->id() < b.second->id();
    });
    auto unique = std::ranges::unique(hadrons);
    hadrons.erase(unique.begin(), unique.end());
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
      ACTS_VERBOSE("Found jet "
                   << i << " with 4-momentum: " << jetFourMomentum(0) << ", "
                   << jetFourMomentum(1) << ", " << jetFourMomentum(2) << ", "
                   << jetFourMomentum(3) << " and " << constituentIndices.size()
                   << " constituents.");

      JetLabel label = JetLabel::Unknown;
      if (m_cfg.doJetLabeling) {
        timer.sample();
        label = classifyJet(jet);
      }

      outfile << ctx.eventNumber << "," << jet.pt() << "," << jet.eta() << ","
              << jet.phi() << "," << static_cast<int>(label);
      outfile << std::endl;

      // Initialize the (track) jet with 4-momentum and jet label
      Acts::FastJet::TruthJetBuilder storedJet(jetFourMomentum, label);

      // Add the jet constituents to the (track)jet
      for (unsigned int j = 0; j < jetConstituents.size(); j++) {
        // Get the index of the constituent in the original input pseudo jets
        constituentIndices.push_back(jetConstituents[j].user_index());
      }

      Acts::FastJet::JetProperties jetProps(storedJet);
    jetProps.setConstituents(constituentIndices);

    outputJets.push_back(storedJet);

      jetLabelCounts[label] += 1;

      ACTS_VERBOSE("-> jet label: " << label);
      ACTS_VERBOSE("-> jet constituents: ");

      // if (logger().doPrint(Acts::Logging::VERBOSE)) {
      //   for (const auto& constituent : constituentIndices) {
      //     const auto& particle = inputParticles.at(constituent);
      //     ACTS_VERBOSE("- " << particle);
      //   }
      // }
    }

    outfile.flush();
    outfile.close();
  }

  ACTS_DEBUG("-> jet label counts: ");
  for (const auto& [label, count] : jetLabelCounts) {
    ACTS_DEBUG("  - " << label << ": " << count);
  }

  m_outputJets(ctx, std::move(outputJets));
  return ProcessCode::SUCCESS;
}
};  // namespace ActsExamples
