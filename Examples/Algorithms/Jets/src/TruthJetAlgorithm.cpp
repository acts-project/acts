// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Jets/TruthJetAlgorithm.hpp"

#include "Acts/Definitions/ParticleData.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <algorithm>
#include <ranges>
#include <stdexcept>

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
  if (m_cfg.doTrackJetMatching) {
    m_inputTracks.initialize(m_cfg.inputTracks);
  }
  m_outputJets.initialize(m_cfg.outputJets);
}

namespace {
ActsPlugins::FastJet::JetLabel jetLabelFromHadronType(
    Acts::HadronType hadronType) {
  using enum Acts::HadronType;
  switch (hadronType) {
    case BBbarMeson:
    case BottomMeson:
    case BottomBaryon:
      return ActsPlugins::FastJet::JetLabel::BJet;
    case CCbarMeson:
    case CharmedMeson:
    case CharmedBaryon:
      return ActsPlugins::FastJet::JetLabel::CJet;
    case StrangeMeson:
    case StrangeBaryon:
    case LightMeson:
    case LightBaryon:
      return ActsPlugins::FastJet::JetLabel::LightJet;
    default:
      return ActsPlugins::FastJet::JetLabel::Unknown;
  }
}

}  // namespace

ProcessCode TruthJetAlgorithm::execute(const AlgorithmContext& ctx) const {
  // Initialize the output container
  std::vector<ActsPlugins::FastJet::TruthJet> outputJetContainer{};

  Acts::ScopedTimer globalTimer("TruthJetAlgorithm", logger(),
                                Acts::Logging::DEBUG);

  const SimParticleContainer& truthParticlesRaw = m_inputTruthParticles(ctx);
  std::vector<const SimParticle*> truthParticles;
  truthParticles.reserve(truthParticlesRaw.size());
  std::ranges::transform(truthParticlesRaw, std::back_inserter(truthParticles),
                         [](const auto& particle) { return &particle; });
  ACTS_DEBUG("Number of truth particles: " << truthParticles.size());

  const fastjet::JetDefinition defaultJetDefinition = fastjet::JetDefinition(
      fastjet::antikt_algorithm, m_cfg.jetClusteringRadius);

  // Get the 4-momentum information from the simulated truth particles
  // and create fastjet::PseudoJet objects
  std::vector<fastjet::PseudoJet> inputPseudoJets;
  {
    Acts::ScopedTimer timer("Input particle building", logger(),
                            Acts::Logging::DEBUG);

    for (unsigned int i = 0; i < truthParticles.size(); ++i) {
      const auto* particle = truthParticles.at(i);

      detail::PrimaryVertexIdGetter primaryVertexIdGetter;

      ACTS_VERBOSE("Primary vertex ID: "
                   << primaryVertexIdGetter(*particle).vertexPrimary()
                   << ", PDG: " << static_cast<int>(particle->pdg()) << ", pT: "
                   << particle->transverseMomentum() / Acts::UnitConstants::GeV
                   << " GeV");

      if (m_cfg.clusterHSParticlesOnly &&
          primaryVertexIdGetter(*particle).vertexPrimary() != 1) {
        continue;
      }

      fastjet::PseudoJet pseudoJet(
          particle->momentum().x(), particle->momentum().y(),
          particle->momentum().z(), particle->energy());

      pseudoJet.set_user_index(i);
      inputPseudoJets.push_back(pseudoJet);
    }
  }

  ACTS_DEBUG("Number of input pseudo jets: " << inputPseudoJets.size());

  std::vector<fastjet::PseudoJet> jets;
  fastjet::ClusterSequence clusterSeq;
  {
    Acts::ScopedTimer timer("Jet clustering", logger(), Acts::Logging::DEBUG);
    // Run the jet clustering
    clusterSeq =
        fastjet::ClusterSequence(inputPseudoJets, defaultJetDefinition);
    // Get the jets above a certain pt threshold
    jets = sorted_by_pt(clusterSeq.inclusive_jets(m_cfg.jetPtMin));
    // Apply eta range cut if specified
    double minEta = m_cfg.jetEtaRange.first;
    double maxEta = m_cfg.jetEtaRange.second;

    std::erase_if(jets, [minEta, maxEta](const auto& jet) {
      return jet.eta() < minEta || jet.eta() > maxEta;
    });
    ACTS_DEBUG("Number of clustered jets: " << jets.size());
  }

  // Find hadrons for jet labeling with sim particles
  std::vector<std::pair<ActsPlugins::FastJet::JetLabel, const SimParticle*>>
      hadrons;

  if (m_cfg.doJetLabeling) {
    Acts::ScopedTimer timer("hadron finding", logger(), Acts::Logging::DEBUG);
    ACTS_DEBUG("Jet labeling is enabled. Finding hadrons for jet labeling.");

    // A lazy view over the simulated particles for hadron finding
    auto hadronView =
        truthParticles | std::views::filter([this](const auto* particle) {
          if (m_cfg.jetLabelingHSHadronsOnly) {
            detail::PrimaryVertexIdGetter primaryVertexIdGetter;
            if (primaryVertexIdGetter(*particle).vertexPrimary() != 1) {
              return false;
            }
          }

          Acts::PdgParticle pdgId{particle->pdg()};
          if (!Acts::ParticleIdHelper::isHadron(pdgId)) {
            return false;
          }

          // Apply a pt cut on B or C hadrons
          auto label =
              jetLabelFromHadronType(Acts::ParticleIdHelper::hadronType(pdgId));
          using enum ActsPlugins::FastJet::JetLabel;

          if (label == BJet || label == CJet) {
            if (particle->transverseMomentum() < m_cfg.jetLabelingHadronPtMin) {
              return false;
            }
          }
          return true;
        }) |
        std::views::transform([](const auto* particle) {
          Acts::PdgParticle pdgId{particle->pdg()};
          auto type = Acts::ParticleIdHelper::hadronType(pdgId);
          auto label = jetLabelFromHadronType(type);
          return std::make_pair(label, particle);
        }) |
        std::views::filter([](const auto& hadron) {
          return hadron.first > ActsPlugins::FastJet::JetLabel::Unknown;
        });

    std::ranges::copy(hadronView, std::back_inserter(hadrons));

    // Deduplicate hadrons based on their pdg id
    std::ranges::sort(hadrons, [](const auto& a, const auto& b) {
      return a.second->pdg() < b.second->pdg();
    });

    auto unique = std::ranges::unique(hadrons);
    hadrons.erase(unique.begin(), unique.end());
  }

  // Jet classification

  auto classifyJet = [&](const fastjet::PseudoJet& jet) {
    auto hadronsInJetView =
        hadrons | std::views::filter([&jet, this](const auto& hadron) {
          const auto& momentum = hadron.second->momentum();
          Acts::Vector3 hadronJetMom{momentum[0], momentum[1], momentum[2]};
          Acts::Vector3 jetMom{jet.px(), jet.py(), jet.pz()};
          return Acts::VectorHelpers::deltaR(jetMom, hadronJetMom) <
                 m_cfg.jetLabelingDeltaR;
        }) |
        std::views::transform([](const auto& hadron) {
          return std::pair{
              hadron.second,
              jetLabelFromHadronType(Acts::ParticleIdHelper::hadronType(
                  Acts::PdgParticle{hadron.second->pdg()}))};
        });

    std::vector<std::pair<const SimParticle*, ActsPlugins::FastJet::JetLabel>>
        hadronsInJet;
    std::ranges::copy(hadronsInJetView, std::back_inserter(hadronsInJet));

    ACTS_VERBOSE("-> hadrons in jet: " << hadronsInJet.size());
    for (const auto& hadron : hadronsInJet) {
      ACTS_VERBOSE(
          "  - " << hadron.first->pdg() << " "
                 << Acts::findName(Acts::PdgParticle{hadron.first->pdg()})
                        .value_or("UNKNOWN")
                 << " label=" << hadron.second);
    }

    auto maxHadronIt = std::ranges::max_element(
        hadronsInJet, [](const auto& a, const auto& b) { return a < b; },
        [](const auto& a) {
          const auto& [hadron, label] = a;
          return label;
        });

    if (maxHadronIt == hadronsInJet.end()) {
      // Now hadronic "jet"
      return ActsPlugins::FastJet::JetLabel::Unknown;
    }

    const auto& [maxHadron, maxHadronLabel] = *maxHadronIt;

    ACTS_VERBOSE("-> max hadron type="
                 << Acts::findName(Acts::PdgParticle{maxHadron->pdg()})
                        .value_or("UNKNOWN")
                 << " label=" << maxHadronLabel);

    return maxHadronLabel;
  };  // jet classification

  boost::container::flat_map<ActsPlugins::FastJet::JetLabel, std::size_t>
      jetLabelCounts;

  {
    Acts::AveragingScopedTimer timer("Jet classification", logger(),
                                     Acts::Logging::DEBUG);

    for (unsigned int i = 0; i < jets.size(); ++i) {
      const auto& jet = jets.at(i);

      // If jet labeling is enabled, classify the jet based on its hadronic
      // content
      ActsPlugins::FastJet::JetLabel jetLabel =
          ActsPlugins::FastJet::JetLabel::Unknown;
      if (m_cfg.doJetLabeling) {
        ACTS_DEBUG("Classifying jet " << i);
        auto sample = timer.sample();
        jetLabel = classifyJet(jet);
      }

      // Initialize truth jet for storing in the output container
      Acts::Vector4 jetFourMom{jet.px(), jet.py(), jet.pz(), jet.e()};
      ActsPlugins::FastJet::TruthJet truthJet(jetFourMom, jetLabel);

      std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();
      std::vector<int> constituentIndices;
      constituentIndices.reserve(jetConstituents.size());

      for (unsigned int j = 0; j < jetConstituents.size(); ++j) {
        constituentIndices.push_back(jetConstituents[j].user_index());
      }

      truthJet.setConstituentIndices(constituentIndices);
      outputJetContainer.push_back(truthJet);

      jetLabelCounts[jetLabel] += 1;

      ACTS_VERBOSE("-> jet label: " << jetLabel);
      ACTS_VERBOSE("-> jet constituents: ");

      if (logger().doPrint(Acts::Logging::VERBOSE)) {
        for (const auto& constituent : constituentIndices) {
          const auto& particle = truthParticles.at(constituent);
          ACTS_VERBOSE("- " << particle);
        }
      }
    }
  }

  ACTS_DEBUG("-> jet label counts: ");
  for (const auto& [label, count] : jetLabelCounts) {
    ACTS_DEBUG("  - " << label << ": " << count);
  }

  if (m_cfg.doTrackJetMatching) {
    const ConstTrackContainer& tracks = m_inputTracks(ctx);
    trackJetMatching(tracks, outputJetContainer);
  }

  m_outputJets(ctx, std::move(outputJetContainer));

  return ProcessCode::SUCCESS;
}

ProcessCode TruthJetAlgorithm::finalize() {
  ACTS_INFO("Finalizing truth jet clustering");
  return ProcessCode::SUCCESS;
}

void TruthJetAlgorithm::trackJetMatching(const ConstTrackContainer& tracks,
                                         TruthJetContainer& jets) const {
  std::unordered_map<std::size_t, std::vector<std::int32_t>>
      jetToTrackIndicesMap;

  for (auto track : tracks) {
    // Find the closest jet to this track and associate if within minDeltaR
    double minDeltaR = 0.4;
    std::int32_t closestJetIndex = -1;

    for (std::size_t ijet = 0; ijet < jets.size(); ++ijet) {
      Acts::Vector4 jet_4mom = jets[ijet].fourMomentum();
      Acts::Vector3 jet_3mom{jet_4mom[0], jet_4mom[1], jet_4mom[2]};

      Acts::Vector4 track_4mom = track.fourMomentum();
      Acts::Vector3 track_3mom{track_4mom[0], track_4mom[1], track_4mom[2]};

      // Calculate deltaR between track and jet
      auto drTrackJet = Acts::VectorHelpers::deltaR(jet_3mom, track_3mom);

      if (drTrackJet < minDeltaR) {
        minDeltaR = drTrackJet;
        closestJetIndex = ijet;
      }
    }
    if (closestJetIndex != -1) {
      jetToTrackIndicesMap[closestJetIndex].push_back(track.index());
    }
  }

  for (std::size_t ijet = 0; ijet < jets.size(); ++ijet) {
    std::size_t nTracksAssociated = 0;
    std::vector<Acts::AnyConstTrackProxy> associatedTracks = {};
    auto search = jetToTrackIndicesMap.find(ijet);
    if (search != jetToTrackIndicesMap.end()) {
      nTracksAssociated = search->second.size();
      ACTS_VERBOSE("Jet " << ijet << " has " << nTracksAssociated
                          << " associated tracks with indices:");
      for (auto trackIdx : search->second) {
        ACTS_VERBOSE("  Track index: " << trackIdx);
        auto constTrack = tracks.getTrack(trackIdx);
        Acts::AnyConstTrackProxy trackProxy(constTrack);
        associatedTracks.push_back(trackProxy);
      }
    }
    jets[ijet].setAssociatedTracks(associatedTracks);
  }
}

}  // namespace ActsExamples
