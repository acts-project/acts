// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Jets/TrackToTruthJetAlgorithm.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <ostream>
#include <stdexcept>

namespace ActsExamples {

TrackToTruthJetAlgorithm::TrackToTruthJetAlgorithm(const Config& cfg,
                                                   Acts::Logging::Level lvl)
    : IAlgorithm("TrackToTruthJetAlgorithm", lvl), m_cfg(cfg) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Input tracks are not configured");
  }
  if (m_cfg.inputJets.empty()) {
    throw std::invalid_argument("Input tracks are not configured");
  }
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputJets.initialize(m_cfg.inputJets);
  m_outputTrackJets.initialize(m_cfg.outputTrackJets);
}

ProcessCode ActsExamples::TrackToTruthJetAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  ACTS_INFO("Executing track to truth jet matching algorithm");

  const auto& tracks = m_inputTracks(ctx);
  const auto& truthJets = m_inputJets(ctx);
  TrackJetContainer jets; // track jets to be filled
  //copy truth jets to jets
  for (const auto& jet : truthJets) {
    jets.emplace_back(jet.getFourMomentum(), jet.getLabel());
    jets.back().setConstituents(jet.getConstituents());
  }

  ACTS_INFO("TrackToTruthJetAlg - Number of tracks: " << tracks.size());
  ACTS_INFO("TrackToTruthJetAlg - Number of truth jets: " << truthJets.size());
  ACTS_INFO("TrackToTruthJetAlg - Number of jets: " << jets.size());

  for (const auto& track : tracks) {
    double minDeltaR = m_cfg.maxDeltaR;
    int closestJetIndex = -1;

    double trackEnergy =
        sqrt(track.absoluteMomentum() *
             track.absoluteMomentum());  // need to add mass here!

    // ACTS_DEBUG("Track index: "
    //            << track.index() << ", momentum: " << track.momentum().x()
    //            << ", " << track.momentum().y() << ", " << track.momentum().z()
    //            << ", energy: " << trackEnergy);

    // Create a fastjet::PseudoJet object for the track
    fastjet::PseudoJet trackJet(track.momentum().x(), track.momentum().y(),
                                track.momentum().z(), trackEnergy);
    trackJet.set_user_index(track.index());

    // Loop over the jets to find the closest one
    for (std::size_t i = 0; i < jets.size(); ++i) {
      const auto& jet = jets[i];

      // Create a fastjet::PseudoJet object for the jet
      fastjet::PseudoJet jetPseudo(jet.getFourMomentum()(0),
                                   jet.getFourMomentum()(1),
                                   jet.getFourMomentum()(2),
                                   jet.getFourMomentum()(3));

      // TODO: instead pseudo jet use directly eta phi from trackjet 
      if (trackJet.delta_R(jetPseudo) < minDeltaR) {
        minDeltaR = trackJet.delta_R(jetPseudo);
        closestJetIndex = i;
      }

      ACTS_DEBUG("Track " << track.index() << " delta R to jet " << i << ": "
                          << trackJet.delta_R(jetPseudo));
      if(closestJetIndex != -1) {
        ACTS_DEBUG("Closest jet so far: " << closestJetIndex
                                          << " with delta R: " << minDeltaR);
      }

    }  // loop over jets

    if (closestJetIndex != -1) {
      ACTS_DEBUG("Adding track " << track.index() << " to jet " << closestJetIndex );
      jets[closestJetIndex].addTrack(track.index());
    }

  }  // loop over tracks
                
  // Write the matched jets to the output
  m_outputTrackJets(ctx, std::move(jets));


  return ProcessCode::SUCCESS;
}

ProcessCode ActsExamples::TrackToTruthJetAlgorithm::finalize() const {
  ACTS_INFO("Finalizing track to truth jet matching algorithm");
  return ProcessCode::SUCCESS;
}

};  // namespace ActsExamples
