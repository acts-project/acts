// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Units.hpp"

#include <optional>
#include <vector>

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

namespace Acts::FastJet {

/// Default jet definition: Anti-kt with a radius of 0.4
const fastjet::JetDefinition DefaultJetDefinition =
    fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4);

template <typename TrackContainer>
class InputTracks {
 public:
  /// Constructor; saves a reference to the container
  /// @param tracks the TrackContainer
  explicit InputTracks(TrackContainer& tracks) : m_tracks{tracks} {}

  /// Get vector for 4-momenta from the track collection
  /// @return vector of fastjet::PseudoJet, one per track
  std::vector<fastjet::PseudoJet> fourMomenta() const;

  /// Get the tracks making up a track-jet
  ///
  /// @param jet the jet from which to get the constituent tracks
  /// @param coreR optional radius inside which to get the tracks
  ///
  /// @return a vector of TrackProxy
  std::vector<typename TrackContainer::TrackProxy> tracksInJet(
      const fastjet::PseudoJet& jet, std::optional<float> coreR = {});

 private:
  TrackContainer& m_tracks;
};

class TrackJetSequence {
 public:
  /// Factory function to create a sequence of track jets
  ///
  /// @param tracks the input tracks
  /// @jetDef the jet definition to use, defaults to "DefaultJetDefinition"
  static TrackJetSequence create(
      std::vector<fastjet::PseudoJet>& tracks,
      const fastjet::JetDefinition& jetDef = DefaultJetDefinition);

  static TrackJetSequence create(
      std::vector<fastjet::PseudoJet>&& tracks,
      const fastjet::JetDefinition& jetDef = DefaultJetDefinition) {
    return create(tracks, jetDef);
  }

  /// Get all the track jets passing the pT & eta cuts
  ///
  /// @param ptMin the minimum jet pT in GeV
  /// @param etaMax the maximum jet absolute eta
  ///
  /// @return a vector of fastjet::PseudoJet objects
  std::vector<fastjet::PseudoJet> jets(float ptMin = 20 *
                                                     Acts::UnitConstants::GeV,
                                       float etaMax = 2.5);

 private:
  /// Main constructor. Users should call "TrackJetSequence::create" instead
  ///
  /// @param clusterSeq the fastjet::ClusterSequence object
  explicit TrackJetSequence(const fastjet::ClusterSequence& clusterSeq)
      : m_clusterSeq{clusterSeq} {}

  fastjet::ClusterSequence m_clusterSeq{};
};

}  // namespace Acts::FastJet

#include "TrackJets.ipp"
