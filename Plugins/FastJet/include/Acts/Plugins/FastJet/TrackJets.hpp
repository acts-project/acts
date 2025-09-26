// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
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

class TruthJetBuilder {
 public:
  explicit TruthJetBuilder(const Acts::Vector4& fm) { m_fourMomentum = fm; }

 private:
  Acts::Vector4 m_fourMomentum{0., 0., 0., 0.};
};

class JetProperties {
 public:
  /// Constructor; saves a reference to the jet
  /// @param jet the jet
  explicit JetProperties(const fastjet::PseudoJet& jet) : m_jet{jet} {}

  /// @brief Set the jet constituents
  /// @param constituents the indices of the constituent tracks
  void setConstituents(const std::vector<int>& constituents) {
    m_constituents = constituents;
  }
  /// @brief Get the jet constituents
  /// @return the indices of the constituent tracks
  const std::vector<int>& getConstituents() const { return m_constituents; }

  /// @brief Get the jet 4-momentum
  /// @return the jet 4-momentum as an Acts::Vector4
  Acts::Vector4 getFourMomentum() const { return m_fourMomentum; }

  /// @brief Add a track to the jet
  /// @param trk_idx the index of the track to add
  void addTrack(const int trackIndex) { m_trackIndices.push_back(trackIndex); }

  /// @brief Get the tracks associated to this jet
  /// @return the indices of the associated tracks
  const std::vector<int>& getTracks() const { return m_trackIndices; }

  /// @brief Print the jet information
  friend std::ostream& operator<<(std::ostream& os,
                                  const JetProperties& jetProps) {
    os << "Jet 4-momentum: " << jetProps.getFourMomentum().transpose()
       << std::endl;
    os << "Constituents: ";
    for (const auto& constituent : jetProps.getConstituents()) {
      os << constituent << " ";
    }
    os << std::endl;
    if (!jetProps.getTracks().empty()) {
      os << "Associated tracks: ";
      for (const auto& trkidx : jetProps.getTracks()) {
        os << trkidx << " ";
      }
      os << std::endl;
    } else {
      os << "No associated tracks." << std::endl;
    }
    return os;
  }

 private:
  const fastjet::PseudoJet& m_jet;
  Acts::Vector4 m_fourMomentum{m_jet.px(), m_jet.py(), m_jet.pz(), m_jet.e()};
  // The indices of the constituents wrt the global container
  std::vector<int> m_constituents{};
  // The indices of the tracks associated to this jet
  std::vector<int> m_trackIndices{};
};

using TrackJetContainer = std::vector<TruthJetBuilder>;

}  // namespace Acts::FastJet

#include "TrackJets.ipp"
