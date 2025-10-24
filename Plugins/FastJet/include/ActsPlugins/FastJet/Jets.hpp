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
#include "ActsFatras/EventData/Barcode.hpp"

#include <optional>
#include <vector>

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

namespace ActsPlugins::FastJet {

/// Common class for jets
class Jet {
 public:
  explicit Jet(const Acts::Vector4& fourMom) : m_fourMomentum{fourMom} {}

  /// @brief Get the jet 4-momentum
  /// @return the jet 4-momentum as an Acts::Vector4
  Acts::Vector4 fourMomentum() const { return m_fourMomentum; }

 private:
  Acts::Vector4 m_fourMomentum{Acts::Vector4::Zero()};
  /// @brief Print the jet information
  friend std::ostream& operator<<(std::ostream& os, const Jet& jet) {
    os << "Jet 4-momentum: " << jet.fourMomentum().transpose() << std::endl;
    return os;
  }
};

template <typename TrackContainer>
class TruthJet : public Jet {
 public:
  explicit TruthJet(const Acts::Vector4& fourMom) : Jet(fourMom) {}

  /// @brief Set the truth particles as constituents of this truth jet from its barcode
  void setConstituents(const std::vector<ActsFatras::Barcode>& constituents) {
    m_constituents = constituents;
  }

  /// @brief Get the truth particles that are truth jet constituents
  const std::vector<ActsFatras::Barcode>& constituents() const {
    return m_constituents;
  }

  /// @brief Set the tracks associated to this truth jet
  void setAssociatedTracks(
      const std::vector<typename TrackContainer::TrackProxy>&
          associatedTracks) {
    m_associatedTracks = associatedTracks;
  }

  /// @brief Get the tracks associated to this truth jet
  const std::vector<typename TrackContainer::TrackProxy>& associatedTracks()
      const {
    return m_associatedTracks;
  }

 private:
  /// @brief  Truth particles as the constituents of the truth jet
  std::vector<ActsFatras::Barcode> m_constituents;
  /// @brief The tracks associated to this truth jet
  std::vector<typename TrackContainer::TrackProxy> m_associatedTracks;
};

template <typename TrackContainer>
class TrackJet : public Jet {
 public:
  explicit TrackJet(const Acts::Vector4& fourMom) : Jet(fourMom) {}

  /// @brief Set the tracks as constituents of this track jet
  void setConstituents(
      const std::vector<typename TrackContainer::TrackProxy>& constituents) {
    m_constituents = constituents;
  }

  /// @brief Get the track jet constituents that are tracks
  const std::vector<typename TrackContainer::TrackProxy>& constituents() const {
    return m_constituents;
  }

 private:
  /// @brief Tracks as the constituents of the track jet
  std::vector<typename TrackContainer::TrackProxy> m_constituents;
};

}  // namespace ActsPlugins::FastJet
