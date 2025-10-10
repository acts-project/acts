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
#include "ActsExamples/EventData/SimParticle.hpp"

#include <optional>
#include <vector>

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

namespace Acts::FastJet {

/// Common class for jets
class Jet {
 public:
  explicit Jet(const Acts::Vector4& fourMom) { m_fourMomentum = fourMom; }

  /// @brief Get the jet 4-momentum
  /// @return the jet 4-momentum as an Acts::Vector4
  Acts::Vector4 getFourMomentum() const { return m_fourMomentum; }

  /// @brief Print the jet information
  friend std::ostream& operator<<(std::ostream& os, const Jet& jet) {
    os << "Jet 4-momentum: " << jet.getFourMomentum().transpose() << std::endl;
    return os;
  }

 private:
  Acts::Vector4 m_fourMomentum{0., 0., 0., 0.};
};

template <typename TrackContainer>
class TruthJet : public Jet {
 public:
  explicit TruthJet(const Acts::Vector4& fourMom) : Jet(fourMom) {}

  /// @brief Set the truth jet constituents
  void setConstituents(
      const std::vector<ActsExamples::SimParticle>& constituents) {
    m_constituents = constituents;
  }

  /// @brief Get the truth jet constituents
  const std::vector<ActsExamples::SimParticle>& getConstituents() const {
    return m_constituents;
  }

 private:
  /// @brief  The constituents of the truth jet are the truth particles
  std::vector<ActsExamples::SimParticle> m_constituents{};
  /// @brief The tracks associated to this truth jet
  std::vector<typename TrackContainer::TrackProxy> m_associatedTracks{};
};
template <typename TrackContainer>
class TrackJet : public Jet {
 public:
  explicit TrackJet(const Acts::Vector4& fourMom) : Jet(fourMom) {}

  /// @brief Set the track jet constituents
  void setConstituents(
      const std::vector<typename TrackContainer::TrackProxy>& constituents) {
    m_constituents = constituents;
  }

  /// @brief Get the track jet constituents
  const std::vector<typename TrackContainer::TrackProxy>& getConstituents()
      const {
    return m_constituents;
  }

 private:
  /// @brief The constituents of the track jet are the tracks
  std::vector<typename TrackContainer::TrackProxy> m_constituents{};
};

}  // namespace Acts::FastJet
