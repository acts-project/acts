// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/AnyTrackProxy.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <vector>

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

namespace ActsPlugins::FastJet {
/// @addtogroup fastjet_plugin
/// @{

/// Jet labeling enum

enum class JetLabel {
  Unknown = 0,
  BJet = 5,
  CJet = 4,
  LightJet = 1,
};

/// Output stream operator for JetLabel
/// @param os The output stream
/// @param label The jet label to output
/// @return The output stream
inline std::ostream& operator<<(std::ostream& os, const JetLabel& label) {
  switch (label) {
    case JetLabel::Unknown:
      os << "Unknown";
      break;
    case JetLabel::BJet:
      os << "BJet";
      break;
    case JetLabel::CJet:
      os << "CJet";
      break;
    case JetLabel::LightJet:
      os << "LightJet";
      break;
    default:
      os << "InvalidJetLabel";
      break;
  }
  return os;
}

/// Common class for jets
class Jet {
 public:
  /// Constructor
  /// @param fourMom The jet 4-momentum
  /// @param label The jet label
  Jet(const Acts::Vector4& fourMom, const JetLabel& label)
      : m_fourMomentum{fourMom}, m_jetLabel{label} {}

  /// @brief Get the jet 4-momentum
  /// @return the jet 4-momentum as an Acts::Vector4
  Acts::Vector4 fourMomentum() const { return m_fourMomentum; }

  /// @brief Get the jet label
  /// @return the jet label as JetLabel enum
  JetLabel jetLabel() const { return m_jetLabel; }

  /// @brief Set the jet label
  /// @param label the jet label as JetLabel enum
  void setJetLabel(const JetLabel& label) { m_jetLabel = label; }

 private:
  Acts::Vector4 m_fourMomentum{Acts::Vector4::Zero()};
  JetLabel m_jetLabel{JetLabel::Unknown};
  /// @brief Print the jet information
  friend std::ostream& operator<<(std::ostream& os, const Jet& jet) {
    os << "Jet 4-momentum: " << jet.fourMomentum().transpose() << std::endl;
    os << "Jet label: " << jet.jetLabel() << std::endl;
    return os;
  }
};

/// @brief Jet constructed from truth particles
class TruthJet : public Jet {
 public:
  /// Constructor
  /// @param fourMom The jet 4-momentum
  /// @param label The jet label
  TruthJet(const Acts::Vector4& fourMom, const JetLabel& label)
      : Jet(fourMom, label) {}

  /// @brief Set the truth particles as constituents of this truth jet from its barcode
  /// @param constituents Vector of barcodes for constituent particles
  void setConstituents(const std::vector<ActsFatras::Barcode>& constituents) {
    m_constituents = constituents;
  }

  /// @brief Set the indices of the truth particles that are constituents of this truth jet
  /// @param indices
  void setConstituentIndices(const std::vector<int>& indices) {
    m_constituentIndices = indices;
  }

  /// @brief Get the truth particles that are truth jet constituents
  /// @return Vector of barcodes for constituent particles
  const std::vector<ActsFatras::Barcode>& constituents() const {
    return m_constituents;
  }

  /// @brief Set the tracks associated to this truth jet
  /// @param associatedTracks Vector of tracks associated to this jet
  void setAssociatedTracks(
      const std::vector<Acts::AnyConstTrackProxy>& associatedTracks) {
    m_associatedTracks = associatedTracks;
  }

  /// @brief Get the tracks associated to this truth jet
  /// @return Vector of tracks associated to this jet
  const std::vector<Acts::AnyConstTrackProxy>& associatedTracks() const {
    return m_associatedTracks;
  }

 private:
  /// @brief  Truth particles as the constituents of the truth jet
  std::vector<ActsFatras::Barcode> m_constituents;
  /// @brief Indices of the constituents in the input collection
  std::vector<int> m_constituentIndices;
  /// @brief The tracks associated to this truth jet
  std::vector<Acts::AnyConstTrackProxy> m_associatedTracks;
};

/// @brief Jet constructed from reconstructed tracks
class TrackJet : public Jet {
 public:
  /// Constructor
  /// @param fourMom The jet 4-momentum
  /// @param label The jet label
  explicit TrackJet(const Acts::Vector4& fourMom, const JetLabel& label)
      : Jet(fourMom, label) {}

  /// @brief Set the tracks as constituents of this track jet
  /// @param constituents Vector of tracks that constitute this jet
  void setConstituents(
      const std::vector<Acts::AnyConstTrackProxy>& constituents) {
    m_constituents = constituents;
  }

  /// @brief Get the track jet constituents that are tracks
  /// @return Vector of tracks that constitute this jet
  const std::vector<Acts::AnyConstTrackProxy>& constituents() const {
    return m_constituents;
  }

 private:
  /// @brief Tracks as the constituents of the track jet
  std::vector<Acts::AnyConstTrackProxy> m_constituents;
};

/// @}
}  // namespace ActsPlugins::FastJet
