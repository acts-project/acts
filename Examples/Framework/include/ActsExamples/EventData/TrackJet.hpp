// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <iostream>
#include <ostream>
#include <vector>

namespace ActsExamples {

enum class JetLabel { Unknown = -99, LightJet = 0, CJet = 4, BJet = 5 };

inline std::ostream& operator<<(std::ostream& os, const JetLabel& label) {
  switch (label) {
    case JetLabel::Unknown:
      os << "Unknown";
      break;
    case JetLabel::LightJet:
      os << "LightJet";
      break;
    case JetLabel::CJet:
      os << "CJet";
      break;
    case JetLabel::BJet:
      os << "BJet";
      break;
  }
  return os;
}

class TrackJet {
 public:
  explicit TrackJet(const Acts::Vector4& fourMomentum)
      : m_fourMomentum(fourMomentum) {}

  TrackJet(const Acts::Vector4& fourMomentum,
           const ActsExamples::JetLabel label)
      : m_fourMomentum(fourMomentum), m_label(label) {}

  void setLabel(ActsExamples::JetLabel jl) { m_label = jl; }

  ActsExamples::JetLabel getLabel() const { return m_label; }

  // TODO::Pass references instead of copies.

  void setConstituents(const std::vector<int>& constituents) {
    m_constituents = constituents;
  }

  std::vector<int> getConstituents() const { return m_constituents; }

  Acts::Vector4 getFourMomentum() const { return m_fourMomentum; }

  Acts::Vector3 getDirection() const {
    return m_fourMomentum.head<3>().normalized();
  }

  void addTrack(const int trk_idx) { m_trackIndices.push_back(trk_idx); }

  std::vector<int> getTracks() const { return m_trackIndices; }

  // @TODO: Turn into an output stream operator
  void print() const {
    std::cout << "Printing Jet information" << std::endl;
    std::cout << "4mom=(" << m_fourMomentum(0) << "," << m_fourMomentum(1)
              << "," << m_fourMomentum(2) << "," << m_fourMomentum(3) << ")"
              << std::endl;

    std::cout << "Formed by " << m_constituents.size() << std::endl;

    for (auto& constituent : m_constituents) {
      std::cout << constituent << " ";
    }
    std::cout << std::endl;

    std::cout << "With  " << m_trackIndices.size() << " associated tracks"
              << std::endl;
    for (auto& trkidx : m_trackIndices) {
      std::cout << trkidx << " ";
    }
    std::cout << std::endl;
  };

 private:
  Acts::Vector4 m_fourMomentum{0., 0., 0., 0.};
  ActsExamples::JetLabel m_label{ActsExamples::JetLabel::Unknown};

  // The indices of the constituexonts wrt the global container
  std::vector<int> m_constituents{};

  // The indices of the tracks associated to this jet
  std::vector<int> m_trackIndices{};
};

using TrackJetContainer = std::vector<TrackJet>;

}  // namespace ActsExamples
