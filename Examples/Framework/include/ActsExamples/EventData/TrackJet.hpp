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
#include <vector>

namespace ActsExamples {

enum jetlabel { unknown = -99, ljet = 0, cjet = 4, bjet = 5 };

enum hadronlabel {
  Hadron = 1,
  BBbarMesonPart = 2,
  CCbarMesonPart = 3,
  BottomMesonPart = 4,
  BottomBaryonPart = 5,
  CharmedMesonPart = 6,
  CharmedBaryonPart = 7,
  StrangeMesonPart = 8,
  StrangeBaryonPart = 9,
  LightMesonPart = 10,
  LightBaryonPart = 11,
  Unknown = 12
};

class TrackJet {
 public:
  TrackJet(const Acts::Vector4& fm) {
    m_fourMomentum = fm;
    m_label = ActsExamples::unknown;
  }

  TrackJet(const Acts::Vector4& fm, const ActsExamples::jetlabel jl) {
    m_fourMomentum = fm;
    m_label = jl;
  }

  void setLabel(const ActsExamples::jetlabel jl) { m_label = jl; }

  ActsExamples::jetlabel getLabel() const { return m_label; }

  // TODO::Pass references instead of copies.

  void setConstituents(const std::vector<int> constituents) {
    m_constituents = constituents;
  }

  std::vector<int> getConstituents() const { return m_constituents; }

  Acts::Vector4 getFourMomentum() const { return m_fourMomentum; }

  void addTrack(const int trk_idx) { m_trk_idxs.push_back(trk_idx); }

  std::vector<int> getTracks() const { return m_trk_idxs; }

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

    std::cout << "With  " << m_trk_idxs.size() << " associated tracks"
              << std::endl;
    for (auto& trkidx : m_trk_idxs) {
      std::cout << trkidx << " ";
    }
    std::cout << std::endl;
  };

 private:
  Acts::Vector4 m_fourMomentum{0., 0., 0., 0.};
  ActsExamples::jetlabel m_label{ActsExamples::unknown};

  // The indices of the constituexonts wrt the global container
  std::vector<int> m_constituents{};

  // The indices of the tracks associated to this jet
  std::vector<int> m_trk_idxs{};
};

using TrackJetContainer = std::vector<TrackJet>;

}  // namespace ActsExamples
