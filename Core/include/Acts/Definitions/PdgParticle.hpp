// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <utility>

namespace Acts {

/// Symbolic values for commonly used PDG particle numbers.
enum PdgParticle : std::int32_t {
  eInvalid = 0,
  eElectron = 11,
  eAntiElectron = -eElectron,
  ePositron = -eElectron,
  eMuon = 13,
  eAntiMuon = -eMuon,
  eTau = 15,
  eAntiTau = -eTau,
  eGamma = 22,
  ePionZero = 111,
  ePionPlus = 211,
  ePionMinus = -ePionPlus,
  eKaonPlus = 321,
  eKaonMinus = -eKaonPlus,
  eNeutron = 2112,
  eAntiNeutron = -eNeutron,
  eProton = 2212,
  eAntiProton = -eProton,
  eLead = 1000822080
};

/// Convert an anti-particle to its particle and leave particles as-is.
static constexpr PdgParticle makeAbsolutePdgParticle(PdgParticle pdg) {
  const auto value = static_cast<std::int32_t>(pdg);
  return static_cast<PdgParticle>((0 <= value) ? value : -value);
}

/// Check if the PDG belongs to a nucleus, i.e. if it has 10 digits.
/// See PDG section "Monte Carlo Particle Numbering Scheme", point 16:
/// https://pdg.lbl.gov/2025/reviews/rpp2024-rev-monte-carlo-numbering.pdf
static constexpr bool isNucleus(PdgParticle pdg) {
  const auto pdgNum = static_cast<std::int32_t>(pdg);
  return std::abs(pdgNum) > 1e9;
}

/// Convert an excited nucleus to its ground state. PDG number of a nucleus has
/// a form 10LZZZAAAI, where I is isomer level; I=0 is the ground state.
/// See PDG section "Monte Carlo Particle Numbering Scheme", point 16:
/// https://pdg.lbl.gov/2025/reviews/rpp2024-rev-monte-carlo-numbering.pdf
static constexpr PdgParticle makeNucleusGroundState(PdgParticle pdg) {
  if (!isNucleus(pdg)) {
    throw std::invalid_argument("PDG must represent a nucleus");
  }
  const auto pdgNum = static_cast<std::int32_t>(pdg);
  // set isomer level to zero
  return static_cast<PdgParticle>((pdgNum / 10) * 10);
}

/// Extract Z and A for a given nucleus. PDG number of a nucleus has a form
/// 10LZZZAAAI, where L is number of lambdas, ZZZ is proton number, AAA is
/// atomic number, I is isomer level. See PDG section "Monte Carlo Particle
/// Numbering Scheme" , point 16:
/// https://pdg.lbl.gov/2025/reviews/rpp2024-rev-monte-carlo-numbering.pdf
static constexpr std::pair<std::int32_t, std::int32_t> extractNucleusZandA(
    PdgParticle pdg) {
  if (!isNucleus(pdg)) {
    throw std::invalid_argument("PDG must represent a nucleus");
  }
  const auto pdgNum = static_cast<std::int32_t>(pdg);
  // proton number respects the charge
  int Z = (pdgNum / 10000) % 1000;
  // atomic number is always positive
  int A = std::abs((pdgNum / 10) % 1000);
  return std::make_pair(Z, A);
}

}  // namespace Acts
