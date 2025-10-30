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
#include <iosfwd>
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
  eLead = 1000822080,
  eUpsilon = 553,        // ϒ
  eUpsilon1S = 553,     // ϒ(1S)
  eUpsilon2S = 100553,  // ϒ(2S)
  eUpsilon3S = 200553,  // ϒ(3S)
  eChiB = 10551,        // χb0
  eJPsi = 443,          // J/ψ
  ePsi2S = 100443,     // ψ(2S)
  eChiC = 10441,       // χc0
  eEta_c = 441,        // ηc
  eB0 = 511,           // B0 meson (bd)
  eBPlus = 521,       // B+ meson (bu)
  eD0 = 421,          // D0 meson (cu)
  eDPlus = 411,      // D+ meson (cd)
  eAntiB0 = -eB0,
  eAntiD0 = -eD0,
  eLambdaB = 5122,  // Λb0 (udb)
  eSigmaB = 5222,   // Σb+ (uub)
  eXiB = 5322,   // Ξb0 (usb)
  eLambdaC = 4122,  // Λc+ (udc)
  eSigmaC = 4222,   // Σc++ (uuc)
  eXiC = 4322,      // Ξc0 (usc)
  eAntiLambdaB = -eLambdaB,
  eAntiLambdaC = -eLambdaC,
  eNeutrinoE = 12,    // electron neutrino
  eNeutrinoMu = 14,   // muon neutrino
  eNeutrinoTau = 16,  // tau neutrino
  eAntiNeutrinoE = -eNeutrinoE,
  eAntiNeutrinoMu = -eNeutrinoMu,
  eAntiNeutrinoTau = -eNeutrinoTau
};

/// Convert an anti-particle to its particle and leave particles as-is.
static constexpr PdgParticle makeAbsolutePdgParticle(PdgParticle pdg) {
  const auto value = static_cast<std::int32_t>(pdg);
  return static_cast<PdgParticle>((0 <= value) ? value : -value);
}

/// Check if the PDG belongs to a nucleus, i.e. if it has 10 digits.
/// See PDG section "Monte Carlo Particle Numbering Scheme", point 16:
/// https://pdg.lbl.gov/2025/reviews/rpp2024-rev-monte-carlo-numbering.pdf
static constexpr bool isNucleus(std::int32_t pdg) {
  return std::abs(pdg) > 1e9;
}

/// Convert an excited nucleus to its ground state. PDG number of a nucleus has
/// a form 10LZZZAAAI, where I is isomer level; I=0 is the ground state.
/// See PDG section "Monte Carlo Particle Numbering Scheme", point 16:
/// https://pdg.lbl.gov/2025/reviews/rpp2024-rev-monte-carlo-numbering.pdf
static constexpr PdgParticle makeNucleusGroundState(std::int32_t pdg) {
  if (!isNucleus(pdg)) {
    throw std::invalid_argument("PDG must represent a nucleus");
  }
  // set isomer level to zero
  return static_cast<PdgParticle>((pdg / 10) * 10);
}

/// Extract Z and A for a given nucleus. PDG number of a nucleus has a form
/// 10LZZZAAAI, where L is number of lambdas, ZZZ is proton number, AAA is
/// atomic number, I is isomer level. See PDG section "Monte Carlo Particle
/// Numbering Scheme" , point 16:
/// https://pdg.lbl.gov/2025/reviews/rpp2024-rev-monte-carlo-numbering.pdf
static constexpr std::pair<std::int32_t, std::int32_t> extractNucleusZandA(
    std::int32_t pdg) {
  if (!isNucleus(pdg)) {
    throw std::invalid_argument("PDG must represent a nucleus");
  }
  // proton number respects the charge
  int Z = (pdg / 10000) % 1000;
  // atomic number is always positive
  int A = std::abs((pdg / 10) % 1000);
  return std::make_pair(Z, A);
}

/// Hadron type classification for B, C, strange and light hadrons.
enum class HadronType {
  Hadron = 1,
  BBbarMeson = 2,
  CCbarMeson = 3,
  BottomMeson = 4,
  BottomBaryon = 5,
  CharmedMeson = 6,
  CharmedBaryon = 7,
  StrangeMeson = 8,
  StrangeBaryon = 9,
  LightMeson = 10,
  LightBaryon = 11,
  Unknown = 12
};

std::ostream& operator<<(std::ostream& os, HadronType hadron);

}  // namespace Acts
