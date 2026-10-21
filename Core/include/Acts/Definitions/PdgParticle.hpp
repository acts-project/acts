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
#include <string>
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
  eKaon0Short = 310,
  eLambda0 = 3122,
  eJPsi = 443,   // J/ψ
  eB0 = 511,     // B0 meson (bd)
  eBPlus = 521,  // B+ meson (bu)
  eD0 = 421,     // D0 meson (cu)
  eDPlus = 411,  // D+ meson (cd)
  eAntiB0 = -eB0,
  eAntiD0 = -eD0,
  eNeutrinoE = 12,    // electron neutrino
  eNeutrinoMu = 14,   // muon neutrino
  eNeutrinoTau = 16,  // tau neutrino
  eAntiNeutrinoE = -eNeutrinoE,
  eAntiNeutrinoMu = -eNeutrinoMu,
  eAntiNeutrinoTau = -eNeutrinoTau
};

/// Convert an anti-particle to its particle and leave particles as-is.
/// @param pdg The PDG particle code
/// @return The absolute PDG particle code
static constexpr PdgParticle makeAbsolutePdgParticle(PdgParticle pdg) {
  const auto value = static_cast<std::int32_t>(pdg);
  return static_cast<PdgParticle>((0 <= value) ? value : -value);
}

/// Check if the PDG belongs to a nucleus, i.e. if it has 10 digits.
/// See PDG section "Monte Carlo Particle Numbering Scheme", point 16:
/// https://pdg.lbl.gov/2025/reviews/rpp2024-rev-monte-carlo-numbering.pdf
/// @param pdg The PDG particle code
/// @return True if the PDG code represents a nucleus
static constexpr bool isNucleus(PdgParticle pdg) {
  const auto value = static_cast<std::int32_t>(pdg);
  return std::abs(value) > 1e9;
}

/// Convert an excited nucleus to its ground state. PDG number of a nucleus has
/// a form 10LZZZAAAI, where I is isomer level; I=0 is the ground state.
/// See PDG section "Monte Carlo Particle Numbering Scheme", point 16:
/// https://pdg.lbl.gov/2025/reviews/rpp2024-rev-monte-carlo-numbering.pdf
/// @param pdg The PDG particle code of the nucleus
/// @return The PDG particle code of the nucleus in its ground state
static constexpr PdgParticle makeNucleusGroundState(PdgParticle pdg) {
  if (!isNucleus(pdg)) {
    throw std::invalid_argument("PDG must represent a nucleus");
  }
  // set isomer level to zero
  const auto value = static_cast<std::int32_t>(pdg);
  return static_cast<PdgParticle>((value / 10) * 10);
}

/// Extract Z and A for a given nucleus. PDG number of a nucleus has a form
/// 10LZZZAAAI, where L is number of lambdas, ZZZ is proton number, AAA is
/// atomic number, I is isomer level. See PDG section "Monte Carlo Particle
/// Numbering Scheme" , point 16:
/// https://pdg.lbl.gov/2025/reviews/rpp2024-rev-monte-carlo-numbering.pdf
/// @param pdg The PDG particle code of the nucleus
/// @return A pair containing the proton number (Z) and atomic number (A)
static constexpr std::pair<std::int32_t, std::int32_t> extractNucleusZandA(
    PdgParticle pdg) {
  if (!isNucleus(pdg)) {
    throw std::invalid_argument("PDG must represent a nucleus");
  }
  const auto value = static_cast<std::int32_t>(pdg);
  // proton number respects the charge
  int Z = (value / 10000) % 1000;
  // atomic number is always positive
  int A = std::abs((value / 10) % 1000);
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

/// Stream operator for HadronType
/// @param os Output stream
/// @param hadron The hadron type to output
/// @return Reference to output stream
std::ostream& operator<<(std::ostream& os, HadronType hadron);

/// Parse a PdgParticle from a particle name string.
/// Supports common particle names like "e-", "e+", "mu-", "mu+", "tau-",
/// "tau+", "gamma", "pi0", "pi+", "pi-", "K+", "K-", "n", "n~", "p", "p~",
/// "Pb".
/// @param name The particle name string
/// @return The corresponding PdgParticle enum value
/// @throws std::invalid_argument if the name is not recognized
PdgParticle parsePdgParticle(const std::string& name);

}  // namespace Acts
