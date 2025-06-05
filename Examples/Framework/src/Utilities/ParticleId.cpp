// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/ParticleId.hpp"

#include <cmath>
#include <ostream>
#include <utility>
#include <vector>

// This file is based on the following files, with modifications:
//
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h
//
// These files are licensed under Apache License 2.0

namespace ActsExamples::ParticleId {

namespace {

class DecodedPID : public std::pair<int, std::vector<int>> {
 public:
  explicit DecodedPID(const int& p) {
    this->first = p;
    this->second.reserve(10);
    int ap = std::abs(p);
    for (; ap; ap /= 10) {
      this->second.push_back(ap % 10);
    }
    std::reverse(this->second.begin(), this->second.end());
  }
  inline DecodedPID shift(const std::size_t n) const {
    return DecodedPID(this->first %
                      static_cast<int>(std::pow(10, ndigits() - n)));
  }
  inline const int& operator()(const std::size_t n) const {
    return this->second.at(n);
  }
  inline const int& last() const { return this->second.back(); }
  inline const int& pid() const { return this->first; }
  inline int max_digit(const int m, const int n) const {
    return *std::max_element(second.rbegin() + m, second.rbegin() + n);
  }
  inline int min_digit(const int m, const int n) const {
    return *std::min_element(second.rbegin() + m, second.rbegin() + n);
  }
  inline std::size_t ndigits() const { return this->second.size(); }
};

static constexpr int DQUARK = 1;
static constexpr int UQUARK = 2;
static constexpr int SQUARK = 3;
static constexpr int CQUARK = 4;
static constexpr int BQUARK = 5;
static constexpr int TQUARK = 6;

// static constexpr int ELECTRON = 11;
// static constexpr int POSITRON = -ELECTRON;
// static constexpr int NU_E = 12;
// static constexpr int MUON = 13;
// static constexpr int NU_MU = 14;
// static constexpr int TAU = 15;
// static constexpr int NU_TAU = 16;

// static constexpr int GLUON = 21;
// APID: 9 rather than 21 is used to denote a gluon/gluino in composite states.
// (From PDG 11g)
// static constexpr int COMPOSITEGLUON = 9;
// static constexpr int PHOTON = 22;
// static constexpr int Z0BOSON = 23;
// static constexpr int WPLUSBOSON = 24;
// static constexpr int HIGGSBOSON = 25;
// static constexpr int ZPRIME = 32;         // Z′/Z^0_2
// static constexpr int ZDBLPRIME = 33;      // Z′′/Z^0_3
// static constexpr int WPLUSPRIME = 34;     // W ′/W^+_2
// static constexpr int HIGGS2 = 35;         // H^0/H^0_2  FIXME Any better
// ideas? static constexpr int HIGGS3 = 36;         // A^0/H^0_3 FIXME Any
// better ideas? static constexpr int HIGGSPLUS = 37;      // H^+ static
// constexpr int HIGGSPLUSPLUS = 38;  // H^++ static constexpr int GRAVITON =
// 39; static constexpr int HIGGS4 = 40;  // a^0/H^0_4 FIXME Any better ideas?
// static constexpr int LEPTOQUARK = 42;

/// PDG Ids for Mavtop madgraph UFO model found under DarkX. The
/// mavtop is a vector-like top partner with coupling to a dark photon.
/// Theory paper: https://arxiv.org/abs/1904.05893
/// Pheno paper: https://arxiv.org/pdf/2112.08425
// static constexpr int DARKPHOTON = 60000;
static constexpr int MAVTOP = 60001;

// static constexpr int PIPLUS = 211;
// static constexpr int PIMINUS = -PIPLUS;
// static constexpr int PI0 = 111;
static constexpr int K0L = 130;

static constexpr int K0S = 310;
static constexpr int K0 = 311;
// static constexpr int KPLUS = 321;
// static constexpr int DPLUS = 411;
// static constexpr int DSTAR = 413;
// static constexpr int D0 = 421;
// static constexpr int DSPLUS = 431;
// static constexpr int JPSI = 443;
// static constexpr int B0 = 511;
// static constexpr int BCPLUS = 541;
// static constexpr int PROTON = 2212;
// static constexpr int NEUTRON = 2112;
// static constexpr int LAMBDA0 = 3122;
// static constexpr int LAMBDACPLUS = 4122;
// static constexpr int LAMBDAB0 = 5122;
// static constexpr int PSI2S = 20443;

/// PDG Rule 12:
/// Generator defined PDG ID values for right handed neutrinos and
/// corresponding W+ boson from a Left-Right symmetric Standard Model
/// extension. (Defined for some MadGraph+Pythia8 samples and
/// referenced in MCTruthClassifierGen.cxx)
// static constexpr int RH_NU_E = 9900012;
// static constexpr int RH_NU_MU = 9900014;
// static constexpr int RH_NU_TAU = 9900016;
// static constexpr int WBOSON_LRSM = 9900024;

// static constexpr int LEAD = 1000822080;
// static constexpr int OXYGEN = 1000080160;
// static constexpr int NEON = 1000100200;

/// PDG rule 8:
/// The pomeron and odderon trajectories and a generic reggeon trajectory
/// of states in QCD areassigned codes 990, 9990, and 110 respectively
// static constexpr int POMERON = 990;
// static constexpr int ODDERON = 9990;
// static constexpr int REGGEON = 110;

/// PDG rule 10:
/// Codes 81–100 are reserved for generator-specific pseudoparticles and
/// concepts. Codes 901–930, 1901–1930, 2901–2930, and 3901–3930 are for
/// additional components of Standard Modelparton distribution functions, where
/// the latter three ranges are intended to distinguish left/right/ longitudinal
/// components. Codes 998 and 999 are reserved for GEANT tracking pur-poses.
// static constexpr int GEANTINOPLUS = 998;
// static constexpr int GEANTINO0 = 999;

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L227-253
bool isMeson(const DecodedPID& p) {
  if (p.ndigits() < 3) {
    return false;
  }
  if (p.ndigits() == 7 && (p(0) == 1 || p(0) == 2)) {
    return false;  // APID don't match SUSY particles
  }
  if (std::abs(p.pid()) == K0S) {
    return true;
  }
  if (std::abs(p.pid()) == K0L) {
    return true;
  }
  if (std::abs(p.pid()) == K0) {
    return true;
  }
  if (p.last() % 2 != 1) {
    return false;
  }
  if (p.max_digit(1, 3) >= 6) {
    return false;
  }
  if (p.max_digit(1, 3) == 0) {
    return false;
  }
  if (p.ndigits() > 3 && *(p.second.rbegin() + 3) != 0) {
    return false;
  }

  if (p.ndigits() == 3 && p(0) == p(1) && p.pid() < 0) {
    return false;
  }
  if (p.ndigits() == 5 && p(2) == p(3) && p.pid() < 0) {
    return false;
  }
  if (p.ndigits() == 7 && p(4) == p(5) && p.pid() < 0) {
    return false;
  }

  if (p.ndigits() == 3 && p(0) >= p(1) && p(1) != 0) {
    return true;
  }
  if (p.ndigits() == 5 && p(2) >= p(3) && p(3) != 0 && p(0) == 1 && p(1) == 0) {
    return true;
  }
  if (p.ndigits() == 5 && p(2) >= p(3) && p(3) != 0 && p(0) == 2 && p(1) == 0 &&
      p.last() > 1) {
    return true;
  }
  if (p.ndigits() == 5 && p(2) >= p(3) && p(3) != 0 && p(0) == 3 && p(1) == 0 &&
      p.last() > 1) {
    return true;
  }

  if (p.ndigits() == 6 && p(3) >= p(4) && p(4) != 0 && p.last() % 2 == 1) {
    return true;
  }

  if (p.ndigits() == 7 && p(0) == 9 && p(1) == 0 && p(4) >= p(5) && p(5) != 0) {
    return true;
  }

  return false;
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L258-291
bool isBaryon(const DecodedPID& p) {
  if (p.ndigits() < 4) {
    return false;
  }
  if (p.max_digit(1, 4) >= 6) {
    return false;
  }
  if (p.min_digit(1, 4) == 0) {
    return false;
  }
  if (p.ndigits() == 4 &&
      (p.last() == 2 || p.last() == 4 || p.last() == 6 || p.last() == 8)) {
    return true;
  }

  if (p.ndigits() == 5 && p(0) == 1 && (p.last() == 2 || p.last() == 4)) {
    return true;
  }
  if (p.ndigits() == 5 && p(0) == 3 && (p.last() == 2 || p.last() == 4)) {
    return true;
  }

  if (p.ndigits() == 6) {
    if (p(0) == 1 && p(1) == 0 && p.last() == 2) {
      return true;
    }
    if (p(0) == 1 && p(1) == 1 && p.last() == 2) {
      return true;
    }
    if (p(0) == 1 && p(1) == 2 && p.last() == 4) {
      return true;
    }

    if (p(0) == 2 && p(1) == 0 && p.last() == 2) {
      return true;
    }
    if (p(0) == 2 && p(1) == 0 && p.last() == 4) {
      return true;
    }
    if (p(0) == 2 && p(1) == 1 && p.last() == 2) {
      return true;
    }

    if (p(0) == 1 && p(1) == 0 && p.last() == 4) {
      return true;
    }
    if (p(0) == 1 && p(1) == 0 && p.last() == 6) {
      return true;
    }
    if (p(0) == 2 && p(1) == 0 && p.last() == 6) {
      return true;
    }
    if (p(0) == 2 && p(1) == 0 && p.last() == 8) {
      return true;
    }
  }

  if (p.ndigits() == 5) {
    if (p(0) == 2 && p.last() == 2) {
      return true;
    }
    if (p(0) == 2 && p.last() == 4) {
      return true;
    }
    if (p(0) == 2 && p.last() == 6) {
      return true;
    }
    if (p(0) == 5 && p.last() == 2) {
      return true;
    }
    if (p(0) == 1 && p.last() == 6) {
      return true;
    }
    if (p(0) == 4 && p.last() == 2) {
      return true;
    }
  }
  return false;
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L301-307
bool isTetraquark(const DecodedPID& p) {
  return (p.ndigits() == 9 && p(0) == 1 && p(5) == 0 &&
          p.max_digit(1, 3) <= 6 && p.min_digit(1, 3) > 0 &&
          p.max_digit(1 + 3, 3 + 3) <= 6 && p.min_digit(1 + 3, 3 + 3) > 0 &&
          (p(3) >= p(4) && p(6) >= p(7)) &&
          ((p(3) > p(6)) || (p(3) == p(6) && (p(4) >= p(7)))));
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L316-320
bool isPentaquark(const DecodedPID& p) {
  return (p.ndigits() == 9 && p(0) == 1 && p.max_digit(1, 6) <= 6 &&
          p.min_digit(1, 6) > 0 &&
          (p(3) >= p(4) && p(4) >= p(5) && p(5) >= p(6)));
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L211-216
bool isDiquark(const DecodedPID& p) {
  if (p.ndigits() == 4 && p(0) >= p(1) && p(2) == 0 && p.last() % 2 == 1 &&
      p.max_digit(2, 4) <= TQUARK) {
    return true;
  }
  return false;
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L399-407
bool isGenSpecific(const int& p) {
  if (p >= 81 && p <= 100) {
    return true;
  }
  if (p >= 901 && p <= 930) {
    return true;
  }
  if (p >= 998 && p <= 999) {
    return true;
  }
  if (p >= 1901 && p <= 1930) {
    return true;
  }
  if (p >= 2901 && p <= 2930) {
    return true;
  }
  if (p >= 3901 && p <= 3930) {
    return true;
  }
  return false;
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L429
bool isSUSY(const DecodedPID& p) {
  return (p.ndigits() == 7 && (p(0) == 1 || p(0) == 2) &&
          !isGenSpecific(p.shift(2).pid()));
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L824-838
int leadingQuark(const DecodedPID& p) {
  if (isQuark(p.pid())) {
    return std::abs(p.pid());
  }
  if (isMeson(p)) {
    return p.max_digit(1, 3);
  }
  if (isDiquark(p)) {
    return p.max_digit(2, 4);
  }
  if (isBaryon(p)) {
    return p.max_digit(1, 4);
  }
  if (isTetraquark(p)) {
    return p.max_digit(1, 5);
  }
  if (isPentaquark(p)) {
    return p.max_digit(1, 6);
  }
  if (isSUSY(p)) {  // APID SUSY case
    auto pp = p.shift(1);
    if (pp.ndigits() == 1) {
      return 0;
    }  // Handle squarks
    if (pp.ndigits() == 3) {
      pp = DecodedPID(pp(1));
    }  // Handle ~q qbar pairs
    if (pp.ndigits() > 3) {
      pp = pp.shift(1);
    }  // Drop gluinos and squarks
    return leadingQuark(pp);
  }
  return 0;
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L853
bool isBottomMeson(const DecodedPID& p) {
  return leadingQuark(p) == BQUARK && isMeson(p);
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L861
bool isBBbarMeson(const DecodedPID& p) {
  return leadingQuark(p) == BQUARK && isMeson(p) &&
         (*(p.second.rbegin() + 2)) == BQUARK &&
         (*(p.second.rbegin() + 1)) == BQUARK;
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L857
bool isCCbarMeson(const DecodedPID& p) {
  return leadingQuark(p) == CQUARK && isMeson(p) &&
         (*(p.second.rbegin() + 2)) == CQUARK &&
         (*(p.second.rbegin() + 1)) == CQUARK;
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L852
bool isCharmMeson(const DecodedPID& p) {
  return leadingQuark(p) == CQUARK && isMeson(p);
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L865-870
bool isLightBaryon(const DecodedPID& p) {
  auto lq = leadingQuark(p);
  return (lq == DQUARK || lq == UQUARK || lq == SQUARK) && isBaryon(p);
}

// bool isHeavyBaryon(const DecodedPID& p) {
//   auto lq = leadingQuark(p);
//   return (lq == CQUARK || lq == BQUARK || lq == TQUARK) && isBaryon(p);
// }

bool isStrangeBaryon(const DecodedPID& p) {
  return leadingQuark(p) == SQUARK && isBaryon(p);
}

bool isCharmBaryon(const DecodedPID& p) {
  return leadingQuark(p) == CQUARK && isBaryon(p);
}

bool isBottomBaryon(const DecodedPID& p) {
  return leadingQuark(p) == BQUARK && isBaryon(p);
}

// bool isTopBaryon(const DecodedPID& p) {
//   return leadingQuark(p) == TQUARK && isBaryon(p);
// }

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L851
bool isStrangeMeson(const DecodedPID& p) {
  return leadingQuark(p) == SQUARK && isMeson(p);
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L849
bool isLightMeson(const DecodedPID& p) {
  auto lq = leadingQuark(p);
  return (lq == DQUARK || lq == UQUARK || lq == SQUARK) && isMeson(p);
}

}  // namespace

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L325
bool isHadron(int pdg) {
  DecodedPID p(pdg);
  return isMeson(p) || isBaryon(p) || isTetraquark(p) || isPentaquark(p);
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L180
bool isLepton(int pdg) {
  auto sp = std::abs(pdg);
  return sp >= 11 && sp <= 18;
}

HadronType hadronType(int pdg) {
  DecodedPID p(pdg);

  using enum HadronType;

  if (isBBbarMeson(p)) {
    return BBbarMeson;
  }
  if (isCCbarMeson(p)) {
    return CCbarMeson;
  }
  if (isBottomMeson(p)) {
    return BottomMeson;
  }
  if (isCharmMeson(p)) {
    return CharmedMeson;
  }
  if (isBottomBaryon(p)) {
    return BottomBaryon;
  }
  if (isCharmBaryon(p)) {
    return CharmedBaryon;
  }
  if (isStrangeBaryon(p)) {
    return StrangeBaryon;
  }
  if (isLightBaryon(p)) {
    return LightBaryon;
  }
  if (isStrangeMeson(p)) {
    return StrangeMeson;
  }
  if (isLightMeson(p)) {
    return LightMeson;
  }

  return Unknown;
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L159
bool isQuark(int pdg) {
  return pdg != 0 && (std::abs(pdg) <= 8 || std::abs(pdg) == MAVTOP);
}

std::ostream& operator<<(std::ostream& os, HadronType hadron) {
  switch (hadron) {
    using enum HadronType;
    case Hadron:
      return os << "Hadron";
    case BBbarMeson:
      return os << "BBbarMeson";
    case CCbarMeson:
      return os << "CCbarMeson";
    case BottomMeson:
      return os << "BottomMeson";
    case BottomBaryon:
      return os << "BottomBaryon";
    case CharmedMeson:
      return os << "CharmedMeson";
    case CharmedBaryon:
      return os << "CharmedBaryon";
    case StrangeMeson:
      return os << "StrangeMeson";
    case StrangeBaryon:
      return os << "StrangeBaryon";
    case LightMeson:
      return os << "LightMeson";
    case LightBaryon:
      return os << "LightBaryon";
    case Unknown:
      return os << "Unknown";
  }
}
}  // namespace ActsExamples::ParticleId
