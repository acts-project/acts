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

// static constexpr int DQUARK = 1;
// static constexpr int UQUARK = 2;
// static constexpr int SQUARK = 3;
// static constexpr int CQUARK = 4;
// static constexpr int BQUARK = 5;
// static constexpr int TQUARK = 6;

// // static constexpr int ELECTRON = 11;
// // static constexpr int POSITRON = -ELECTRON;
// // static constexpr int NU_E = 12;
// // static constexpr int MUON = 13;
// // static constexpr int NU_MU = 14;
// // static constexpr int TAU = 15;
// // static constexpr int NU_TAU = 16;

// static constexpr int GLUON = 21;
// // APID: 9 rather than 21 is used to denote a gluon/gluino in composite
// states.
// // (From PDG 11g)
// static constexpr int COMPOSITEGLUON = 9;
// // static constexpr int PHOTON = 22;
// // static constexpr int Z0BOSON = 23;
// // static constexpr int WPLUSBOSON = 24;
// // static constexpr int HIGGSBOSON = 25;
// // static constexpr int ZPRIME = 32;         // Z′/Z^0_2
// // static constexpr int ZDBLPRIME = 33;      // Z′′/Z^0_3
// // static constexpr int WPLUSPRIME = 34;     // W ′/W^+_2
// // static constexpr int HIGGS2 = 35;         // H^0/H^0_2  FIXME Any better
// // ideas? static constexpr int HIGGS3 = 36;         // A^0/H^0_3 FIXME Any
// // better ideas? static constexpr int HIGGSPLUS = 37;      // H^+ static
// // constexpr int HIGGSPLUSPLUS = 38;  // H^++ static constexpr int GRAVITON =
// // 39; static constexpr int HIGGS4 = 40;  // a^0/H^0_4 FIXME Any better
// ideas? static constexpr int LEPTOQUARK = 42;

// /// PDG Ids for Mavtop madgraph UFO model found under DarkX. The
// /// mavtop is a vector-like top partner with coupling to a dark photon.
// /// Theory paper: https://arxiv.org/abs/1904.05893
// /// Pheno paper: https://arxiv.org/pdf/2112.08425
// // static constexpr int DARKPHOTON = 60000;
// static constexpr int MAVTOP = 60001;

// // static constexpr int PIPLUS = 211;
// // static constexpr int PIMINUS = -PIPLUS;
// // static constexpr int PI0 = 111;
// static constexpr int K0L = 130;

// static constexpr int K0S = 310;
// static constexpr int K0 = 311;
// // static constexpr int KPLUS = 321;
// // static constexpr int DPLUS = 411;
// // static constexpr int DSTAR = 413;
// // static constexpr int D0 = 421;
// // static constexpr int DSPLUS = 431;
// // static constexpr int JPSI = 443;
// // static constexpr int B0 = 511;
// // static constexpr int BCPLUS = 541;
// // static constexpr int PROTON = 2212;
// // static constexpr int NEUTRON = 2112;
// // static constexpr int LAMBDA0 = 3122;
// // static constexpr int LAMBDACPLUS = 4122;
// // static constexpr int LAMBDAB0 = 5122;
// // static constexpr int PSI2S = 20443;

// /// PDG Rule 12:
// /// Generator defined PDG ID values for right handed neutrinos and
// /// corresponding W+ boson from a Left-Right symmetric Standard Model
// /// extension. (Defined for some MadGraph+Pythia8 samples and
// /// referenced in MCTruthClassifierGen.cxx)
// // static constexpr int RH_NU_E = 9900012;
// // static constexpr int RH_NU_MU = 9900014;
// // static constexpr int RH_NU_TAU = 9900016;
// // static constexpr int WBOSON_LRSM = 9900024;

// // static constexpr int LEAD = 1000822080;
// // static constexpr int OXYGEN = 1000080160;
// // static constexpr int NEON = 1000100200;

// /// PDG rule 8:
// /// The pomeron and odderon trajectories and a generic reggeon trajectory
// /// of states in QCD areassigned codes 990, 9990, and 110 respectively
// // static constexpr int POMERON = 990;
// // static constexpr int ODDERON = 9990;
// // static constexpr int REGGEON = 110;

// /// PDG rule 10:
// /// Codes 81–100 are reserved for generator-specific pseudoparticles and
// /// concepts. Codes 901–930, 1901–1930, 2901–2930, and 3901–3930 are for
// /// additional components of Standard Modelparton distribution functions,
// where
// /// the latter three ranges are intended to distinguish left/right/
// longitudinal
// /// components. Codes 998 and 999 are reserved for GEANT tracking pur-poses.
// // static constexpr int GEANTINOPLUS = 998;
// // static constexpr int GEANTINO0 = 999;

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L227-253
// bool isMeson(const DecodedPID& p) {
//   if (p.ndigits() < 3) {
//     return false;
//   }
//   if (p.ndigits() == 7 && (p(0) == 1 || p(0) == 2)) {
//     return false;  // APID don't match SUSY particles
//   }
//   if (std::abs(p.pid()) == K0S) {
//     return true;
//   }
//   if (std::abs(p.pid()) == K0L) {
//     return true;
//   }
//   if (std::abs(p.pid()) == K0) {
//     return true;
//   }
//   if (p.last() % 2 != 1) {
//     return false;
//   }
//   if (p.max_digit(1, 3) >= 6) {
//     return false;
//   }
//   if (p.max_digit(1, 3) == 0) {
//     return false;
//   }
//   if (p.ndigits() > 3 && *(p.second.rbegin() + 3) != 0) {
//     return false;
//   }

//   if (p.ndigits() == 3 && p(0) == p(1) && p.pid() < 0) {
//     return false;
//   }
//   if (p.ndigits() == 5 && p(2) == p(3) && p.pid() < 0) {
//     return false;
//   }
//   if (p.ndigits() == 7 && p(4) == p(5) && p.pid() < 0) {
//     return false;
//   }

//   if (p.ndigits() == 3 && p(0) >= p(1) && p(1) != 0) {
//     return true;
//   }
//   if (p.ndigits() == 5 && p(2) >= p(3) && p(3) != 0 && p(0) == 1 && p(1) ==
//   0) {
//     return true;
//   }
//   if (p.ndigits() == 5 && p(2) >= p(3) && p(3) != 0 && p(0) == 2 && p(1) == 0
//   &&
//       p.last() > 1) {
//     return true;
//   }
//   if (p.ndigits() == 5 && p(2) >= p(3) && p(3) != 0 && p(0) == 3 && p(1) == 0
//   &&
//       p.last() > 1) {
//     return true;
//   }

//   if (p.ndigits() == 6 && p(3) >= p(4) && p(4) != 0 && p.last() % 2 == 1) {
//     return true;
//   }

//   if (p.ndigits() == 7 && p(0) == 9 && p(1) == 0 && p(4) >= p(5) && p(5) !=
//   0) {
//     return true;
//   }

//   return false;
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L258-291
// bool isBaryon(const DecodedPID& p) {
//   if (p.ndigits() < 4) {
//     return false;
//   }
//   if (p.max_digit(1, 4) >= 6) {
//     return false;
//   }
//   if (p.min_digit(1, 4) == 0) {
//     return false;
//   }
//   if (p.ndigits() == 4 &&
//       (p.last() == 2 || p.last() == 4 || p.last() == 6 || p.last() == 8)) {
//     return true;
//   }

//   if (p.ndigits() == 5 && p(0) == 1 && (p.last() == 2 || p.last() == 4)) {
//     return true;
//   }
//   if (p.ndigits() == 5 && p(0) == 3 && (p.last() == 2 || p.last() == 4)) {
//     return true;
//   }

//   if (p.ndigits() == 6) {
//     if (p(0) == 1 && p(1) == 0 && p.last() == 2) {
//       return true;
//     }
//     if (p(0) == 1 && p(1) == 1 && p.last() == 2) {
//       return true;
//     }
//     if (p(0) == 1 && p(1) == 2 && p.last() == 4) {
//       return true;
//     }

//     if (p(0) == 2 && p(1) == 0 && p.last() == 2) {
//       return true;
//     }
//     if (p(0) == 2 && p(1) == 0 && p.last() == 4) {
//       return true;
//     }
//     if (p(0) == 2 && p(1) == 1 && p.last() == 2) {
//       return true;
//     }

//     if (p(0) == 1 && p(1) == 0 && p.last() == 4) {
//       return true;
//     }
//     if (p(0) == 1 && p(1) == 0 && p.last() == 6) {
//       return true;
//     }
//     if (p(0) == 2 && p(1) == 0 && p.last() == 6) {
//       return true;
//     }
//     if (p(0) == 2 && p(1) == 0 && p.last() == 8) {
//       return true;
//     }
//   }

//   if (p.ndigits() == 5) {
//     if (p(0) == 2 && p.last() == 2) {
//       return true;
//     }
//     if (p(0) == 2 && p.last() == 4) {
//       return true;
//     }
//     if (p(0) == 2 && p.last() == 6) {
//       return true;
//     }
//     if (p(0) == 5 && p.last() == 2) {
//       return true;
//     }
//     if (p(0) == 1 && p.last() == 6) {
//       return true;
//     }
//     if (p(0) == 4 && p.last() == 2) {
//       return true;
//     }
//   }
//   return false;
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L301-307
// bool isTetraquark(const DecodedPID& p) {
//   return (p.ndigits() == 9 && p(0) == 1 && p(5) == 0 &&
//           p.max_digit(1, 3) <= 6 && p.min_digit(1, 3) > 0 &&
//           p.max_digit(1 + 3, 3 + 3) <= 6 && p.min_digit(1 + 3, 3 + 3) > 0 &&
//           (p(3) >= p(4) && p(6) >= p(7)) &&
//           ((p(3) > p(6)) || (p(3) == p(6) && (p(4) >= p(7)))));
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L316-320
// bool isPentaquark(const DecodedPID& p) {
//   return (p.ndigits() == 9 && p(0) == 1 && p.max_digit(1, 6) <= 6 &&
//           p.min_digit(1, 6) > 0 &&
//           (p(3) >= p(4) && p(4) >= p(5) && p(5) >= p(6)));
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L211-216
// bool isDiquark(const DecodedPID& p) {
//   if (p.ndigits() == 4 && p(0) >= p(1) && p(2) == 0 && p.last() % 2 == 1 &&
//       p.max_digit(2, 4) <= TQUARK) {
//     return true;
//   }
//   return false;
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L399-407
// bool isGenSpecific(const int& p) {
//   if (p >= 81 && p <= 100) {
//     return true;
//   }
//   if (p >= 901 && p <= 930) {
//     return true;
//   }
//   if (p >= 998 && p <= 999) {
//     return true;
//   }
//   if (p >= 1901 && p <= 1930) {
//     return true;
//   }
//   if (p >= 2901 && p <= 2930) {
//     return true;
//   }
//   if (p >= 3901 && p <= 3930) {
//     return true;
//   }
//   return false;
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L429
// bool isSUSY(const DecodedPID& p) {
//   return (p.ndigits() == 7 && (p(0) == 1 || p(0) == 2) &&
//           !isGenSpecific(p.shift(2).pid()));
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L824-838
// int leadingQuark(const DecodedPID& p) {
//   if (isQuark(p.pid())) {
//     return std::abs(p.pid());
//   }
//   if (isMeson(p)) {
//     return p.max_digit(1, 3);
//   }
//   if (isDiquark(p)) {
//     return p.max_digit(2, 4);
//   }
//   if (isBaryon(p)) {
//     return p.max_digit(1, 4);
//   }
//   if (isTetraquark(p)) {
//     return p.max_digit(1, 5);
//   }
//   if (isPentaquark(p)) {
//     return p.max_digit(1, 6);
//   }
//   if (isSUSY(p)) {  // APID SUSY case
//     auto pp = p.shift(1);
//     if (pp.ndigits() == 1) {
//       return 0;
//     }  // Handle squarks
//     if (pp.ndigits() == 3) {
//       pp = DecodedPID(pp(1));
//     }  // Handle ~q qbar pairs
//     if (pp.ndigits() > 3) {
//       pp = pp.shift(1);
//     }  // Drop gluinos and squarks
//     return leadingQuark(pp);
//   }
//   return 0;
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L853
// bool isBottomMeson(const DecodedPID& p) {
//   return leadingQuark(p) == BQUARK && isMeson(p);
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L861
// bool isBBbarMeson(const DecodedPID& p) {
//   return leadingQuark(p) == BQUARK && isMeson(p) &&
//          (*(p.second.rbegin() + 2)) == BQUARK &&
//          (*(p.second.rbegin() + 1)) == BQUARK;
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L857
// bool isCCbarMeson(const DecodedPID& p) {
//   return leadingQuark(p) == CQUARK && isMeson(p) &&
//          (*(p.second.rbegin() + 2)) == CQUARK &&
//          (*(p.second.rbegin() + 1)) == CQUARK;
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L852
// bool isCharmMeson(const DecodedPID& p) {
//   return leadingQuark(p) == CQUARK && isMeson(p);
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L865-870
// bool isLightBaryon(const DecodedPID& p) {
//   auto lq = leadingQuark(p);
//   return (lq == DQUARK || lq == UQUARK || lq == SQUARK) && isBaryon(p);
// }

// // bool isHeavyBaryon(const DecodedPID& p) {
// //   auto lq = leadingQuark(p);
// //   return (lq == CQUARK || lq == BQUARK || lq == TQUARK) && isBaryon(p);
// // }

// bool isStrangeBaryon(const DecodedPID& p) {
//   return leadingQuark(p) == SQUARK && isBaryon(p);
// }

// bool isCharmBaryon(const DecodedPID& p) {
//   return leadingQuark(p) == CQUARK && isBaryon(p);
// }

// bool isBottomBaryon(const DecodedPID& p) {
//   return leadingQuark(p) == BQUARK && isBaryon(p);
// }

// // bool isTopBaryon(const DecodedPID& p) {
// //   return leadingQuark(p) == TQUARK && isBaryon(p);
// // }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L851
// bool isStrangeMeson(const DecodedPID& p) {
//   return leadingQuark(p) == SQUARK && isMeson(p);
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L849
// bool isLightMeson(const DecodedPID& p) {
//   auto lq = leadingQuark(p);
//   return (lq == DQUARK || lq == UQUARK || lq == SQUARK) && isMeson(p);
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L347
// bool isGluon(const int& p) {
//   return p == GLUON;
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L414-419
// bool isGlueball(const DecodedPID& p) {
//   if (p.ndigits() > 4) {
//     return false;  // APID avoid classifying R-Glueballs as SM Glueballs
//   }
//   return ((p.ndigits() == 3 && p(0) == COMPOSITEGLUON &&
//            p(1) == COMPOSITEGLUON && (p(2) == 1 || p(2) == 5)) ||
//           (p.ndigits() == 4 && p(0) == COMPOSITEGLUON &&
//            p(1) == COMPOSITEGLUON && p(2) == COMPOSITEGLUON &&
//            (p(3) == 3 || p(3) == 7)));
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L382
// bool isLeptoQuark(int p) {
//   return std::abs(p) == LEPTOQUARK;
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L1048-1050
// bool isRHadron(const DecodedPID& p) {
//   return (isRBaryon(p) || isRMeson(p) || isRGlueball(p));
// }

// //
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L1080
// bool isStrongInteracting(const DecodedPID& p) {
//   return (isGluon(p.pid()) || isQuark(p.pid()) || isDiquark(p) ||
//           isGlueball(p) || isLeptoQuark(p.pid()) || isHadron(p.pid()) ||
//           isRHadron(p));
// }

static const int TABLESIZE = 100;
static const std::array<int, TABLESIZE> triple_charge = {
    +0, -1, +2, -1, +2, -1, +2, -1, +2, +0, +0, -3, +0, -3, +0, -3, +0,
    -3, +0, +0, +0, +0, +0, +0, +3, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +3, +0, +0, +3, +6, +0, +0, +0, -1, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0};
static const std::array<int, TABLESIZE> double_spin = {
    +0, +1, +1, +1, +1, +1, +1, +1, +1, +0, +0, +1, +1, +1, +1, +1, +1,
    +1, +1, +0, +2, +2, +2, +2, +2, +0, +0, +0, +0, +0, +0, +0, +2, +2,
    +2, +0, +0, +0, +0, +4, +0, +0, -1, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +1, +2, +0, +2, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0};

static const int DQUARK = 1;
static const int UQUARK = 2;
static const int SQUARK = 3;
static const int CQUARK = 4;
static const int BQUARK = 5;
static const int TQUARK = 6;

static const int ELECTRON = 11;
static const int POSITRON = -ELECTRON;
static const int NU_E = 12;
static const int MUON = 13;
static const int NU_MU = 14;
static const int TAU = 15;
static const int NU_TAU = 16;

static const int GLUON = 21;
// APID: 9 rather than 21 is used to denote a gluon/gluino in composite states.
// (From PDG 11g)
static const int COMPOSITEGLUON = 9;
static const int PHOTON = 22;
static const int Z0BOSON = 23;
static const int WPLUSBOSON = 24;
static const int HIGGSBOSON = 25;
static const int ZPRIME = 32;         // Z′/Z^0_2
static const int ZDBLPRIME = 33;      // Z′′/Z^0_3
static const int WPLUSPRIME = 34;     // W ′/W^+_2
static const int HIGGS2 = 35;         // H^0/H^0_2  FIXME Any better ideas?
static const int HIGGS3 = 36;         // A^0/H^0_3 FIXME Any better ideas?
static const int HIGGSPLUS = 37;      // H^+
static const int HIGGSPLUSPLUS = 38;  // H^++
static const int GRAVITON = 39;
static const int HIGGS4 = 40;  // a^0/H^0_4 FIXME Any better ideas?
static const int LEPTOQUARK = 42;

/// PDG Ids for Mavtop madgraph UFO model found under DarkX. The
/// mavtop is a vector-like top partner with coupling to a dark photon.
/// Theory paper: https://arxiv.org/abs/1904.05893
/// Pheno paper: https://arxiv.org/pdf/2112.08425
static const int DARKPHOTON = 60000;
static const int MAVTOP = 60001;

static const int PIPLUS = 211;
static const int PIMINUS = -PIPLUS;
static const int PI0 = 111;
static const int K0L = 130;

static const int K0S = 310;
static const int K0 = 311;
static const int KPLUS = 321;
static const int DPLUS = 411;
static const int DSTAR = 413;
static const int D0 = 421;
static const int DSPLUS = 431;
static const int JPSI = 443;
static const int B0 = 511;
static const int BCPLUS = 541;
static const int PROTON = 2212;
static const int NEUTRON = 2112;
static const int LAMBDA0 = 3122;
static const int LAMBDACPLUS = 4122;
static const int LAMBDAB0 = 5122;
static const int PSI2S = 20443;

/// PDG Rule 12:
/// Generator defined PDG ID values for right handed neutrinos and
/// corresponding W+ boson from a Left-Right symmetric Standard Model
/// extension. (Defined for some MadGraph+Pythia8 samples and
/// referenced in MCTruthClassifierGen.cxx)
static const int RH_NU_E = 9900012;
static const int RH_NU_MU = 9900014;
static const int RH_NU_TAU = 9900016;
static const int WBOSON_LRSM = 9900024;

static const int LEAD = 1000822080;
static const int OXYGEN = 1000080160;
static const int NEON = 1000100200;

/// PDG rule 8:
/// The pomeron and odderon trajectories and a generic reggeon trajectory
/// of states in QCD areassigned codes 990, 9990, and 110 respectively
static const int POMERON = 990;
static const int ODDERON = 9990;
static const int REGGEON = 110;

/// PDG rule 10:
/// Codes 81–100 are reserved for generator-specific pseudoparticles and
/// concepts. Codes 901–930, 1901–1930, 2901–2930, and 3901–3930 are for
/// additional components of Standard Modelparton distribution functions, where
/// the latter three ranges are intended to distinguish left/right/ longitudinal
/// components. Codes 998 and 999 are reserved for GEANT tracking pur-poses.
static const int GEANTINOPLUS = 998;
static const int GEANTINO0 = 999;

/// PDG rule 2:
/// Quarks and leptons are numbered consecutively starting from 1 and 11
/// respectively; to dothis they are first ordered by family and within
/// families by weak isospin.
/// APID: the fourth generation quarks are quarks.
template <class T>
inline bool isQuark(const T& p) {
  return isQuark(p->pdg_id());
}
template <>
inline bool isQuark(const int& p) {
  return p != 0 && (std::abs(p) <= 8 || std::abs(p) == MAVTOP);
}
template <>
inline bool isQuark(const DecodedPID& p) {
  return isQuark(p.pid());
}

template <class T>
inline bool isSMQuark(const T& p) {
  return isSMQuark(p->pdg_id());
}
template <>
inline bool isSMQuark(const int& p) {
  return p != 0 && std::abs(p) <= TQUARK;
}
template <>
inline bool isSMQuark(const DecodedPID& p) {
  return isSMQuark(p.pid());
}

template <class T>
inline bool isStrange(const T& p) {
  return isStrange(p->pdg_id());
}
template <>
inline bool isStrange(const int& p) {
  return std::abs(p) == 3;
}

template <class T>
inline bool isCharm(const T& p) {
  return isCharm(p->pdg_id());
}
template <>
inline bool isCharm(const int& p) {
  return std::abs(p) == 4;
}

template <class T>
inline bool isBottom(const T& p) {
  return isBottom(p->pdg_id());
}
template <>
inline bool isBottom(const int& p) {
  return std::abs(p) == 5;
}

template <class T>
inline bool isTop(const T& p) {
  return isTop(p->pdg_id());
}
template <>
inline bool isTop(const int& p) {
  return std::abs(p) == 6;
}

/// APID: the fourth generation leptons are leptons.
template <class T>
inline bool isLepton(const T& p) {
  return isLepton(p->pdg_id());
}
template <>
inline bool isLepton(const int& p) {
  auto sp = std::abs(p);
  return sp >= 11 && sp <= 18;
}
template <>
inline bool isLepton(const DecodedPID& p) {
  return isLepton(p.pid());
}

template <class T>
inline bool isSMLepton(const T& p) {
  return isSMLepton(p->pdg_id());
}
template <>
inline bool isSMLepton(const int& p) {
  auto sp = std::abs(p);
  return sp >= 11 && sp <= 16;
}
template <>
inline bool isSMLepton(const DecodedPID& p) {
  return isSMLepton(p.pid());
}

/// APID: the fourth generation leptons are leptons.
template <class T>
inline bool isChLepton(const T& p) {
  return isChLepton(p->pdg_id());
}
template <>
inline bool isChLepton(const int& p) {
  auto sp = std::abs(p);
  return sp >= 11 && sp <= 18 && sp % 2 == 1;
}

template <class T>
inline bool isElectron(const T& p) {
  return isElectron(p->pdg_id());
}
template <>
inline bool isElectron(const int& p) {
  return std::abs(p) == ELECTRON;
}

template <class T>
inline bool isMuon(const T& p) {
  return isMuon(p->pdg_id());
}
template <>
inline bool isMuon(const int& p) {
  return std::abs(p) == MUON;
}

template <class T>
inline bool isTau(const T& p) {
  return isTau(p->pdg_id());
}
template <>
inline bool isTau(const int& p) {
  return std::abs(p) == TAU;
}

/// APID: the fourth generation neutrinos are neutrinos.
template <class T>
inline bool isNeutrino(const T& p) {
  return isNeutrino(p->pdg_id());
}
template <>
inline bool isNeutrino(const int& p) {
  auto sp = std::abs(p);
  return sp == NU_E || sp == NU_MU || sp == NU_TAU || sp == 18;
}

template <class T>
inline bool isSMNeutrino(const T& p) {
  return isSMNeutrino(p->pdg_id());
}
template <>
inline bool isSMNeutrino(const int& p) {
  auto sp = std::abs(p);
  return sp == NU_E || sp == NU_MU || sp == NU_TAU;
}

/// PDG rule 4
/// Diquarks have 4-digit numbers with nq1 >= nq2 and nq3 = 0
/// APID: the diquarks with fourth generation are not diquarks
template <class T>
inline bool isDiquark(const T& p) {
  return isDiquark(p->pdg_id());
}
template <>
inline bool isDiquark(const DecodedPID& p) {
  if (p.ndigits() == 4 && p(0) >= p(1) && p(2) == 0 && p.last() % 2 == 1 &&
      p.max_digit(2, 4) <= TQUARK)
    return true;
  return false;
}
template <>
inline bool isDiquark(const int& p) {
  auto value_digits = DecodedPID(p);
  return isDiquark(value_digits);
}

/// Table 43.1
///  PDG rule 5a:
///  The numbers specifying the meson’s quark content conform to the convention
///  nq1= 0 and nq2 >= nq3. The special case K0L is the sole exception to this
///  rule. PDG rule 5C: The special numbers 310 and 130 are given to the K0S and
///  K0L respectively. APID: The special code K0 is used when a generator uses
///  K0S/K0L
template <class T>
inline bool isMeson(const T& p) {
  return isMeson(p->pdg_id());
}
template <>
inline bool isMeson(const DecodedPID& p) {
  if (p.ndigits() < 3)
    return false;
  if (p.ndigits() == 7 && (p(0) == 1 || p(0) == 2))
    return false;  // APID don't match SUSY particles
  if (std::abs(p.pid()) == K0S)
    return true;
  if (std::abs(p.pid()) == K0L)
    return true;
  if (std::abs(p.pid()) == K0)
    return true;
  if (p.last() % 2 != 1)
    return false;
  if (p.max_digit(1, 3) >= 6)
    return false;
  if (p.max_digit(1, 3) == 0)
    return false;
  if (p.ndigits() > 3 && *(p.second.rbegin() + 3) != 0)
    return false;

  if (p.ndigits() == 3 && p(0) == p(1) && p.pid() < 0)
    return false;
  if (p.ndigits() == 5 && p(2) == p(3) && p.pid() < 0)
    return false;
  if (p.ndigits() == 7 && p(4) == p(5) && p.pid() < 0)
    return false;

  if (p.ndigits() == 3 && p(0) >= p(1) && p(1) != 0)
    return true;
  if (p.ndigits() == 5 && p(2) >= p(3) && p(3) != 0 && p(0) == 1 && p(1) == 0)
    return true;
  if (p.ndigits() == 5 && p(2) >= p(3) && p(3) != 0 && p(0) == 2 && p(1) == 0 &&
      p.last() > 1)
    return true;
  if (p.ndigits() == 5 && p(2) >= p(3) && p(3) != 0 && p(0) == 3 && p(1) == 0 &&
      p.last() > 1)
    return true;

  if (p.ndigits() == 6 && p(3) >= p(4) && p(4) != 0 && p.last() % 2 == 1)
    return true;

  if (p.ndigits() == 7 && p(0) == 9 && p(1) == 0 && p(4) >= p(5) && p(5) != 0)
    return true;

  return false;
}
template <>
inline bool isMeson(const int& p) {
  auto value_digits = DecodedPID(p);
  return isMeson(value_digits);
}

/// Table 43.2
template <class T>
inline bool isBaryon(const T& p) {
  return isBaryon(p->pdg_id());
}
template <>
inline bool isBaryon(const DecodedPID& p) {
  if (p.ndigits() < 4)
    return false;
  if (p.max_digit(1, 4) >= 6)
    return false;
  if (p.min_digit(1, 4) == 0)
    return false;
  if (p.ndigits() == 4 &&
      (p.last() == 2 || p.last() == 4 || p.last() == 6 || p.last() == 8))
    return true;

  if (p.ndigits() == 5 && p(0) == 1 && (p.last() == 2 || p.last() == 4))
    return true;
  if (p.ndigits() == 5 && p(0) == 3 && (p.last() == 2 || p.last() == 4))
    return true;

  if (p.ndigits() == 6) {
    if (p(0) == 1 && p(1) == 0 && p.last() == 2)
      return true;
    if (p(0) == 1 && p(1) == 1 && p.last() == 2)
      return true;
    if (p(0) == 1 && p(1) == 2 && p.last() == 4)
      return true;

    if (p(0) == 2 && p(1) == 0 && p.last() == 2)
      return true;
    if (p(0) == 2 && p(1) == 0 && p.last() == 4)
      return true;
    if (p(0) == 2 && p(1) == 1 && p.last() == 2)
      return true;

    if (p(0) == 1 && p(1) == 0 && p.last() == 4)
      return true;
    if (p(0) == 1 && p(1) == 0 && p.last() == 6)
      return true;
    if (p(0) == 2 && p(1) == 0 && p.last() == 6)
      return true;
    if (p(0) == 2 && p(1) == 0 && p.last() == 8)
      return true;
  }

  if (p.ndigits() == 5) {
    if (p(0) == 2 && p.last() == 2)
      return true;
    if (p(0) == 2 && p.last() == 4)
      return true;
    if (p(0) == 2 && p.last() == 6)
      return true;
    if (p(0) == 5 && p.last() == 2)
      return true;
    if (p(0) == 1 && p.last() == 6)
      return true;
    if (p(0) == 4 && p.last() == 2)
      return true;
  }
  return false;
}
template <>
inline bool isBaryon(const int& p) {
  auto value_digits = DecodedPID(p);
  return isBaryon(value_digits);
}

/// PDG rule 14
/// The 9-digit tetra-quark codes are±1nrnLnq1nq20nq3nq4nJ. For the
/// particleq1q2is a diquarkand
/// ̄q3 ̄q4an antidiquark, sorted such thatnq1≥nq2,nq3≥nq4,nq1≥nq3,
/// andnq2≥nq4ifnq1=nq3.
/// For the antiparticle, given with a negative sign, ̄q1 ̄q2is an antidiquark
/// andq3q4a diquark,
/// with the same sorting except that eithernq1> nq3ornq2> nq4(so
/// thatflavour-diagonal states are particles). Thenr,nL, andnJnumbers have the
/// same meaningas for ordinary hadrons.
template <class T>
inline bool isTetraquark(const T& p) {
  return isTetraquark(p->pdg_id());
}
template <>
inline bool isTetraquark(const DecodedPID& p) {
  return (p.ndigits() == 9 && p(0) == 1 && p(5) == 0 &&
          p.max_digit(1, 3) <= 6 && p.min_digit(1, 3) > 0 &&
          p.max_digit(1 + 3, 3 + 3) <= 6 && p.min_digit(1 + 3, 3 + 3) > 0 &&
          (p(3) >= p(4) && p(6) >= p(7)) &&
          ((p(3) > p(6)) || (p(3) == p(6) && (p(4) >= p(7)))));
}
template <>
inline bool isTetraquark(const int& p) {
  auto value_digits = DecodedPID(p);
  return isTetraquark(value_digits);
}

/// PDG rule 15
/// The 9-digit penta-quark codes are±1nrnLnq1nq2nq3nq4nq5nJ, sorted such
/// thatnq1≥nq2≥nq3≥nq4. In the particle the first four are quarks and the fifth
/// an antiquark while t
/// heopposite holds in the antiparticle, which is given with a negative sign.
/// Thenr,nL, andnJnumbers have the same meaning as for ordinary hadrons.
template <class T>
inline bool isPentaquark(const T& p) {
  return isPentaquark(p->pdg_id());
}
template <>
inline bool isPentaquark(const DecodedPID& p) {
  return (p.ndigits() == 9 && p(0) == 1 && p.max_digit(1, 6) <= 6 &&
          p.min_digit(1, 6) > 0 &&
          (p(3) >= p(4) && p(4) >= p(5) && p(5) >= p(6)));
}
template <>
inline bool isPentaquark(const int& p) {
  auto value_digits = DecodedPID(p);
  return isPentaquark(value_digits);
}

// APID Mesons, Baryons, Tetraquarks and Pentaquarks are Hadrons
template <class T>
inline bool isHadron(const T& p) {
  return isHadron(p->pdg_id());
}
template <>
inline bool isHadron(const DecodedPID& p) {
  return isMeson(p) || isBaryon(p) || isTetraquark(p) || isPentaquark(p);
}
template <>
inline bool isHadron(const int& p) {
  auto value_digits = DecodedPID(p);
  return isHadron(value_digits);
}

/// PDG rule 8:
/// The pomeron and odderon trajectories and a generic reggeon trajectory
/// of states in QCD areassigned codes 990, 9990, and 110 respectively
template <class T>
inline bool isTrajectory(const T& p) {
  return isTrajectory(p->pdg_id());
}
template <>
inline bool isTrajectory(const int& p) {
  return std::abs(p) == POMERON || std::abs(p) == ODDERON ||
         std::abs(p) == REGGEON;
}

/// PDG rule 9:
/// Two-digit numbers in the range 21–30 are provided for the Standard
/// Model gauge and Higgs bosons.
/// PDG rule 11b:
/// The graviton and the boson content of a two-Higgs-doublet scenario
/// and of additional SU(2)×U(1) groups are found in the range 31–40.
template <class T>
inline bool isBoson(const T& p) {
  return isBoson(p->pdg_id());
}
template <>
inline bool isBoson(const int& p) {
  auto sp = std::abs(p);
  return sp > 20 && sp < 41;
}
template <>
inline bool isBoson(const DecodedPID& p) {
  return isBoson(p.pid());
}

template <class T>
inline bool isGluon(const T& p) {
  return isGluon(p->pdg_id());
}
template <>
inline bool isGluon(const int& p) {
  return p == GLUON;
}

template <class T>
inline bool isPhoton(const T& p) {
  return isPhoton(p->pdg_id());
}
template <>
inline bool isPhoton(const int& p) {
  return p == PHOTON;
}

template <class T>
inline bool isZ(const T& p) {
  return isZ(p->pdg_id());
}
template <>
inline bool isZ(const int& p) {
  return p == Z0BOSON;
}

template <class T>
inline bool isW(const T& p) {
  return isW(p->pdg_id());
}
template <>
inline bool isW(const int& p) {
  return std::abs(p) == WPLUSBOSON;
}

/// APID: Additional "Heavy"/"prime" versions of W and Z bosons (Used in
/// MCTruthClassifier)
template <class T>
inline bool isHeavyBoson(const T& p) {
  return isHeavyBoson(p->pdg_id());
}
template <>
inline bool isHeavyBoson(const int& p) {
  return p == ZPRIME || p == ZDBLPRIME || std::abs(p) == WPLUSPRIME;
}

/// APID: HIGGS boson is only one particle.
template <class T>
inline bool isHiggs(const T& p) {
  return isHiggs(p->pdg_id());
}
template <>
inline bool isHiggs(const int& p) {
  return p == HIGGSBOSON;
}

/// APID: Additional Higgs bosons for MSSM (Used in MCTruthClassifier)
template <class T>
inline bool isMSSMHiggs(const T& p) {
  return isMSSMHiggs(p->pdg_id());
}
template <>
inline bool isMSSMHiggs(const int& p) {
  return p == HIGGS2 || p == HIGGS3 || std::abs(p) == HIGGSPLUS;
}

template <class T>
inline bool isGraviton(const T& p) {
  return isGraviton(p->pdg_id());
}
template <>
inline bool isGraviton(const int& p) {
  return p == GRAVITON;
}

template <class T>
inline bool isResonance(const T& p) {
  return isZ(p) || isW(p) || isHiggs(p) || isTop(p);
}  // APID: not including t' (pdg_id=8), Z', Z'' and W'+ or BSM Higgs bosons

/// PDG rule 11c:
/// “One-of-a-kind” exotic particles are assigned numbers in the range
/// 41–80. The subrange 61-80 can be used for new heavier fermions in
/// generic models, where partners to the SM fermions would have codes
/// oﬀset by 60. If required, however, other assignments could be
/// made.
template <class T>
inline bool isLeptoQuark(const T& p) {
  return isLeptoQuark(p->pdg_id());
}
template <>
inline bool isLeptoQuark(const int& p) {
  return std::abs(p) == LEPTOQUARK;
}

template <class T>
inline bool isPythia8Specific(const T& p) {
  return isPythia8Specific(p->pdg_id());
}
template <>
inline bool isPythia8Specific(const DecodedPID& p) {
  return (p.ndigits() == 7 && p(0) == 9 && p(1) == 9);
}
template <>
inline bool isPythia8Specific(const int& p) {
  auto value_digits = DecodedPID(p);
  return isPythia8Specific(value_digits);
}

/// PDG Rule 12:
/// APID: Helper function for right-handed neutrino states
/// These are generator defined PDG ID values for right handed
/// neutrinos. (Defined for some MadGraph+Pythia8 samples and
/// referenced in MCTruthClassifierGen.cxx)
template <class T>
inline bool isNeutrinoRH(const T& p) {
  return isNeutrinoRH(p->pdg_id());
}
template <>
inline bool isNeutrinoRH(const int& p) {
  return (std::abs(p) == RH_NU_E || std::abs(p) == RH_NU_MU ||
          std::abs(p) == RH_NU_TAU);
}

/// Main Table
/// for MC internal use 81–100,901–930,998-999,1901–1930,2901–2930, and
/// 3901–3930
template <class T>
inline bool isGenSpecific(const T& p) {
  return isGenSpecific(p->pdg_id());
}
template <>
inline bool isGenSpecific(const int& p) {
  if (p >= 81 && p <= 100)
    return true;
  if (p >= 901 && p <= 930)
    return true;
  if (p >= 998 && p <= 999)
    return true;
  if (p >= 1901 && p <= 1930)
    return true;
  if (p >= 2901 && p <= 2930)
    return true;
  if (p >= 3901 && p <= 3930)
    return true;
  return false;
}

template <class T>
inline bool isGeantino(const T& p) {
  return isGeantino(p->pdg_id());
}
template <>
inline bool isGeantino(const int& p) {
  return (std::abs(p) == GEANTINO0 || std::abs(p) == GEANTINOPLUS);
}

/// APID: Definition of Glueballs: SM glueballs 99X (X=1,5), 999Y (Y=3,7)
template <class T>
inline bool isGlueball(const T& p) {
  return isGlueball(p->pdg_id());
}
template <>
inline bool isGlueball(const DecodedPID& p) {
  if (p.ndigits() > 4)
    return false;  // APID avoid classifying R-Glueballs as SM Glueballs
  return ((p.ndigits() == 3 && p(0) == COMPOSITEGLUON &&
           p(1) == COMPOSITEGLUON && (p(2) == 1 || p(2) == 5)) ||
          (p.ndigits() == 4 && p(0) == COMPOSITEGLUON &&
           p(1) == COMPOSITEGLUON && p(2) == COMPOSITEGLUON &&
           (p(3) == 3 || p(3) == 7)));
}
template <>
inline bool isGlueball(const int& p) {
  auto value_digits = DecodedPID(p);
  return isGlueball(value_digits);
}

/// PDG rule 11d
/// Fundamental supersymmetric particles are identified by adding a nonzero n to
/// the particle number. The superpartner of a boson or a left-handed fermion
/// has n = 1 while the superpartner of a right-handed fermion has n = 2. When
/// mixing occurs, such as between the winos and charged Higgsinos to give
/// charginos, or between left and right sfermions, the lighter physical state
/// is given the smaller basis state number.
template <class T>
inline bool isSUSY(const T& p) {
  return isSUSY(p->pdg_id());
}
template <>
inline bool isSUSY(const DecodedPID& p) {
  return (p.ndigits() == 7 && (p(0) == 1 || p(0) == 2) &&
          !isGenSpecific(p.shift(2).pid()));
}
template <>
inline bool isSUSY(const int& p) {
  auto value_digits = DecodedPID(p);
  return isSUSY(value_digits);
}

// APID: Super-partners of standard model quarks only
template <class T>
inline bool isSquark(const T& p) {
  return isSquark(p->pdg_id());
}
template <>
inline bool isSquark(const DecodedPID& p) {
  auto pp = p.shift(1);
  return isSUSY(p) && isSMQuark(pp);
}
template <>
inline bool isSquark(const int& p) {
  auto value_digits = DecodedPID(p);
  return isSquark(value_digits);
}

// APID: Super-partners of left-handed standard model quarks only
template <class T>
inline bool isSquarkLH(const T& p) {
  return isSquarkLH(p->pdg_id());
}
template <>
inline bool isSquarkLH(const DecodedPID& p) {
  return isSquark(p) && (p(0) == 1);
}
template <>
inline bool isSquarkLH(const int& p) {
  auto value_digits = DecodedPID(p);
  return isSquarkLH(value_digits);
}

// APID: Super-partners of right-handed standard model quarks only
template <class T>
inline bool isSquarkRH(const T& p) {
  return isSquarkRH(p->pdg_id());
}
template <>
inline bool isSquarkRH(const DecodedPID& p) {
  return isSquark(p) && (p(0) == 2);
}
template <>
inline bool isSquarkRH(const int& p) {
  auto value_digits = DecodedPID(p);
  return isSquarkRH(value_digits);
}

template <class T>
inline bool hasSquark(const T& p, const int& q) {
  return hasSquark(p->pdg_id(), q);
}
template <>
inline bool hasSquark(const DecodedPID& p, const int& q) {
  auto pp = p.shift(1);
  return isSUSY(p) && pp.ndigits() != 2 &&
         pp(0) ==
             q;  // skip lepton and boson super-partners by vetoing ndigits==2
}
template <>
inline bool hasSquark(const int& p, const int& q) {
  auto value_digits = DecodedPID(p);
  return hasSquark(value_digits, q);
}

// APID: Super-partners of standard model leptons only
template <class T>
inline bool isSlepton(const T& p) {
  return isSlepton(p->pdg_id());
}
template <>
inline bool isSlepton(const DecodedPID& p) {
  auto pp = p.shift(1);
  return isSUSY(p) && isSMLepton(pp);
}
template <>
inline bool isSlepton(const int& p) {
  auto value_digits = DecodedPID(p);
  return isSlepton(value_digits);
}

// APID: Super-partners of left-handed standard model leptons only
template <class T>
inline bool isSleptonLH(const T& p) {
  return isSleptonLH(p->pdg_id());
}
template <>
inline bool isSleptonLH(const DecodedPID& p) {
  return isSlepton(p) && (p(0) == 1);
}
template <>
inline bool isSleptonLH(const int& p) {
  auto value_digits = DecodedPID(p);
  return isSleptonLH(value_digits);
}

// APID: Super-partners of right-handed standard model leptons only
template <class T>
inline bool isSleptonRH(const T& p) {
  return isSleptonRH(p->pdg_id());
}
template <>
inline bool isSleptonRH(const DecodedPID& p) {
  return isSlepton(p) && (p(0) == 2);
}
template <>
inline bool isSleptonRH(const int& p) {
  auto value_digits = DecodedPID(p);
  return isSleptonRH(value_digits);
}

// APID: Super-partners of gauge bosons including gravitons
template <class T>
inline bool isGaugino(const T& p) {
  return isGaugino(p->pdg_id());
}
template <>
inline bool isGaugino(const DecodedPID& p) {
  auto pp = p.shift(1);
  return isSUSY(p) && isBoson(pp.pid());
}
template <>
inline bool isGaugino(const int& p) {
  auto value_digits = DecodedPID(p);
  return isGaugino(value_digits);
}

/// PDG rule 11e
/// Technicolor states have n = 3, with technifermions treated like ordinary
/// fermions. States which are ordinary color singlets have n_r = 0. Color
/// octets have n_r = 1. If a state has non-trivial quantum numbers under the
/// topcolor groups SU(3)1×SU(3)2, the quantum numbers are specified by tech,
/// ij, where i and j are 1 or 2. nLis then 2i+j. The coloron V8, is a heavy
/// gluon color octet and thus is 3100021
template <class T>
inline bool isTechnicolor(const T& p) {
  return isTechnicolor(p->pdg_id());
}
template <>
inline bool isTechnicolor(const DecodedPID& p) {
  const auto& pp = (p.ndigits() == 7) ? p.shift(2) : DecodedPID(0);
  return (p.ndigits() == 7 && p(0) == 3 && (p(1) == 0 || p(0) == 1) &&
          (isQuark(pp) || isLepton(pp) || isBoson(pp) || isGlueball(pp) ||
           isDiquark(pp) || isHadron(pp)));
}
template <>
inline bool isTechnicolor(const int& p) {
  auto value_digits = DecodedPID(p);
  return isTechnicolor(value_digits);
}

/// PDG rule 11f
/// Excited (composite) quarks and leptons are identified by setting n= 4 and
/// nr= 0
template <class T>
inline bool isExcited(const T& p) {
  return isExcited(p->pdg_id());
}
template <>
inline bool isExcited(const DecodedPID& p) {
  const auto& pp = (p.ndigits() == 7) ? p.shift(2) : DecodedPID(0);
  return (p.ndigits() == 7 && (p(0) == 4 && p(1) == 0) &&
          (isLepton(pp) || isQuark(pp)));
}
template <>
inline bool isExcited(const int& p) {
  auto value_digits = DecodedPID(p);
  return isExcited(value_digits);
}

/// PDG rule 11g:
/// Within several scenarios of new physics, it is possible to have colored
/// particles suﬃciently long-lived for color-singlet hadronic states to form
/// around them. In the context of supersymmetric scenarios, these states are
/// called R-hadrons, since they carry odd R- parity. R-hadron codes, deﬁned
/// here, should be viewed as templates for corresponding codes also in other
/// scenarios, for any long-lived particle that is either an unﬂavored color
/// octet or a ﬂavored color triplet. The R-hadron code is obtained by combining
/// the SUSY particle code with a code for the light degrees of freedom, with as
/// many intermediate zeros removed from the former as required to make place
/// for the latter at the end. (To exemplify, a sparticle n00000n˜q combined
/// with quarks q1 and q2 obtains code n00n˜qnq1 nq2 nJ .) Speciﬁcally, the
/// new-particle spin decouples in the limit of large masses, so that the ﬁnal
/// nJ digit is deﬁned by the spin state of the light-quark system alone. An
/// appropriate number of nq digits is used to deﬁne the ordinary-quark content.
/// As usual, 9 rather than 21 is used to denote a gluon/gluino in composite
/// states. The sign of the hadron agrees with that of the constituent new
/// particle (a color triplet) where there is a distinct new antiparticle, and
/// else is deﬁned as for normal hadrons. Particle names are R with the ﬂavor
/// content as lower index.

/// APID: Definition of R-Glueballs: 100099X (X=1,3), 100999Y (Y=1,5)
/// APID: NB In the current numbering scheme, some states with 2
/// gluinos + gluon or 2 gluons + gluino could have degenerate
/// PDG_IDs.
template <class T>
inline bool isRGlueball(const T& p) {
  return isRGlueball(p->pdg_id());
}
template <>
inline bool isRGlueball(const DecodedPID& p) {
  if (p.ndigits() != 7 || p(0) != 1)
    return false;
  auto pp = p.shift(1);
  return ((pp.ndigits() == 3 && pp(0) == COMPOSITEGLUON &&
           pp(1) == COMPOSITEGLUON && (pp(2) == 1 || pp(2) == 3)) ||
          (pp.ndigits() == 4 && pp(0) == COMPOSITEGLUON &&
           pp(1) == COMPOSITEGLUON && pp(2) == COMPOSITEGLUON &&
           (pp(3) == 1 || pp(3) == 5)));
}
template <>
inline bool isRGlueball(const int& p) {
  auto value_digits = DecodedPID(p);
  return isRGlueball(value_digits);
}

// APID Define R-Mesons as gluino-quark-antiquark and squark-antiquark bound
// states (ignore 4th generation squarks/quarks) NB Current models only allow
// gluino-quark-antiquark, stop-antiquark and sbottom-antiquark states
template <class T>
inline bool isRMeson(const T& p) {
  return isRMeson(p->pdg_id());
}
template <>
inline bool isRMeson(const DecodedPID& p) {
  if (!(p.ndigits() == 7 && (p(0) == 1 || p(0) == 2)))
    return false;
  auto pp = p.shift(1);
  return (
      // Handle ~gluino-quark-antiquark states
      (pp.ndigits() == 4 && pp(0) == COMPOSITEGLUON &&
       pp.max_digit(1, 3) < COMPOSITEGLUON && pp(2) <= pp(1) &&
       isSMQuark(pp(1)) && isSMQuark(pp(2)) &&
       (pp.last() == 1 || pp.last() == 3)) ||
      // Handle squark-antiquark states (previously called Smeson/mesoninos)
      (pp.ndigits() == 3 && pp.max_digit(1, 3) < COMPOSITEGLUON &&
       pp(1) <= pp(0) && isSMQuark(pp(0)) && isSMQuark(pp(1)) &&
       pp.last() == 2));
}
template <>
inline bool isRMeson(const int& p) {
  auto value_digits = DecodedPID(p);
  return isRMeson(value_digits);
}

// APID Define R-Baryons as gluino-quark-quark-quark and squark-quark-quark
// bound states (ignore 4th generation squarks/quarks) NB Current models only
// allow gluino-quark-quark-quark, stop-quark-quark and sbottom-quark-quark
// states
template <class T>
inline bool isRBaryon(const T& p) {
  return isRBaryon(p->pdg_id());
}
template <>
inline bool isRBaryon(const DecodedPID& p) {
  if (!(p.ndigits() == 7 && (p(0) == 1 || p(0) == 2)))
    return false;
  auto pp = p.shift(1);
  return (
      // Handle ~gluino-quark-quark-quark states
      (pp.ndigits() == 5 && pp(0) == COMPOSITEGLUON &&
       pp.max_digit(1, 4) < COMPOSITEGLUON && pp(2) <= pp(1) &&
       pp(3) <= pp(2) && isSMQuark(pp(1)) && isSMQuark(pp(2)) &&
       isSMQuark(pp(3)) && (pp.last() == 2 || pp.last() == 4)) ||
      // Handle squark-quark-quark states (previously called Sbaryons)
      (pp.ndigits() == 4 && pp.max_digit(1, 4) < COMPOSITEGLUON &&
       pp(1) <= pp(0) && pp(2) <= pp(1) && isSMQuark(pp(0)) &&
       isSMQuark(pp(1)) && isSMQuark(pp(2)) &&
       (pp.last() == 1 || pp.last() == 3)));
}
template <>
inline bool isRBaryon(const int& p) {
  auto value_digits = DecodedPID(p);
  return isRBaryon(value_digits);
}

/// PDG rule 11h
/// A black hole in models with extra dimensions has code 5000040. Kaluza-Klein
/// excitations in models with extra dimensions have n = 5 or n = 6, to
/// distinguish excitations of left-or right-handed fermions or, in case of
/// mixing, the lighter or heavier state (cf. 11d). The non zero nr digit gives
/// the radial excitation number, in scenarios where the level spacing allows
/// these to be
///  distinguished. Should the model also contain supersymmetry, excited SUSY
///  states would be denoted by a nn_r > 0, with n = 1 or 2 as usual.
/// Should some colored states be long-lived enough that hadrons would form
/// around them, the coding strategy of 11g applies, with the initial two nnr
/// digits preserved in the combined code.
template <class T>
inline bool isKK(const T& p) {
  return isKK(p->pdg_id());
}
template <>
inline bool isKK(const DecodedPID& p) {
  return (p.ndigits() == 7 && (p(0) == 5 || p(0) == 6));
}
template <>
inline bool isKK(const int& p) {
  auto value_digits = DecodedPID(p);
  return isKK(value_digits);
}

/// PDG rule 11i
/// Magnetic monopoles and dyons are assumed to have one unit of Dirac monopole
/// charge and a variable integer number nq1nq2 nq3 units of electric charge.
/// Codes 411nq1nq2 nq3 0 are then used when the magnetic and electrical charge
/// sign agree and 412nq1nq2 nq3 0 when they disagree, with the overall sign of
/// the particle set by the magnetic charge. For now no spin information is
/// provided.
template <class T>
inline bool isMonopole(const T& p) {
  return isMonopole(p->pdg_id());
}
template <>
inline bool isMonopole(const DecodedPID& p) {
  return (p.ndigits() == 7 && p(0) == 4 && p(1) == 1 &&
          (p(2) == 1 || p(2) == 2) && p(6) == 0);
}
template <>
inline bool isMonopole(const int& p) {
  auto value_digits = DecodedPID(p);
  return isMonopole(value_digits);
}

/// PDG rule 11j:
/// The nature of Dark Matter (DM) is not known, and therefore a definitive
/// classificationis too early. Candidates within specific scenarios are
/// classified therein, such as 1000022 for the lightest neutralino.
/// Generic fundamental states can be given temporary codes in the range 51 -
/// 60, with 51, 52 and 53 reserved for spin 0, 1/2 and 1 ones (this could also
/// be an axion state). Generic mediators of s-channel DM pair creation of
/// annihilation can be given codes 54 and 55 for spin 0 or 1 ones. Separate
/// antiparticles, with negativecodes, may or may not exist. More elaborate new
/// scenarios should be constructed with n= 5 and nr = 9. APID: Only the 51-60
/// range is considered DM. The antiparticles are assumed to exist.
template <class T>
inline bool isDM(const T& p) {
  return isDM(p->pdg_id());
}
template <>
inline bool isDM(const int& p) {
  auto sp = std::abs(p);
  return (sp >= 51 && sp <= 60) || sp == DARKPHOTON;
}

/// PDG rule 11k
/// Hidden Valley particles have n = 4 and n_r = 9, and trailing numbers in
/// agreement with their nearest-analog standard particles, as far as possible.
/// Thus 4900021 is the gauge boson g_v of a confining gauge field,
/// 490000n_{q_v} and 490001n_{l_v} fundamental constituents charged or not
/// under this, 4900022 is the γ_v of a non-confining field, and
/// 4900n_{q_{v1}}n_{q_{v2}}n_J a Hidden Valley meson.
template <class T>
inline bool isHiddenValley(const T& p) {
  return isHiddenValley(p->pdg_id());
}
template <>
inline bool isHiddenValley(const DecodedPID& p) {
  const auto& pp = (p.ndigits() == 7) ? p.shift(2) : DecodedPID(0);
  return (p.ndigits() == 7 && p(0) == 4 && p(1) == 9 &&
          (isQuark(pp) || isLepton(pp) || isBoson(pp) || isGlueball(pp) ||
           isDiquark(pp) || isHadron(pp)));
}
template <>
inline bool isHiddenValley(const int& p) {
  auto value_digits = DecodedPID(p);
  return isHiddenValley(value_digits);
}

/// In addition, there is a need to identify ”Q-ball” and similar very exotic
/// (multi-charged) particles which may have large, non-integer charge. These
/// particles are assigned the ad-hoc numbering +/-100XXXY0, where the charge is
/// XXX.Y. or +/-200XXYY0, where the charge is XX/YY. The case of +/-200XXYY0 is
/// legacy, see https://gitlab.cern.ch/atlas/athena/-/merge_requests/25862 Note
/// that no other quantum numbers besides the charge are considered for these
/// generic multi-charged particles (e.g. isSUSY() is false for them). Such a
/// model was used in previous Run-1 (1301.5272,1504.04188) and Run-2
/// (1812.03673,2303.13613) ATLAS searches.
template <class T>
inline bool isGenericMultichargedParticle(const T& p) {
  return isGenericMultichargedParticle(p->pdg_id());
}
template <>
inline bool isGenericMultichargedParticle(const DecodedPID& p) {
  return (p.ndigits() == 8 && (p(0) == 1 || p(0) == 2) && p(1) == 0 &&
          p(2) == 0 && p(7) == 0);
}
template <>
inline bool isGenericMultichargedParticle(const int& p) {
  auto value_digits = DecodedPID(p);
  return isGenericMultichargedParticle(value_digits);
}

/// PDG rule 16
/// Nuclear codes are given as 10-digit numbers ±10LZZZAAAI.
/// For a (hyper)nucleus consisting of n_p protons, n_n neutrons and
/// n_Λ Λ’s:
/// A = n_p + n_n + n_Λ gives the total baryon number,
/// Z = n_p gives the total charge,
/// L = n_Λ gives the total number of strange quarks.
/// I gives the isomer level, with I= 0 corresponding to the ground
/// state and I > 0 to excitations, see
/// [http://www.nndc.bnl.gov/amdc/web/nubase en.html], where states
/// denoted m, n, p ,q translate to I= 1–4. As examples, the deuteron
/// is 1000010020 and 235U is 1000922350. To avoid ambiguities,
/// nuclear codes should not be applied to a single hadron, like p, n or
/// Λ^0, where quark-contents-based codes already exist.
template <class T>
inline bool isNucleus(const T& p) {
  return isNucleus(p->pdg_id());
}
template <>
inline bool isNucleus(const DecodedPID& p) {
  if (std::abs(p.pid()) == PROTON)
    return true;
  return (p.ndigits() == 10 && p(0) == 1 && p(1) == 0);
}
template <>
inline bool isNucleus(const int& p) {
  auto value_digits = DecodedPID(p);
  return isNucleus(value_digits);
}

template <class T>
inline bool hasQuark(const T& p, const int& q);
template <>
inline bool hasQuark(const DecodedPID& p, const int& q) {
  if (isQuark(p.pid())) {
    return (std::abs(p.pid()) == q);
  }
  if (isMeson(p)) {
    return *(p.second.rbegin() + 1) == q || *(p.second.rbegin() + 2) == q;
  }
  if (isDiquark(p)) {
    auto i = std::find(p.second.rbegin() + 2, p.second.rbegin() + 4, q);
    return (i != p.second.rbegin() + 4);
  }
  if (isBaryon(p)) {
    auto i = std::find(p.second.rbegin() + 1, p.second.rbegin() + 4, q);
    return (i != p.second.rbegin() + 4);
  }
  if (isTetraquark(p)) {
    auto i = std::find(p.second.rbegin() + 1, p.second.rbegin() + 5, q);
    return (i != p.second.rbegin() + 5);
  }
  if (isPentaquark(p)) {
    auto i = std::find(p.second.rbegin() + 1, p.second.rbegin() + 6, q);
    return (i != p.second.rbegin() + 6);
  }
  if (isNucleus(p) && std::abs(p.pid()) != PROTON) {
    return (q == 1 || q == 2 || (q == 3 && p(2) > 0));
  }
  if (isSUSY(p)) {  // APID SUSY case
    auto pp = p.shift(1);
    if (pp.ndigits() == 1) {
      return false;
    }  // Handle squarks
    if (pp.ndigits() == 3) {
      return (pp(1) == q);
    }  // Handle ~q qbar pairs
    if (pp.ndigits() == 4) {
      return (pp(1) == q || pp(2) == q);
    }  // Ignore gluinos and squarks
    if (pp.ndigits() == 5) {
      return (pp(1) == q || pp(2) == q || pp(3) == q);
    }  // Ignore gluinos and squarks
    if (pp.ndigits() > 5) {
      pp = pp.shift(1);
    }  // Drop gluinos and squarks
    return hasQuark(pp, q);
  }
  return false;
}
template <>
inline bool hasQuark(const int& p, const int& q) {
  auto value_digits = DecodedPID(p);
  return hasQuark(value_digits, q);
}

template <class T>
inline bool hasStrange(const T& p) {
  return hasQuark(p, SQUARK);
}
template <class T>
inline bool hasCharm(const T& p) {
  return hasQuark(p, CQUARK);
}
template <class T>
inline bool hasBottom(const T& p) {
  return hasQuark(p, BQUARK);
}
template <class T>
inline bool hasTop(const T& p) {
  return hasQuark(p, TQUARK);
}

// APID: The baryon number is defined as:
// B = (1/3)*( n_q - n_{qbar} )
// where n_q⁠ is the number of quarks, and ⁠n_{qbar} is the number of
// antiquarks. By convention, squarks have the same quantum numbers as
// the corresponding quarks (modulo spin and R), so have baryon number
// 1/3.
template <class T>
inline int baryonNumber3(const T& p) {
  return baryonNumber3(p->pdg_id());
}
template <>
inline int baryonNumber3(const DecodedPID& p) {
  if (isQuark(p.pid())) {
    return (p.pid() > 0) ? 1 : -1;
  }
  if (isDiquark(p)) {
    return (p.pid() > 0) ? 2 : -2;
  }
  if (isMeson(p) || isTetraquark(p)) {
    return 0;
  }
  if (isBaryon(p) || isPentaquark(p)) {
    return (p.pid() > 0) ? 3 : -3;
  }
  if (isNucleus(p)) {
    const int result = 3 * p(8) + 30 * p(7) + 300 * p(6);
    return (p.pid() > 0) ? result : -result;
  }
  if (isSUSY(p)) {
    auto pp = p.shift(1);
    if (pp.ndigits() < 3) {
      return baryonNumber3(pp);
    }  // super-partners of fundamental particles
    if (pp(0) == COMPOSITEGLUON) {
      if (pp(1) == COMPOSITEGLUON) {
        return 0;
      }  // R-Glueballs
      if (pp.ndigits() == 4) {
        return 0;
      }  // states with gluino-quark-antiquark
      if (pp.ndigits() == 5) {
        return (p.pid() > 0) ? 3 : -3;
      }  // states with gluino-quark-quark-quark
    }
    if (pp.ndigits() == 3) {
      return 0;
    }  // squark-antiquark
    if (pp.ndigits() == 4) {
      return (p.pid() > 0) ? 3 : -3;
    }  // states with squark-quark-quark
  }
  return 0;
}
template <>
inline int baryonNumber3(const int& p) {
  auto value_digits = DecodedPID(p);
  return baryonNumber3(value_digits);
}

template <class T>
inline double baryonNumber(const T& p) {
  return baryonNumber(p->pdg_id());
}
template <>
inline double baryonNumber(const DecodedPID& p) {
  return static_cast<double>(baryonNumber3(p)) / 3.0;
}
template <>
inline double baryonNumber(const int& p) {
  auto value_digits = DecodedPID(p);
  return static_cast<double>(baryonNumber3(value_digits)) / 3.0;
}

// APID: The strangeness of a particle is defined as:
// S = − ( n_s − n_{sbar} )
// where n_s represents the number of strange quarks and n_{sbar}
// represents the number of strange antiquarks. By convention, strange
// squarks have the same quantum numbers as strange quarks (modulo
// spin and R), so have strangeness -1.
static const std::array<int, 10> is_strange = {+0, +0, +0, -1, +0,
                                               +0, +0, +0, +0, +0};
template <class T>
inline int strangeness(const T& p) {
  return strangeness(p->pdg_id());
}
template <>
inline int strangeness(const DecodedPID& p) {
  if (isNucleus(p) && p.ndigits() == 10) {
    return (p.pid() > 0) ? -p(2) : p(2);
  }
  if (isStrange(p.pid())) {
    return (p.pid() > 0) ? -1 : 1;
  }
  if (!hasStrange(p) && !hasSquark(p, SQUARK)) {
    return 0;
  }
  if (std::abs(p.pid()) == K0) {
    return (p.pid() > 0) ? 1 : -1;
  }
  std::size_t nq = 0;
  int sign = 1;
  int signmult = 1;
  int result = 0;
  bool classified = false;
  if (!classified && isMeson(p)) {
    classified = true;
    nq = 2;
    if ((*(p.second.rbegin() + 2)) == 2 || (*(p.second.rbegin() + 2)) == 4) {
      sign = -1;
    }
    signmult = -1;
  }
  if (!classified && isDiquark(p)) {
    return is_strange.at(p(0)) + is_strange.at(p(1));
  }
  if (!classified && isBaryon(p)) {
    classified = true;
    nq = 3;
  }
  if (!classified && isTetraquark(p)) {
    return is_strange.at(p(3)) + is_strange.at(p(4)) - is_strange.at(p(6)) -
           is_strange.at(p(7));
  }
  if (!classified && isPentaquark(p)) {
    return is_strange.at(p(3)) + is_strange.at(p(4)) + is_strange.at(p(5)) +
           is_strange.at(p(6)) - is_strange.at(p(7));
  }
  if (!classified && isSUSY(p)) {
    nq = 0;
    auto pp = p.shift(1);
    if (pp.ndigits() < 3) {
      return strangeness(pp);
    }  // super-partners of fundamental particles
    if (pp(0) == COMPOSITEGLUON) {
      if (pp(1) == COMPOSITEGLUON) {
        return 0;
      }  // R-Glueballs
      if (pp.ndigits() == 4 || pp.ndigits() == 5) {
        pp = pp.shift(1);  // Remove gluino
      }
    }
    if (pp.ndigits() == 3) {
      classified = true;
      nq = 2;
      if (p.last() % 2 == 0) {
        sign = -1;
      }
      signmult = -1;
    }  // states with quark-antiquark or squark-antiquark
    if (pp.ndigits() == 4) {
      classified = true;
      nq = 3;
    }  // states with quark-quark-quark or squark-quark-quark
  }
  for (auto r = p.second.rbegin() + 1; r != p.second.rbegin() + 1 + nq; ++r) {
    result += is_strange.at(*r) * sign;
    sign *= signmult;
  }
  return p.pid() > 0 ? result : -result;
}
template <>
inline int strangeness(const int& p) {
  auto value_digits = DecodedPID(p);
  return strangeness(value_digits);
}

template <class T>
inline int numberOfLambdas(const T& p) {
  return numberOfLambdas(p->pdg_id());
}
template <>
inline int numberOfLambdas(const DecodedPID& p) {
  if (std::abs(p.pid()) == LAMBDA0) {
    return (p.pid() > 0) ? 1 : -1;
  }
  if (isNucleus(p) && p.ndigits() == 10) {
    return (p.pid() > 0) ? p(2) : -p(2);
  }
  return 0;
}
template <>
inline int numberOfLambdas(const int& p) {
  auto value_digits = DecodedPID(p);
  return numberOfLambdas(value_digits);
}

template <class T>
inline int numberOfProtons(const T& p) {
  return numberOfProtons(p->pdg_id());
}
template <>
inline int numberOfProtons(const DecodedPID& p) {
  if (std::abs(p.pid()) == PROTON) {
    return (p.pid() > 0) ? 1 : -1;
  }
  if (isNucleus(p)) {
    const int result = p(5) + 10 * p(4) + 100 * p(3);
    return (p.pid() > 0) ? result : -result;
  }
  return 0;
}
template <>
inline int numberOfProtons(const int& p) {
  auto value_digits = DecodedPID(p);
  return numberOfProtons(value_digits);
}

/// APID: graviton and all Higgs extensions are BSM
template <class T>
inline bool isBSM(const T& p) {
  return isBSM(p->pdg_id());
}
template <>
inline bool isBSM(const DecodedPID& p) {
  if (p.pid() == GRAVITON || std::abs(p.pid()) == MAVTOP ||
      p.pid() == DARKPHOTON)
    return true;
  if (std::abs(p.pid()) > 16 && std::abs(p.pid()) < 19)
    return true;
  if (std::abs(p.pid()) > 31 && std::abs(p.pid()) < 39)
    return true;
  if (std::abs(p.pid()) > 39 && std::abs(p.pid()) < 81)
    return true;
  if (std::abs(p.pid()) > 6 && std::abs(p.pid()) < 9)
    return true;
  if (isSUSY(p))
    return true;
  if (isGenericMultichargedParticle(p))
    return true;
  if (isTechnicolor(p))
    return true;
  if (isExcited(p))
    return true;
  if (isKK(p))
    return true;
  if (isHiddenValley(p))
    return true;
  return false;
}
template <>
inline bool isBSM(const int& p) {
  if (p == GRAVITON || std::abs(p) == MAVTOP || p == DARKPHOTON)
    return true;
  if (std::abs(p) > 16 && std::abs(p) < 19)
    return true;
  if (std::abs(p) > 31 && std::abs(p) < 38)
    return true;
  if (std::abs(p) > 39 && std::abs(p) < 81)
    return true;
  if (std::abs(p) > 6 && std::abs(p) < 9)
    return true;
  auto value_digits = DecodedPID(p);
  return isBSM(value_digits);
}

template <class T>
inline bool isTransportable(const T& p) {
  return isTransportable(p->pdg_id());
}
template <>
inline bool isTransportable(const DecodedPID& p) {
  return isPhoton(p.pid()) || isGeantino(p.pid()) || isHadron(p) ||
         isLepton(p.pid()) || p.pid() == DARKPHOTON;
}
template <>
inline bool isTransportable(const int& p) {
  auto value_digits = DecodedPID(p);
  return isTransportable(value_digits);
}

/// Av: we implement here an ATLAS-sepcific convention: all particles which are
/// 99xxxxx are fine.
template <class T>
inline bool isValid(const T& p) {
  return isValid(p->pdg_id());
}
template <>
inline bool isValid(const DecodedPID& p) {
  return p.pid() != 0 &&
         (isQuark(p) || isLepton(p) || isBoson(p) || isGlueball(p) ||
          isTrajectory(p.pid()) || isGenSpecific(p.pid()) || isDiquark(p) ||
          isBSM(p) || isHadron(p) || isNucleus(p) || isGeantino(p.pid()) ||
          isPythia8Specific(p));
}
template <>
inline bool isValid(const int& p) {
  if (!p)
    return false;
  if (std::abs(p) < 42)
    return true;
  if (isGenSpecific(p))
    return true;
  auto value_digits = DecodedPID(p);
  return isValid(value_digits);
}

template <class T>
inline int leadingQuark(const T& p) {
  return leadingQuark(p->pdg_id());
}
template <>
inline int leadingQuark(const DecodedPID& p) {
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

template <>
inline int leadingQuark(const int& p) {
  auto value_digits = DecodedPID(p);
  return leadingQuark(value_digits);
}

template <class T>
inline bool isLightHadron(const T& p) {
  auto lq = leadingQuark(p);
  return (lq == DQUARK || lq == UQUARK || lq == SQUARK) && isHadron(p);
}
template <class T>
inline bool isHeavyHadron(const T& p) {
  auto lq = leadingQuark(p);
  return (lq == CQUARK || lq == BQUARK || lq == TQUARK) && isHadron(p);
}
template <class T>
inline bool isStrangeHadron(const T& p) {
  return leadingQuark(p) == SQUARK && isHadron(p);
}
template <class T>
inline bool isCharmHadron(const T& p) {
  return leadingQuark(p) == CQUARK && isHadron(p);
}
template <class T>
inline bool isBottomHadron(const T& p) {
  return leadingQuark(p) == BQUARK && isHadron(p);
}
template <class T>
inline bool isTopHadron(const T& p) {
  return leadingQuark(p) == TQUARK && isHadron(p);
}

template <class T>
inline bool isLightMeson(const T& p) {
  auto lq = leadingQuark(p);
  return (lq == DQUARK || lq == UQUARK || lq == SQUARK) && isMeson(p);
}
template <class T>
inline bool isHeavyMeson(const T& p) {
  auto lq = leadingQuark(p);
  return (lq == CQUARK || lq == BQUARK || lq == TQUARK) && isMeson(p);
}
template <class T>
inline bool isStrangeMeson(const T& p) {
  return leadingQuark(p) == SQUARK && isMeson(p);
}
template <class T>
inline bool isCharmMeson(const T& p) {
  return leadingQuark(p) == CQUARK && isMeson(p);
}
template <class T>
inline bool isBottomMeson(const T& p) {
  return leadingQuark(p) == BQUARK && isMeson(p);
}
template <class T>
inline bool isTopMeson(const T& p) {
  return leadingQuark(p) == TQUARK && isMeson(p);
}

template <class T>
inline bool isCCbarMeson(const T& p) {
  return isCCbarMeson(p->pdg_id());
}
template <>
inline bool isCCbarMeson(const DecodedPID& p) {
  return leadingQuark(p) == CQUARK && isMeson(p) &&
         (*(p.second.rbegin() + 2)) == CQUARK &&
         (*(p.second.rbegin() + 1)) == CQUARK;
}
template <>
inline bool isCCbarMeson(const int& p) {
  return isCCbarMeson(DecodedPID(p));
}

template <class T>
inline bool isBBbarMeson(const T& p) {
  return isBBbarMeson(p->pdg_id());
}
template <>
inline bool isBBbarMeson(const DecodedPID& p) {
  return leadingQuark(p) == BQUARK && isMeson(p) &&
         (*(p.second.rbegin() + 2)) == BQUARK &&
         (*(p.second.rbegin() + 1)) == BQUARK;
}
template <>
inline bool isBBbarMeson(const int& p) {
  return isBBbarMeson(DecodedPID(p));
}

template <class T>
inline bool isLightBaryon(const T& p) {
  auto lq = leadingQuark(p);
  return (lq == DQUARK || lq == UQUARK || lq == SQUARK) && isBaryon(p);
}
template <class T>
inline bool isHeavyBaryon(const T& p) {
  auto lq = leadingQuark(p);
  return (lq == CQUARK || lq == BQUARK || lq == TQUARK) && isBaryon(p);
}
template <class T>
inline bool isStrangeBaryon(const T& p) {
  return leadingQuark(p) == SQUARK && isBaryon(p);
}
template <class T>
inline bool isCharmBaryon(const T& p) {
  return leadingQuark(p) == CQUARK && isBaryon(p);
}
template <class T>
inline bool isBottomBaryon(const T& p) {
  return leadingQuark(p) == BQUARK && isBaryon(p);
}
template <class T>
inline bool isTopBaryon(const T& p) {
  return leadingQuark(p) == TQUARK && isBaryon(p);
}

// APID: This function selects B-Hadrons which predominantly decay weakly.
// (Commonly used definition in GeneratorFilters package.) 5[1-4]1 L = J = 0, S
// = 0 5[1-5][1-4]2 J = 1/2, n_r = 0, n_L =0
template <class T>
inline bool isWeaklyDecayingBHadron(const T& p) {
  return isWeaklyDecayingBHadron(p->pdg_id());
}
template <>
inline bool isWeaklyDecayingBHadron(const int& p) {
  const int pid = std::abs(p);
  return (pid == 511 ||   // B0
          pid == 521 ||   // B+
          pid == 531 ||   // B_s0
          pid == 541 ||   // B_c+
          pid == 5122 ||  // Lambda_b0
          pid == 5132 ||  // Xi_b-
          pid == 5232 ||  // Xi_b0
          pid == 5112 ||  // Sigma_b-
          pid == 5212 ||  // Sigma_b0
          pid == 5222 ||  // Sigma_b+
          pid == 5332 ||  // Omega_b-
          pid == 5142 ||  // Xi_bc0
          pid == 5242 ||  // Xi_bc+
          pid == 5412 ||  // Xi'_bc0
          pid == 5422 ||  // Xi'_bc+
          pid == 5342 ||  // Omega_bc0
          pid == 5432 ||  // Omega'_bc0
          pid == 5442 ||  // Omega_bcc+
          pid == 5512 ||  // Xi_bb-
          pid == 5522 ||  // Xi_bb0
          pid == 5532 ||  // Omega_bb-
          pid == 5542);   // Omega_bbc0
}
template <>
inline bool isWeaklyDecayingBHadron(const DecodedPID& p) {
  return isWeaklyDecayingBHadron(p.pid());
}

// APID: This function selects C-Hadrons which predominantly decay weakly.
// (Commonly used definition in GeneratorFilters package.) 4[1-3]1 L = J = 0, S
// = 0 4[1-4][1-3]2 J = 1/2, n_r = 0, n_L =0 NB Omitting pid = 4322 (Xi'_C+) a
// this undergoes an EM rather than weak decay.  (There was an old version of
// Herwig that decayed it weakly, but this was fixed in Herwig 7.)
template <class T>
inline bool isWeaklyDecayingCHadron(const T& p) {
  return isWeaklyDecayingCHadron(p->pdg_id());
}
template <>
inline bool isWeaklyDecayingCHadron(const int& p) {
  const int pid = std::abs(p);
  return (pid == 411 ||   // D+
          pid == 421 ||   // D0
          pid == 431 ||   // Ds+
          pid == 4122 ||  // Lambda_c+
          pid == 4132 ||  // Xi_c0
          pid == 4232 ||  // Xi_c+
          pid == 4212 ||  // Xi_c0
          pid == 4332 ||  // Omega_c0
          pid == 4412 ||  // Xi_cc+
          pid == 4422 ||  // Xi_cc++
          pid == 4432);   // Omega_cc+
}
template <>
inline bool isWeaklyDecayingCHadron(const DecodedPID& p) {
  return isWeaklyDecayingCHadron(p.pid());
}

template <class T>
inline int charge3(const T& p) {
  return charge3(p->pdg_id());
}
template <class T>
inline double fractionalCharge(const T& p) {
  return fractionalCharge(p->pdg_id());
}
template <class T>
inline double charge(const T& p) {
  if (isGenericMultichargedParticle(
          p))  // BSM multi-charged particles might have a fractional charge
               // that's not a multiple of 1/3
    return fractionalCharge(p);
  else
    return 1.0 * charge3(p) / 3.0;
}
template <class T>
inline double threeCharge(const T& p) {
  return charge3(p);
}
template <class T>
inline bool isCharged(const T& p) {
  return charge3(p) != 0;
}

template <>
inline int charge3(const DecodedPID& p) {
  auto ap = std::abs(p.pid());
  if (ap < TABLESIZE)
    return p.pid() > 0 ? triple_charge.at(ap) : -triple_charge.at(ap);
  if (ap == K0)
    return 0;
  if (ap == GEANTINO0)
    return 0;
  if (ap == GEANTINOPLUS)
    return p.pid() > 0 ? 3 : -3;
  if (ap == MAVTOP)
    return p.pid() > 0 ? 2 : -2;
  std::size_t nq = 0;
  int sign = 1;
  int signmult = 1;
  int result = 0;
  bool classified = false;
  if (!classified && isMeson(p)) {
    classified = true;
    nq = 2;
    if ((*(p.second.rbegin() + 2)) == 2 || (*(p.second.rbegin() + 2)) == 4) {
      sign = -1;
    }
    signmult = -1;
  }
  if (!classified && isDiquark(p)) {
    return triple_charge.at(p(0)) + triple_charge.at(p(1));
  }
  if (!classified && isBaryon(p)) {
    classified = true;
    nq = 3;
  }
  if (!classified && isTetraquark(p)) {
    return triple_charge.at(p(3)) + triple_charge.at(p(4)) -
           triple_charge.at(p(6)) - triple_charge.at(p(7));
  }
  if (!classified && isPentaquark(p)) {
    return triple_charge.at(p(3)) + triple_charge.at(p(4)) +
           triple_charge.at(p(5)) + triple_charge.at(p(6)) -
           triple_charge.at(p(7));
  }
  if (!classified && isNucleus(p)) {
    return 3 * numberOfProtons(p);
  }
  if (!classified && isSUSY(p)) {
    nq = 0;
    auto pp = p.shift(1);
    if (pp.ndigits() < 3) {
      return charge3(pp);
    }  // super-partners of fundamental particles
    if (pp(0) == COMPOSITEGLUON) {
      if (pp(1) == COMPOSITEGLUON) {
        return 0;
      }  // R-Glueballs
      if (pp.ndigits() == 4 || pp.ndigits() == 5) {
        pp = pp.shift(1);  // Remove gluino
      }
    }
    if (pp.ndigits() == 3) {
      classified = true;
      nq = 2;
      if (p.last() % 2 == 0) {
        sign = -1;
      }
      signmult = -1;
    }  // states with squark-antiquark or quark-anti-quark
    if (pp.ndigits() == 4) {
      classified = true;
      nq = 3;
    }  // states with squark-quark-quark or quark-quark-quark
  }
  if (!classified && isMonopole(p)) {
    /// Codes 411nq1nq2 nq3 0  are then used when the magnetic and electrical
    /// charge sign agree and 412nq1nq2 nq3 0
    ///  when they disagree, with the overall sign of the particle set by the
    ///  magnetic charge.
    result = 3 * (p(3) * 100 + p(4) * 10 + p(5));
    return ((p.pid() > 0 && p(2) == 1) || (p.pid() < 0 && p(2) == 2)) ? result
                                                                      : -result;
  }
  if (!classified && isGenericMultichargedParticle(p)) {
    double abs_charge = 0.0;
    if (p(0) == 1)
      abs_charge = p(3) * 100. + p(4) * 10. + p(5) * 1 +
                   p(6) * 0.1;  // multi-charged particle PDG ID is +/-100XXXY0,
                                // where the charge is XXX.Y
    if (p(0) == 2)
      abs_charge =
          (p(3) * 10. + p(4)) /
          (p(5) * 10.0 + p(6));  // multi-charged particle PDG ID is
                                 // +/-200XXYY0, where the charge is XX/YY
    int abs_threecharge = static_cast<int>(std::round(
        abs_charge *
        3.));  // the multi-charged particles might have a fractional charge
               // that's not a multiple of 1/3, in that case round to the
               // closest multiple of 1/3 for charge3 and threecharge
    return p.pid() > 0 ? abs_threecharge : -1 * abs_threecharge;
  }
  for (auto r = p.second.rbegin() + 1; r != p.second.rbegin() + 1 + nq; ++r) {
    result += triple_charge.at(*r) * sign;
    sign *= signmult;
  }
  return p.pid() > 0 ? result : -result;
}
template <>
inline int charge3(const int& p) {
  int ap = std::abs(p);
  if (ap < TABLESIZE)
    return p > 0 ? triple_charge.at(ap) : -triple_charge.at(ap);
  auto value_digits = DecodedPID(p);
  return charge3(value_digits);
}

template <class T>
inline bool isNeutral(const T& p) {
  return p->pdg_id() != 0 && charge3(p) == 0;
}
template <>
inline bool isNeutral(const DecodedPID& p) {
  return p.pid() != 0 && charge3(p) == 0;
}
template <>
inline bool isNeutral(const int& p) {
  auto value_digits = DecodedPID(p);
  return isNeutral(value_digits);
}

template <>
inline double fractionalCharge(const DecodedPID& p) {
  if (!isGenericMultichargedParticle(p))
    return 1.0 * charge3(p) /
           3.0;  // this method is written for multi-charged particles, still
                 // make sure other cases are handled properly
  double abs_charge = 0;
  if (p(0) == 1)
    abs_charge = p(3) * 100. + p(4) * 10. + p(5) * 1 +
                 p(6) * 0.1;  // multi-charged particle PDG ID is +/-100XXXY0,
                              // where the charge is XXX.Y
  if (p(0) == 2)
    abs_charge =
        (p(3) * 10. + p(4)) /
        (p(5) * 10.0 + p(6));  // multi-charged particle PDG ID is +/-200XXYY0,
                               // where the charge is XX/YY
  return p.pid() > 0 ? abs_charge : -1 * abs_charge;
}
template <>
inline double fractionalCharge(const int& p) {
  auto value_digits = DecodedPID(p);
  return fractionalCharge(value_digits);
}

// APID: Including Z' and Z'' as EM interacting.
template <class T>
inline bool isEMInteracting(const T& p) {
  return isEMInteracting(p->pdg_id());
}
template <>
inline bool isEMInteracting(const int& p) {
  return (isPhoton(p) || isZ(p) || p == ZPRIME || p == ZDBLPRIME ||
          std::abs(charge(p)) > std::numeric_limits<double>::epsilon() ||
          isMonopole(p));
}

template <class T>
inline bool isParton(const T& p) {
  return isQuark(p) || isGluon(p);
}

// APID: Intended to return 2J
// Useful for G4ParticleDefinition constructor
template <class T>
inline int spin2(const T& p) {
  return spin2(p->pdg_id());
}
template <>
inline int spin2(const DecodedPID& p) {
  if (isSUSY(p)) {
    auto pp = p.shift(1);
    auto ap = std::abs(pp.pid());
    if (ap < TABLESIZE) {
      return std::abs(double_spin.at(ap) - 1);
    }  // sparticles (0->1, 1 -> 0,  2->1,  4->3)
    return p.last() - 1;  // R-Hadrons (p.last() == 2J +1)
  }
  auto ap = std::abs(p.pid());
  if (ap == K0S) {
    return 0;
  }
  if (ap == K0L) {
    return 0;
  }
  if (ap == MAVTOP) {
    return 1;
  }  // TODO check this
  if (ap == DARKPHOTON) {
    return 2;
  }  // TODO check this
  if (ap < TABLESIZE) {
    return double_spin.at(ap);
  }  // fundamental particles
  if (isHadron(p)) {
    return p.last() - 1;
  }  // Hadrons (p.last == 2J+1 - special cases handled above)
  if (isMonopole(p)) {
    return 0;
  }  // PDG 11i - For now no spin information is provided. Also matches the
     // definition in the G4Extensions/Monopole package.
  if (isGenericMultichargedParticle(p)) {
    return 0;
  }  // APID Matches the definition in the G4Extensions/Monopole package.
  if (isNucleus(p)) {
    return 1;
  }  // TODO need to explicitly deal with nuclei
  return p.last() > 0 ? 1 : 0;  //  Anything else - best guess
}
template <>
inline int spin2(const int& p) {
  auto value_digits = DecodedPID(p);
  return spin2(value_digits);
}

template <class T>
inline double spin(const T& p) {
  return spin(p->pdg_id());
}
template <>
inline double spin(const DecodedPID& p) {
  return 1.0 * spin2(p) / 2.0;
}
template <>
inline double spin(const int& p) {
  auto value_digits = DecodedPID(p);
  return spin(value_digits);
}

template <class T>
inline bool isRHadron(const T& p) {
  return isRHadron(p->pdg_id());
}
template <>
inline bool isRHadron(const DecodedPID& p) {
  return (isRBaryon(p) || isRMeson(p) || isRGlueball(p));
}
template <>
inline bool isRHadron(const int& p) {
  auto value_digits = DecodedPID(p);
  return isRHadron(value_digits);
}

// APID: Returns an unordered list of the quarks contained by the current
// particle
template <class T>
inline std::vector<int> containedQuarks(const T& p) {
  return containedQuarks(p->pdg_id());
}
template <>
inline std::vector<int> containedQuarks(const int& p) {
  auto pp = DecodedPID(p);
  std::vector<int> quarks;
  if (isQuark(pp.pid())) {
    quarks.push_back(std::abs(pp.pid()));
  } else if (isDiquark(pp)) {
    quarks.push_back(pp(0));
    quarks.push_back(pp(1));
  } else if (isMeson(pp)) {
    quarks.push_back(*(pp.second.rbegin() + 1));
    quarks.push_back(*(pp.second.rbegin() + 2));
  } else if (isBaryon(pp)) {
    for (std::size_t digit = 1; digit < 4; ++digit) {
      quarks.push_back(*(pp.second.rbegin() + digit));
    }
  } else if (isTetraquark(pp)) {
    for (std::size_t digit = 1; digit < 5; ++digit) {
      quarks.push_back(*(pp.second.rbegin() + digit));
    }
  } else if (isPentaquark(pp)) {
    for (std::size_t digit = 1; digit < 6; ++digit) {
      quarks.push_back(*(pp.second.rbegin() + digit));
    }
  } else if (isNucleus(pp)) {
    const int A = std::abs(baryonNumber3(pp) / 3);
    const int Z = std::abs(numberOfProtons(pp));
    const int L = std::abs(numberOfLambdas(pp));
    const int n_uquarks = A + Z;
    const int n_dquarks = 2 * A - Z - L;
    const int n_squarks = L;
    quarks.reserve(3 * A);
    quarks.insert(quarks.end(), n_dquarks, 1);
    quarks.insert(quarks.end(), n_uquarks, 2);
    quarks.insert(quarks.end(), n_squarks, 3);
  } else if (isSUSY(pp)) {  // APID SUSY case
    pp = pp.shift(1);
    if (pp.ndigits() > 1) {  // skip squarks
      if (pp.ndigits() == 3) {
        pp = DecodedPID(pp(1));
      }  // Handle ~q qbar pairs
      if (pp.ndigits() > 3) {
        pp = pp.shift(1);
      }  // Drop gluinos and squarks
      return containedQuarks(pp.pid());
    }
  }
  return quarks;
}
template <>
inline std::vector<int> containedQuarks(const DecodedPID& p) {
  return containedQuarks(p.pid());
}

template <class T>
inline bool isStrongInteracting(const T& p) {
  return isStrongInteracting(p->pdg_id());
}
template <>
inline bool isStrongInteracting(const int& p) {
  return (isGluon(p) || isQuark(p) || isDiquark(p) || isGlueball(p) ||
          isLeptoQuark(p) || isHadron(p) || isRHadron(p));
}  // APID: Glueballs and R-Hadrons are also strong-interacting

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

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/HepMCHelpers.h#L33
bool isInteracting(int pdg) {
  return isStrongInteracting(pdg) || isEMInteracting(pdg) || isGeantino(pdg);
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
