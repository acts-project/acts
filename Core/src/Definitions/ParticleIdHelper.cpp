// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/ParticleData.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <ostream>
#include <ranges>
#include <utility>
#include <vector>

// This file is based on the following files, with modifications:
//
// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h
//
// These files are licensed under Apache License 2.0

namespace Acts {
namespace {

class DecodedPID : public std::pair<int, std::vector<int>> {
 public:
  explicit DecodedPID(const int& p) {
    this->first = p;
    this->second.reserve(10);
    int ap = std::abs(p);
    for (; ap != 0; ap /= 10) {
      this->second.push_back(ap % 10);
    }
    std::ranges::reverse(this->second);
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

const int TABLESIZE = 100;
const std::array<int, TABLESIZE> triple_charge = {
    +0, -1, +2, -1, +2, -1, +2, -1, +2, +0, +0, -3, +0, -3, +0, -3, +0,
    -3, +0, +0, +0, +0, +0, +0, +3, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +3, +0, +0, +3, +6, +0, +0, +0, -1, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0};

const int DQUARK = 1;
const int UQUARK = 2;
const int SQUARK = 3;
const int CQUARK = 4;
const int BQUARK = 5;
const int TQUARK = 6;

const int ELECTRON = 11;
const int MUON = 13;
const int TAU = 15;

const int GLUON = 21;
// APID: 9 rather than 21 is used to denote a gluon/gluino in composite states.
// (From PDG 11g)
const int COMPOSITEGLUON = 9;
const int PHOTON = 22;
const int Z0BOSON = 23;
const int ZPRIME = 32;     // Z′/Z^0_2
const int ZDBLPRIME = 33;  // Z′′/Z^0_3
const int LEPTOQUARK = 42;

/// PDG Ids for Mavtop madgraph UFO model found under DarkX. The
/// mavtop is a vector-like top partner with coupling to a dark photon.
/// Theory paper: https://arxiv.org/abs/1904.05893
/// Pheno paper: https://arxiv.org/pdf/2112.08425
const int MAVTOP = 60001;

const int K0L = 130;

const int K0S = 310;
const int K0 = 311;
const int PROTON = 2212;

/// PDG rule 10:
/// Codes 81–100 are reserved for generator-specific pseudoparticles and
/// concepts. Codes 901–930, 1901–1930, 2901–2930, and 3901–3930 are for
/// additional components of Standard Modelparton distribution functions, where
/// the latter three ranges are intended to distinguish left/right/ longitudinal
/// components. Codes 998 and 999 are reserved for GEANT tracking pur-poses.
const int GEANTINOPLUS = 998;
const int GEANTINO0 = 999;

/// PDG rule 2:
/// Quarks and leptons are numbered consecutively starting from 1 and 11
/// respectively; to dothis they are first ordered by family and within
/// families by weak isospin.
/// APID: the fourth generation quarks are quarks.
bool isQuark(int p) {
  return p != 0 && (std::abs(p) <= 8 || std::abs(p) == MAVTOP);
}

bool isSMQuark(int p) {
  return p != 0 && std::abs(p) <= TQUARK;
}

/// PDG rule 4
/// Diquarks have 4-digit numbers with nq1 >= nq2 and nq3 = 0
/// APID: the diquarks with fourth generation are not diquarks
bool isDiquark(const DecodedPID& p) {
  if (p.ndigits() == 4 && p(0) >= p(1) && p(2) == 0 && p.last() % 2 == 1 &&
      p.max_digit(2, 4) <= TQUARK) {
    return true;
  }
  return false;
}

/// Table 43.1
///  PDG rule 5a:
///  The numbers specifying the meson's quark content conform to the convention
///  nq1= 0 and nq2 >= nq3. The special case K0L is the sole exception to this
///  rule. PDG rule 5C: The special numbers 310 and 130 are given to the K0S and
///  K0L respectively. APID: The special code K0 is used when a generator uses
///  K0S/K0L
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

/// Table 43.2
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
bool isTetraquark(const DecodedPID& p) {
  return (p.ndigits() == 9 && p(0) == 1 && p(5) == 0 &&
          p.max_digit(1, 3) <= 6 && p.min_digit(1, 3) > 0 &&
          p.max_digit(1 + 3, 3 + 3) <= 6 && p.min_digit(1 + 3, 3 + 3) > 0 &&
          (p(3) >= p(4) && p(6) >= p(7)) &&
          ((p(3) > p(6)) || (p(3) == p(6) && (p(4) >= p(7)))));
}

/// PDG rule 15
/// The 9-digit penta-quark codes are±1nrnLnq1nq2nq3nq4nq5nJ, sorted such
/// thatnq1≥nq2≥nq3≥nq4. In the particle the first four are quarks and the fifth
/// an antiquark while t
/// heopposite holds in the antiparticle, which is given with a negative sign.
/// Thenr,nL, andnJnumbers have the same meaning as for ordinary hadrons.
bool isPentaquark(const DecodedPID& p) {
  return (p.ndigits() == 9 && p(0) == 1 && p.max_digit(1, 6) <= 6 &&
          p.min_digit(1, 6) > 0 &&
          (p(3) >= p(4) && p(4) >= p(5) && p(5) >= p(6)));
}

// APID Mesons, Baryons, Tetraquarks and Pentaquarks are Hadrons
bool isHadron(const DecodedPID& p) {
  return isMeson(p) || isBaryon(p) || isTetraquark(p) || isPentaquark(p);
}

/// PDG rule 9:
/// Two-digit numbers in the range 21–30 are provided for the Standard
/// Model gauge and Higgs bosons.
/// PDG rule 11b:
/// The graviton and the boson content of a two-Higgs-doublet scenario
/// and of additional SU(2)×U(1) groups are found in the range 31–40.

bool isGluon(int p) {
  return p == GLUON;
}

bool isPhoton(int p) {
  return p == PHOTON;
}

bool isZ(int p) {
  return p == Z0BOSON;
}

/// PDG rule 11c:
/// "One-of-a-kind" exotic particles are assigned numbers in the range
/// 41–80. The subrange 61-80 can be used for new heavier fermions in
/// generic models, where partners to the SM fermions would have codes
/// oﬀset by 60. If required, however, other assignments could be
/// made.
bool isLeptoQuark(int p) {
  return std::abs(p) == LEPTOQUARK;
}

/// Main Table
/// for MC internal use 81–100,901–930,998-999,1901–1930,2901–2930, and
/// 3901–3930
bool isGenSpecific(int p) {
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

bool isGeantino(int p) {
  return (std::abs(p) == GEANTINO0 || std::abs(p) == GEANTINOPLUS);
}

/// APID: Definition of Glueballs: SM glueballs 99X (X=1,5), 999Y (Y=3,7)
bool isGlueball(const DecodedPID& p) {
  if (p.ndigits() > 4) {
    return false;  // APID avoid classifying R-Glueballs as SM Glueballs
  }
  return ((p.ndigits() == 3 && p(0) == COMPOSITEGLUON &&
           p(1) == COMPOSITEGLUON && (p(2) == 1 || p(2) == 5)) ||
          (p.ndigits() == 4 && p(0) == COMPOSITEGLUON &&
           p(1) == COMPOSITEGLUON && p(2) == COMPOSITEGLUON &&
           (p(3) == 3 || p(3) == 7)));
}

/// PDG rule 11d
/// Fundamental supersymmetric particles are identified by adding a nonzero n to
/// the particle number. The superpartner of a boson or a left-handed fermion
/// has n = 1 while the superpartner of a right-handed fermion has n = 2. When
/// mixing occurs, such as between the winos and charged Higgsinos to give
/// charginos, or between left and right sfermions, the lighter physical state
/// is given the smaller basis state number.
bool isSUSY(const DecodedPID& p) {
  return (p.ndigits() == 7 && (p(0) == 1 || p(0) == 2) &&
          !isGenSpecific(p.shift(2).pid()));
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
bool isRGlueball(const DecodedPID& p) {
  if (p.ndigits() != 7 || p(0) != 1) {
    return false;
  }
  auto pp = p.shift(1);
  return ((pp.ndigits() == 3 && pp(0) == COMPOSITEGLUON &&
           pp(1) == COMPOSITEGLUON && (pp(2) == 1 || pp(2) == 3)) ||
          (pp.ndigits() == 4 && pp(0) == COMPOSITEGLUON &&
           pp(1) == COMPOSITEGLUON && pp(2) == COMPOSITEGLUON &&
           (pp(3) == 1 || pp(3) == 5)));
}

// APID Define R-Mesons as gluino-quark-antiquark and squark-antiquark bound
// states (ignore 4th generation squarks/quarks) NB Current models only allow
// gluino-quark-antiquark, stop-antiquark and sbottom-antiquark states
bool isRMeson(const DecodedPID& p) {
  if (!(p.ndigits() == 7 && (p(0) == 1 || p(0) == 2))) {
    return false;
  }
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

// APID Define R-Baryons as gluino-quark-quark-quark and squark-quark-quark
// bound states (ignore 4th generation squarks/quarks) NB Current models only
// allow gluino-quark-quark-quark, stop-quark-quark and sbottom-quark-quark
// states
bool isRBaryon(const DecodedPID& p) {
  if (!(p.ndigits() == 7 && (p(0) == 1 || p(0) == 2))) {
    return false;
  }
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

/// PDG rule 11i
/// Magnetic monopoles and dyons are assumed to have one unit of Dirac monopole
/// charge and a variable integer number nq1nq2 nq3 units of electric charge.
/// Codes 411nq1nq2 nq3 0 are then used when the magnetic and electrical charge
/// sign agree and 412nq1nq2 nq3 0 when they disagree, with the overall sign of
/// the particle set by the magnetic charge. For now no spin information is
/// provided.
bool isMonopole(const DecodedPID& p) {
  return (p.ndigits() == 7 && p(0) == 4 && p(1) == 1 &&
          (p(2) == 1 || p(2) == 2) && p(6) == 0);
}

/// In addition, there is a need to identify "Q-ball" and similar very exotic
/// (multi-charged) particles which may have large, non-integer charge. These
/// particles are assigned the ad-hoc numbering +/-100XXXY0, where the charge is
/// XXX.Y. or +/-200XXYY0, where the charge is XX/YY. The case of +/-200XXYY0 is
/// legacy, see https://gitlab.cern.ch/atlas/athena/-/merge_requests/25862 Note
/// that no other quantum numbers besides the charge are considered for these
/// generic multi-charged particles (e.g. isSUSY() is false for them). Such a
/// model was used in previous Run-1 (1301.5272,1504.04188) and Run-2
/// (1812.03673,2303.13613) ATLAS searches.
bool isGenericMultichargedParticle(const DecodedPID& p) {
  return (p.ndigits() == 8 && (p(0) == 1 || p(0) == 2) && p(1) == 0 &&
          p(2) == 0 && p(7) == 0);
}

/// PDG rule 16
/// Nuclear codes are given as 10-digit numbers ±10LZZZAAAI.
/// For a (hyper)nucleus consisting of n_p protons, n_n neutrons and
/// n_Λ Λ's:
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
bool isNucleus(const DecodedPID& p) {
  if (std::abs(p.pid()) == PROTON) {
    return true;
  }
  return (p.ndigits() == 10 && p(0) == 1 && p(1) == 0);
}

int numberOfProtons(const DecodedPID& p) {
  if (std::abs(p.pid()) == PROTON) {
    return (p.pid() > 0) ? 1 : -1;
  }
  if (isNucleus(p)) {
    const int result = p(5) + 10 * p(4) + 100 * p(3);
    return (p.pid() > 0) ? result : -result;
  }
  return 0;
}

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

bool isLightMeson(const DecodedPID& p) {
  auto lq = leadingQuark(p);
  return (lq == DQUARK || lq == UQUARK || lq == SQUARK) && isMeson(p);
}

bool isStrangeMeson(const DecodedPID& p) {
  return leadingQuark(p) == SQUARK && isMeson(p);
}

bool isCharmMeson(const DecodedPID& p) {
  return leadingQuark(p) == CQUARK && isMeson(p);
}

bool isBottomMeson(const DecodedPID& p) {
  return leadingQuark(p) == BQUARK && isMeson(p);
}

inline bool isCCbarMeson(const DecodedPID& p) {
  return leadingQuark(p) == CQUARK && isMeson(p) &&
         (*(p.second.rbegin() + 2)) == CQUARK &&
         (*(p.second.rbegin() + 1)) == CQUARK;
}

inline bool isBBbarMeson(const DecodedPID& p) {
  return leadingQuark(p) == BQUARK && isMeson(p) &&
         (*(p.second.rbegin() + 2)) == BQUARK &&
         (*(p.second.rbegin() + 1)) == BQUARK;
}

bool isLightBaryon(const DecodedPID& p) {
  auto lq = leadingQuark(p);
  return (lq == DQUARK || lq == UQUARK || lq == SQUARK) && isBaryon(p);
}

bool isStrangeBaryon(const DecodedPID& p) {
  return leadingQuark(p) == SQUARK && isBaryon(p);
}

bool isCharmBaryon(const DecodedPID& p) {
  return leadingQuark(p) == CQUARK && isBaryon(p);
}

bool isBottomBaryon(const DecodedPID& p) {
  return leadingQuark(p) == BQUARK && isBaryon(p);
}

double fractionalCharge(const DecodedPID& p);
int charge3(const DecodedPID& p);

double charge(const DecodedPID& p) {
  if (isGenericMultichargedParticle(
          p)) {  // BSM multi-charged particles might have a fractional charge
                 // that's not a multiple of 1/3
    return fractionalCharge(p);
  } else {
    return 1.0 * charge3(p) / 3.0;
  }
}

int charge3(const DecodedPID& p) {
  auto ap = std::abs(p.pid());
  if (ap < TABLESIZE) {
    return copySign(triple_charge.at(ap), p.pid());
  }
  if (ap == K0) {
    return 0;
  }
  if (ap == GEANTINO0) {
    return 0;
  }
  if (ap == GEANTINOPLUS) {
    return copySign(3, p.pid());
  }
  if (ap == MAVTOP) {
    return copySign(2, p.pid());
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
    if (p(0) == 1) {
      abs_charge = p(3) * 100. + p(4) * 10. + p(5) * 1 +
                   p(6) * 0.1;  // multi-charged particle PDG ID is
                                // +/-100XXXY0, where the charge is XXX.Y
    }
    if (p(0) == 2) {
      abs_charge = (p(3) * 10. + p(4)) /
                   (p(5) * 10.0 + p(6));  // multi-charged particle PDG ID is
      // +/-200XXYY0, where the charge is XX/YY
    }
    int abs_threecharge = static_cast<int>(std::round(
        abs_charge *
        3.));  // the multi-charged particles might have a fractional charge
               // that's not a multiple of 1/3, in that case round to the
               // closest multiple of 1/3 for charge3 and threecharge
    return copySign(abs_threecharge, p.pid());
  }
  for (auto r = p.second.rbegin() + 1; r != p.second.rbegin() + 1 + nq; ++r) {
    result += triple_charge.at(*r) * sign;
    sign *= signmult;
  }
  return copySign(result, p.pid());
}
double fractionalCharge(const DecodedPID& p) {
  if (!isGenericMultichargedParticle(p)) {
    return 1.0 * charge3(p) /
           3.0;  // this method is written for multi-charged particles, still
                 // make sure other cases are handled properly
  }
  double abs_charge = 0;
  if (p(0) == 1) {
    abs_charge = p(3) * 100. + p(4) * 10. + p(5) * 1 +
                 p(6) * 0.1;  // multi-charged particle PDG ID is +/-100XXXY0,
  }
  // where the charge is XXX.Y
  if (p(0) == 2) {
    abs_charge =
        (p(3) * 10. + p(4)) /
        (p(5) * 10.0 + p(6));  // multi-charged particle PDG ID is +/-200XXYY0,
  }
  // where the charge is XX/YY
  return copySign(abs_charge, p.pid());
}

// APID: Including Z' and Z'' as EM interacting.
bool isEMInteracting(const DecodedPID& p) {
  return (isPhoton(p.pid()) || isZ(p.pid()) || p.pid() == ZPRIME ||
          p.pid() == ZDBLPRIME ||
          std::abs(charge(p)) > std::numeric_limits<double>::epsilon() ||
          isMonopole(p));
}

bool isRHadron(const DecodedPID& p) {
  return (isRBaryon(p) || isRMeson(p) || isRGlueball(p));
}

bool isStrongInteracting(const DecodedPID& p) {
  return (isGluon(p.pid()) || isQuark(p.pid()) || isDiquark(p) ||
          isGlueball(p) || isLeptoQuark(p.pid()) || isHadron(p) ||
          isRHadron(p));
}  // APID: Glueballs and R-Hadrons are also strong-interacting

}  // namespace

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L325
bool ParticleIdHelper::isHadron(PdgParticle pdg) {
  DecodedPID p(pdg);
  return isMeson(p) || isBaryon(p) || isTetraquark(p) || isPentaquark(p);
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/AtlasPID.h#L180
bool ParticleIdHelper::isLepton(PdgParticle pdg) {
  auto sp = std::abs(pdg);
  return sp >= 11 && sp <= 18;
}

bool ParticleIdHelper::isMuon(PdgParticle pdg) {
  return std::abs(pdg) == MUON;
}

bool ParticleIdHelper::isElectron(PdgParticle pdg) {
  return std::abs(pdg) == ELECTRON;
}

bool ParticleIdHelper::isPhoton(PdgParticle pdg) {
  return std::abs(pdg) == PHOTON;
}

bool ParticleIdHelper::isTau(PdgParticle pdg) {
  return std::abs(pdg) == TAU;
}

HadronType ParticleIdHelper::hadronType(PdgParticle pdg) {
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
bool ParticleIdHelper::isQuark(PdgParticle pdg) {
  return pdg != 0 && (std::abs(pdg) <= 8 || std::abs(pdg) == MAVTOP);
}

// https://gitlab.cern.ch/atlas/athena/-/blob/b0898f93585c4eec97550a9e0a16f5b6e6b6b973/Generators/TruthUtils/TruthUtils/HepMCHelpers.h#L33
bool ParticleIdHelper::isInteracting(PdgParticle pdg) {
  DecodedPID p(pdg);
  return isStrongInteracting(p) || isEMInteracting(p) || isGeantino(pdg);
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
  return os;
}

PdgParticle parsePdgParticle(const std::string& name) {
  if (name == "e-") {
    return eElectron;
  } else if (name == "e+") {
    return ePositron;
  } else if (name == "mu-") {
    return eMuon;
  } else if (name == "mu+") {
    return eAntiMuon;
  } else if (name == "tau-") {
    return eTau;
  } else if (name == "tau+") {
    return eAntiTau;
  } else if (name == "gamma") {
    return eGamma;
  } else if (name == "pi0") {
    return ePionZero;
  } else if (name == "pi+") {
    return ePionPlus;
  } else if (name == "pi-") {
    return ePionMinus;
  } else if (name == "K+") {
    return eKaonPlus;
  } else if (name == "K-") {
    return eKaonMinus;
  } else if (name == "n") {
    return eNeutron;
  } else if (name == "n~") {
    return eAntiNeutron;
  } else if (name == "p") {
    return eProton;
  } else if (name == "p~") {
    return eAntiProton;
  } else if (name == "Pb") {
    return eLead;
  } else {
    throw std::invalid_argument("Unknown particle name: " + name);
  }
}
}  // namespace Acts
