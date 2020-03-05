// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>

namespace Acts {

/// Symbolic values for commonly used PDG particle numbers.
enum PdgParticle : int32_t {
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
  eNeutron = 2112,
  eAntiNeutron = -eNeutron,
  eProton = 2212,
  eAntiProton = -eProton,
};

/// Convert an anti-particle to its particle and leave particles as-is.
constexpr inline PdgParticle makeAbsolutePdgParticle(PdgParticle pdg) {
  const auto value = static_cast<int32_t>(pdg);
  return static_cast<PdgParticle>((0 <= value) ? value : -value);
}

}  // namespace Acts
