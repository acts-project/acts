// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

#include "Acts/Utilities/Units.hpp"

namespace Acts {
namespace constants {

using namespace ::Acts::UnitLiterals;

// See (1) table 32.1
// K/A*Z = 0.5 * 30.7075MeV/(g/mm2) * Z/A * rho[g/mm3]
constexpr double ka_BetheBloch = 30.7075_MeV * 1_mm * 1_mm;

/// Ionisation potential
/// Ionization - Bethe-Bloch
/// See ATL-SOFT-PUB-2008-003 equation (4)
/// 16 eV * Z**0.9
constexpr double eionisation = 16_eV;

/// Plasma energy [ eV ]
constexpr double eplasma = 28.816_eV;

/// Fine structure constexprant
constexpr double alpha = 1. / 137.;

/// Multiple scattering parameters for MIPS [ in MeV]
constexpr double main_RutherfordScott = 13.6_MeV;
constexpr double log_RutherfordScott = 0.038;

/// Multiple scattering parameters for electrons [ in MeV ]
constexpr double main_RossiGreisen = 17.5_MeV;
constexpr double log_RossiGreisen = 0.125;

/// Electron mass [ in MeV ]
constexpr double me = 0.51099891_MeV;

/// 4 * scaling factor for Landau FWHM to gaussian sigma
const double landau2gauss = 2. / (std::sqrt(2. * std::log(2.)));

}  // end of namespace constants
}  // end of namespace Acts
