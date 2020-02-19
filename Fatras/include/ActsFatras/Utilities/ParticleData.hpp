// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iosfwd>
#include <string_view>

#include "Acts/Utilities/PdgParticle.hpp"

namespace ActsFatras {

/// Find the charge for a given PDG particle number.
///
/// @return Charge in native units or NaN if not available.
float findCharge(Acts::PdgParticle pdg);

/// Find the mass for a given PDG particle number.
///
/// @return Mass in native units or zero if not available.
float findMass(Acts::PdgParticle pdg);

/// Find a descriptive particle name for a given PDG particle number.
///
/// @return Particle name or empty if not available.
std::string_view findName(Acts::PdgParticle pdg);

}  // namespace ActsFatras

namespace Acts {

/// Print PDG particle numbers with a descriptive name.
///
/// @note This is a bit hacky since it modifies the namespace of
///       another library (ActsCore). Since the core library does not need
///       to contain the full particle data table, it can only be defined
///       here. It also only extends the output and should not have any
///       side effects. We are probably fine.
std::ostream& operator<<(std::ostream& os, PdgParticle pdg);

}  // namespace Acts
