// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <iosfwd>

namespace ActsFatras {

/// Simulation outcome identifier.
///
/// Encodes the outcome of the particle after the simulation
enum class SimulationOutcome : std::uint32_t {
  Alive = 0,
  KilledInteraction = 1,
  KilledVolumeExit = 2,
  KilledTime = 3,
  KilledSecondaryParticle = 4,
};

/// Print simulation outcome to output stream
/// @param os Output stream
/// @param outcome Simulation outcome to print
/// @return Output stream
std::ostream &operator<<(std::ostream &os, SimulationOutcome outcome);

}  // namespace ActsFatras
