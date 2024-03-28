// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <iosfwd>

namespace ActsExamples {

/// Particle outcome identifier.
///
/// Encodes the status of the particle during the Geant4 simulation
enum class Geant4ParticleStatus : uint32_t {
  Undefined = 0,
  KilledAbsorbed = 1,
  KilledVolume = 2,
  KilledTime = 3,
  KilledSecondary = 4,
};

std::ostream &operator<<(std::ostream &os, Geant4ParticleStatus particleStatus);

}  // namespace ActsExamples
