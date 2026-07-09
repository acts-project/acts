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

/// Generation process type identifier.
///
/// Encodes the type of process that generated a particle.
///
/// Fatras' own physics modules only ever set eDecay, ePhotonConversion and
/// eBremsstrahlung. The remaining values (eNuclearInteraction, eIonisation,
/// eOther) are not produced by Fatras (as of June 2026); they exist to
/// "translate" the finer production-process information available from the full
/// Geant4 simulation (mapped from the G4 creator process in the Geant4
/// integration's ParticleTrackingAction) so it can be stored on the particle
/// and written out.
enum class GenerationProcess : std::uint32_t {
  eUndefined = 0,
  eDecay = 1,
  ePhotonConversion = 2,
  eBremsstrahlung = 3,
  eNuclearInteraction = 4,
  // Geant4-fed only (see note above):
  eIonisation = 5,  // delta-ray / knock-on electron from ionisation
  eOther = 6,       // other material interaction (Compton, photoelectric, ...)
};

/// Print generation process type to output stream
/// @param os Output stream
/// @param processType Generation process type to print
/// @return Output stream
std::ostream &operator<<(std::ostream &os, GenerationProcess processType);

}  // namespace ActsFatras
