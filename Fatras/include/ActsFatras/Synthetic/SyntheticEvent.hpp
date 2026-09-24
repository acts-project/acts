// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"

#include <cstdint>
#include <vector>

namespace ActsFatras::Synthetic {

/// A generated particle, kept so that a caller can report the truth content of
/// an event without a truth-matching stage.
struct GeneratedParticle {
  /// Transverse momentum
  float pt{};
  /// Pseudorapidity
  float eta{};
  /// Azimuth of the momentum at the perigee, in [-pi, pi]
  float phi{};
  /// Transverse impact parameter with respect to the beam axis
  float d0{};
  /// Longitudinal position at the perigee
  float z0{};
  /// Charge, in units of the elementary one
  float charge{};
  /// Radius at which the particle was produced, zero for a primary
  float productionRadius{};
  /// Longitudinal position at which the particle was produced, zero for a
  /// primary
  float productionZ{};
  /// Index into `DetectorLayout::surfaces` of the surface it was produced on,
  /// `kNoSurface` for one from the luminous region or from a decay, both of
  /// which happen away from any surface
  std::uint32_t productionSurface{kNoSurface};
  /// Number of space points the particle left in the detector
  std::uint32_t numHits{};
  /// How many interactions back the luminous region is: zero for a particle
  /// out of it, one for what such a particle made, and so on. See
  /// `SecondaryConfig::maxGenerations`.
  std::uint8_t generation{};

  /// Whether the particle comes from the luminous region rather than from a
  /// material interaction or a decay.
  /// @return whether it is a primary
  bool primary() const { return generation == 0; }
};

/// A generated event.
///
/// Beyond the standard columns and the variances a triplet seeder needs, the
/// space points carry `layerId` into `DetectorLayout::layers`, `particleId`
/// into `Event::particles`, and the `clusterWidth` and `localPositionY` the
/// GBTS seeder expects.
struct Event {
  /// The space points
  Acts::SpacePointContainer spacePoints;
  /// The generating particles, indexed by the `particleId` column
  std::vector<GeneratedParticle> particles;
};

}  // namespace ActsFatras::Synthetic
