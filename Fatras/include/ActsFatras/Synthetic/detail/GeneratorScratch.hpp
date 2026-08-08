// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/Synthetic/SyntheticEvent.hpp"
#include "ActsFatras/Synthetic/detail/Helix.hpp"
#include "ActsFatras/Synthetic/detail/Propagation.hpp"

#include <cmath>
#include <cstdint>
#include <numbers>
#include <vector>

namespace ActsFatras::Synthetic::detail {

/// One space point, in the cylindrical coordinates it is measured in.
struct SmearedHit {
  float r{};
  float phi{};
  float z{};
  /// Index into `DetectorLayout::layers`
  std::uint32_t layer{};
  /// Index into the particles the tracks are propagated alongside
  std::uint32_t particle{};
};

/// A secondary too soft to leave the surface it was made on: not propagated,
/// only measured where it was born.
struct PendingStub {
  /// Radius of the production point
  float r{};
  /// Azimuth
  float phi{};
  /// Longitudinal position
  float z{};
  /// Transverse momentum, below `SecondarySamplingConfig::minPt`
  float pt{};
  /// Charge, in units of the elementary one
  float charge{};
  /// Index into `DetectorLayout::layers` of the layer it was made on
  std::uint32_t layer{};
  /// Index into `DetectorLayout::surfaces` of the surface it was made on
  std::uint32_t surface{};
  /// One more than the generation of whatever made it
  std::uint8_t generation{};
};

/// The particle a track is recorded as, in the coordinates a consumer wants.
/// Where it came from is the caller's to fill in.
inline GeneratedParticle makeParticle(const Helix& helix, const float pt) {
  GeneratedParticle particle;
  particle.pt = pt;
  particle.eta = std::asinh(helix.cotTheta);
  // the helix azimuth is not wrapped while it is propagated
  particle.phi = std::remainder(helix.phi0, 2.f * std::numbers::pi_v<float>);
  particle.d0 = helix.d0;
  particle.z0 = helix.z0;
  particle.charge = helix.charge;
  return particle;
}

/// Per-event working storage, reused across events rather than reallocated.
struct GeneratorScratch {
  /// Tracks to propagate, growing as secondaries are queued. A track is its
  /// helix and nothing else: what it is recorded as sits at the same index of
  /// the particles it is propagated alongside, and whatever appends to one
  /// appends to the other.
  std::vector<Helix> tracks;
  /// Hits collected so far this event
  std::vector<SmearedHit> hits;
  /// Crossings of the track currently being propagated
  std::vector<SurfaceCrossing> crossings;
  /// Secondaries too soft to be propagated, resolved into space points once the
  /// particle indices are known
  std::vector<PendingStub> stubs;
};

}  // namespace ActsFatras::Synthetic::detail
