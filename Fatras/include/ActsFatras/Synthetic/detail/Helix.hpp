// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Closed-form helix propagation in a solenoid field along z.

#include "Acts/Utilities/detail/periodic.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <optional>

namespace ActsFatras::Synthetic::detail {

/// A helix in a solenoid field along z, at its perigee to the beam axis.
struct Helix {
  /// Azimuth of the momentum at the perigee
  float phi0{};
  /// Longitudinal position at the perigee
  float z0{};
  /// Ratio of the longitudinal to the transverse momentum
  float cotTheta{};
  /// Radius of curvature
  float radius{};
  /// Charge sign, +1 or -1
  float charge{};
  /// Distance of the centre of the circle from the beam axis
  float rCentre{};
  /// Signed transverse impact parameter with respect to the beam axis. Stored
  /// rather than taken from `rCentre`, a difference that cancels in float once
  /// it is much smaller than `radius`.
  float d0{};
  /// Turning angle the track starts at; below it the helix does not exist
  float minGamma{};
};

/// The radius of curvature of a unit-charge track.
inline float helixRadius(const float pt, const float bFieldZ) {
  return pt / bFieldZ;
}

/// Build a helix from perigee parameters.
inline Helix makeHelix(const float d0, const float phi0, const float z0,
                       const float cotTheta, const float pt, const float charge,
                       const float bFieldZ) {
  Helix helix;
  helix.phi0 = phi0;
  helix.z0 = z0;
  helix.cotTheta = cotTheta;
  helix.radius = helixRadius(pt, bFieldZ);
  helix.charge = charge;
  helix.d0 = d0;
  helix.rCentre = helix.radius - charge * d0;
  return helix;
}

/// The smallest turning angle at which a helix reaches a radius, ignoring
/// `Helix::minGamma`. Nothing if it never reaches it.
inline std::optional<float> helixBaseGammaAtRadius(const Helix& helix,
                                                   const float r) {
  // The law of cosines solved for `sin^2(gamma / 2)` rather than for
  // `cos(gamma)`, which rounds to one in float at the innermost radii above
  // about 30 GeV in a 2 T field.
  const float denominator = 4.f * helix.rCentre * helix.radius;
  if (denominator <= 0.f) {
    return std::nullopt;
  }
  const float halfSinSquared = (r * r - helix.d0 * helix.d0) / denominator;
  if (halfSinSquared < 0.f) {
    // the circle never comes as close to the beam axis as this radius
    return std::nullopt;
  }
  if (halfSinSquared > 1.f) {
    // the particle curls up before this radius, or never comes back to it
    return std::nullopt;
  }
  return 2.f * std::asin(std::sqrt(halfSinSquared));
}

/// Carry a turning angle forward by whole turns into `[from, from + 2 * pi)`.
inline float helixGammaAfter(const float gamma, const float from) {
  return Acts::detail::wrap_periodic(gamma, from,
                                     2.f * std::numbers::pi_v<float>);
}

/// Whether a turning angle falls on the half turn heading back inwards, where
/// the azimuth is the mirror of the outgoing one.
inline bool helixIsInbound(const float gamma) {
  return Acts::detail::radian_pos(gamma) > std::numbers::pi_v<float>;
}

/// Mirror an azimuth across the line through the beam axis and the centre of
/// the circle, i.e. turn an outgoing crossing into the inbound one at the same
/// radius.
inline float helixPhiMirrored(const Helix& helix, const float phi) {
  constexpr float pi = std::numbers::pi_v<float>;
  return 2.f * helix.phi0 - helix.charge * pi - phi;
}

/// The azimuth of a helix at a radius.
inline float helixPhiAtRadius(const Helix& helix, const float r) {
  const float sinPhi =
      (r * r - helix.charge * helix.d0 * (helix.rCentre + helix.radius)) /
      (2.f * r * helix.rCentre);
  return helix.phi0 - helix.charge * std::asin(std::clamp(sinPhi, -1.f, 1.f));
}

/// The radius a helix has reached after a turning angle, the inverse of
/// `helixBaseGammaAtRadius`.
inline float helixRadiusAtGamma(const Helix& helix, const float gamma) {
  const float halfSin = std::sin(0.5f * gamma);
  const float squared = helix.d0 * helix.d0 +
                        4.f * helix.rCentre * helix.radius * halfSin * halfSin;
  return std::sqrt(std::max(0.f, squared));
}

/// Rebuild a helix from the production point `(r, phi, z)` and the momentum
/// there, `direction` being its azimuth. `Helix::minGamma` restricts it to the
/// part beyond that point.
inline Helix makeHelixFromPoint(const float r, const float phi, const float z,
                                const float direction, const float cotTheta,
                                const float pt, const float charge,
                                const float bFieldZ) {
  constexpr float pi = std::numbers::pi_v<float>;

  Helix helix;
  helix.cotTheta = cotTheta;
  helix.radius = helixRadius(pt, bFieldZ);
  helix.charge = charge;

  const float deltaPhi = direction - phi;
  const float centreX =
      r * std::cos(phi) + charge * helix.radius * std::sin(direction);
  const float centreY =
      r * std::sin(phi) - charge * helix.radius * std::cos(direction);
  helix.rCentre = std::hypot(centreX, centreY);
  // `radius - rCentre` cancels for a track pointing near the beam axis; the
  // same difference taken over their sum does not.
  helix.d0 = -(charge * r * r + 2.f * helix.radius * r * std::sin(deltaPhi)) /
             (helix.radius + helix.rCentre);
  helix.phi0 = std::atan2(centreY, centreX) + charge * pi / 2.f;

  // Which leg the point is on: the radius cannot say, the momentum can.
  float gamma = helixBaseGammaAtRadius(helix, r).value_or(0.f);
  if (std::cos(deltaPhi) < 0.f) {
    gamma = 2.f * pi - gamma;
  }
  helix.z0 = z - cotTheta * helix.radius * gamma;
  helix.minGamma = gamma;
  return helix;
}

/// Turning angle past which a track has left the tracker for good, radially or
/// along z, capped at `maxGamma`. One that never leaves is left there.
inline float helixEscapeGamma(const Helix& helix, const float escapeRadius,
                              const float escapeHalfZ, const float maxGamma) {
  float gamma = maxGamma;
  if (helix.rCentre + helix.radius > escapeRadius) {
    if (const std::optional<float> exit =
            helixBaseGammaAtRadius(helix, escapeRadius);
        exit.has_value()) {
      // the first crossing after the track was produced, not the first of all
      gamma = std::min(gamma, helixGammaAfter(*exit, helix.minGamma));
    }
  }
  // z grows at a fixed rate with the turning angle, in the sign of the drift
  const float drift = helix.cotTheta * helix.radius;
  if (drift != 0.f) {
    const float room = escapeHalfZ - (drift > 0.f ? helix.z0 : -helix.z0);
    gamma = std::min(gamma, std::max(room / std::abs(drift), 0.f));
  }
  return gamma;
}

}  // namespace ActsFatras::Synthetic::detail
