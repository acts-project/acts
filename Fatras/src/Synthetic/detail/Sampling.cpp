// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/detail/Sampling.hpp"

#include "ActsFatras/Synthetic/EventConfig.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <optional>
#include <random>

using namespace Acts::UnitLiterals;

namespace ActsFatras {

float Synthetic::detail::samplePt(const float u, const float minPt,
                                  const float ptScale, const float n) {
  // the root solve runs in double throughout: it is five fixed steps with no
  // convergence test, so it is given the margin rather than the parameters
  // being trusted to stay in the range that was tried
  const double scale = ptScale;
  // below two the spectrum does not normalise; clamped so that a fit walking
  // through the boundary gets a bad value and not a NaN
  const double exponent = std::max(static_cast<double>(n), 2.0001);
  const double a = 1.0 / (exponent - 2.0);
  const double b = 1.0 / (exponent - 1.0);
  const auto survival = [&](const double y) {
    const double v = 1.0 + y;
    return a * std::pow(v, 2.0 - exponent) - b * std::pow(v, 1.0 - exponent);
  };

  const double y0 = minPt / scale;
  const double total = survival(y0);
  const double target = total * (1.0 - static_cast<double>(u));
  constexpr int steps = 5;

  // the mode of the density, where the two conditionings meet
  const double yMode = 1.0 / (exponent - 1.0);
  if (target >= survival(yMode)) {
    // the density is at most `y`, which puts the root above this
    double lo = std::sqrt(std::max(y0 * y0 + 2.0 * (total - target), y0 * y0));
    double hi = yMode;
    double y = lo;
    for (int i = 0; i < steps; ++i) {
      const double s = survival(y);
      (s > target ? lo : hi) = y;
      const double w = 0.5 * y * y + (s - target) * std::pow(1.0 + y, exponent);
      const double next = std::sqrt(2.0 * std::max(w, 0.0));
      y = (lo <= next && next <= hi) ? next : 0.5 * (lo + hi);
    }
    return static_cast<float>(scale * y);
  }

  // `m` decreases with `y`, so `S = A m - B m^p` increases with `m`
  const double p = (1.0 - exponent) / (2.0 - exponent);
  double lo = 0.0;
  double hi = std::pow(1.0 + yMode, 2.0 - exponent);
  double m = std::min(target / a, hi);
  for (int i = 0; i < steps; ++i) {
    const double s = a * m - b * std::pow(m, p);
    (s > target ? hi : lo) = m;
    const double slope = a - b * p * std::pow(m, p - 1.0);
    const double next =
        slope > 0.0 ? m - (s - target) / slope : 0.5 * (lo + hi);
    m = (lo <= next && next <= hi) ? next : 0.5 * (lo + hi);
  }
  return static_cast<float>(scale *
                            (std::pow(m, 1.0 / (2.0 - exponent)) - 1.0));
}

float Synthetic::detail::sampleRapidity(std::mt19937& rng, const float lo,
                                        const float hi, const float edge,
                                        const float width) {
  if (width <= 0.f) {
    return sampleUniform(rng, lo, hi);
  }
  while (true) {
    const float rapidity = sampleUniform(rng, lo, hi);
    const float accept =
        1.f / (1.f + std::exp((std::abs(rapidity) - edge) / width));
    if (sampleUniform01(rng) < accept) {
      return rapidity;
    }
  }
}

std::optional<Synthetic::detail::SecondaryKinematics>
Synthetic::detail::sampleSecondary(
    std::mt19937& rng, const SecondarySamplingConfig& cfg, const float trackPhi,
    const float trackCotTheta, const float radialPhi,
    const float radialCotTheta, const float parentP,
    const float electronShare) {
  constexpr float pi = std::numbers::pi_v<float>;
  constexpr float pivot = 1_GeV;

  // The channel picks the axis: an electron follows a photon in from the beam
  // line, a nuclear product follows whatever crossed the surface.
  const bool electron = sampleUniform01(rng) < electronShare;
  const float axisPhi = electron ? radialPhi : trackPhi;
  const float axisCotTheta = electron ? radialCotTheta : trackCotTheta;

  // Longitudinal momentum and transverse kick, drawn independently; their
  // ratio is the opening angle.
  float alongP = 0.f;
  float acrossP = 0.f;
  if (electron || sampleUniform01(rng) >= cfg.evaporationFraction) {
    const float scale = electron ? cfg.electronScale : cfg.momentumScale;
    const float exponent =
        electron ? cfg.electronExponent : cfg.momentumExponent;
    const float width = electron ? cfg.electronSpread : cfg.momentumSpread;
    const float kt = electron ? cfg.electronKt : cfg.kt;
    alongP = sampleLogNormal(rng, scale * std::pow(parentP / pivot, exponent),
                             width);
    acrossP = sampleRayleigh(rng, kt);
  } else {
    // The struck nucleus breaking up: isotropic, so this is the only channel
    // that reaches the backward hemisphere.
    const float p = sampleRayleigh(rng, cfg.evaporationScale);
    const float cosOpen = sampleIsotropicCosine(rng);
    alongP = p * cosOpen;
    acrossP = p * std::sqrt(std::max(0.f, 1.f - cosOpen * cosOpen));
  }

  // A daughter cannot outrun its parent; scaling both parts keeps the angle.
  const float momentum = std::hypot(alongP, acrossP);
  if (momentum > parentP && momentum > 0.f) {
    const float shrink = parentP / momentum;
    alongP *= shrink;
    acrossP *= shrink;
  }
  const float bounded = std::min(momentum, parentP);
  const float cosOpening = bounded > 0.f ? alongP / bounded : 1.f;
  const float sinOpening = bounded > 0.f ? acrossP / bounded : 0.f;

  // Rotate the axis by that opening angle, about a direction drawn uniformly
  // around it, in the azimuthal and polar frame across the axis.
  const float axisSinTheta = 1.f / std::hypot(1.f, axisCotTheta);
  const float axisCosTheta = axisCotTheta * axisSinTheta;
  const float cosPhi = std::cos(axisPhi);
  const float sinPhi = std::sin(axisPhi);
  const float psi = sampleUniform(rng, -pi, pi);
  const float along = cosOpening;
  const float acrossPhi = sinOpening * std::cos(psi);
  const float acrossTheta = sinOpening * std::sin(psi);
  const float dx = along * axisSinTheta * cosPhi - acrossPhi * sinPhi -
                   acrossTheta * axisCosTheta * cosPhi;
  const float dy = along * axisSinTheta * sinPhi + acrossPhi * cosPhi -
                   acrossTheta * axisCosTheta * sinPhi;
  const float dz = along * axisCosTheta + acrossTheta * axisSinTheta;

  // the direction is a unit vector, so its transverse part is the daughter's
  // own sin(theta) and turns its momentum back into a pT
  const float sinThetaDaughter = std::max(std::hypot(dx, dy), 1e-6f);
  SecondaryKinematics out;
  out.phi = std::atan2(dy, dx);
  out.cotTheta = dz / sinThetaDaughter;
  out.pt = bounded * sinThetaDaughter;
  if (out.pt < cfg.minPt) {
    // too soft to leave the surface it was made on; the caller turns it into a
    // stub there rather than dropping it
    return std::nullopt;
  }
  return out;
}

}  // namespace ActsFatras
