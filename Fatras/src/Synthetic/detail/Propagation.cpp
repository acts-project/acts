// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/detail/Propagation.hpp"

#include <cmath>
#include <numbers>
#include <optional>
#include <utility>

namespace ActsFatras::Synthetic::detail {

namespace {

/// Whether a crossing landed on the surface at all, `along` being `|z|` on a
/// cylinder and `r` on a disc. The precondition of `acceptCrossing`, tested
/// first because most crossings of a surface's plane are off its end, where
/// there is neither a layer nor material.
bool onSurface(const DetectorSurface& surface, const float along) {
  return surface.minBound <= along && along <= surface.maxBound;
}

/// Fill in what a crossing needs beyond its position, and keep it if it either
/// lands on a layer or meets material. Called once the crossing is known to be
/// `onSurface`, which is what `along` was tested against.
bool acceptCrossing(const DetectorLayout& layout,
                    const DetectorSurface& surface, const float along,
                    SurfaceCrossing& crossing) {
  crossing.material = &surface.material.at(along);
  if (surface.layers.empty()) {
    // A passive surface is bounded by its own extent alone, having no layers;
    // without one a curling track keeps crossing the beam pipe for ever.
    return true;
  }
  if (const std::optional<std::uint32_t> layer =
          findLayer(layout, surface, crossing.r, crossing.z);
      layer.has_value()) {
    crossing.layer = *layer;
    return true;
  }
  // Material off every layer is still material -- the gap between two rings is
  // the support carrying them. Only the hit needs a layer.
  return !crossing.material->isVacuum();
}

}  // namespace

}  // namespace ActsFatras::Synthetic::detail

namespace ActsFatras {

std::optional<std::uint32_t> Synthetic::detail::findLayer(
    const DetectorLayout& layout, const DetectorSurface& surface, const float r,
    const float z) {
  const float selector = surface.shape == SurfaceShape::Cylinder ? z : r;
  for (const std::uint32_t l : surface.layers) {
    const DetectorLayer& layer = layout.layers[l];
    if (selector >= layer.minBound && selector < layer.maxBound) {
      return l;
    }
  }
  return std::nullopt;
}

void Synthetic::detail::propagateHelix(const DetectorLayout& layout,
                                       const Helix& helix,
                                       std::vector<SurfaceCrossing>& crossings,
                                       const std::uint32_t productionSurface,
                                       const float maxGamma) {
  constexpr float pi = std::numbers::pi_v<float>;
  // 1 / sin(theta), the transverse path length of a unit step along the track
  const float invSinTheta = std::hypot(1.f, helix.cotTheta);
  // How far along z the track gets: z grows with the turning angle at a fixed
  // rate, so its two ends bound every crossing between them and a cylinder the
  // track never reaches is skipped before it is intersected.
  const float drift = helix.cotTheta * helix.radius;
  const float zFirst = helix.z0 + drift * helix.minGamma;
  const float zLast = helix.z0 + drift * maxGamma;
  const float minAbsZ = (zFirst < 0.f) == (zLast < 0.f)
                            ? std::min(std::abs(zFirst), std::abs(zLast))
                            : 0.f;
  const float maxAbsZ = std::max(std::abs(zFirst), std::abs(zLast));

  const std::uint32_t numSurfaces =
      static_cast<std::uint32_t>(layout.surfaces.size());
  // Filled surface by surface rather than declared inside the loop, which would
  // clear every field of it for each of them.
  SurfaceCrossing crossing;
  for (std::uint32_t s = 0; s < numSurfaces; ++s) {
    const DetectorSurface& surface = layout.surfaces[s];

    if (surface.shape == SurfaceShape::Cylinder) {
      if (surface.maxBound < minAbsZ || surface.minBound > maxAbsZ) {
        // the track is past either end of it for the whole of its path
        continue;
      }
      const std::optional<float> base =
          helixBaseGammaAtRadius(helix, surface.refCoord);
      if (!base.has_value()) {
        continue;
      }
      // A circle crosses a radius twice per turn, outwards at `base` and back
      // inwards at `2 * pi - base`, the two azimuths related by a reflection
      // across the line through the centre. Drawn on the first crossing the
      // surface keeps, since most tracks pass its radius beyond either end.
      std::optional<std::pair<float, float>> phiOutIn;
      for (std::uint32_t leg = 0; leg < 2; ++leg) {
        if (leg == 1 && *base >= pi) {
          // the radius is the turning point of the helix, so the two legs are
          // the same crossing
          continue;
        }
        // this leg's first turn after the production point, then one a turn
        float first = helixGammaAfter(leg == 0 ? *base : 2.f * pi - *base,
                                      helix.minGamma);
        // A departure is not a crossing, and it sits exactly on `minGamma`:
        // the same intersection of the same helix, so they agree to the bit.
        if (s == productionSurface && first == helix.minGamma) {
          first += 2.f * pi;
        }
        for (float gamma = first; gamma <= maxGamma; gamma += 2.f * pi) {
          crossing.surface = s;
          crossing.layer = kNoLayer;
          crossing.r = surface.refCoord;
          crossing.gamma = gamma;
          crossing.z = helix.z0 + drift * gamma;
          const float absZ = std::abs(crossing.z);
          if (!onSurface(surface, absZ) ||
              !acceptCrossing(layout, surface, absZ, crossing)) {
            continue;
          }
          if (!phiOutIn.has_value()) {
            const float phiOut = helixPhiAtRadius(helix, surface.refCoord);
            phiOutIn = {phiOut, helixPhiMirrored(helix, phiOut)};
          }
          crossing.phi =
              helixIsInbound(gamma) ? phiOutIn->second : phiOutIn->first;
          // the angle between the momentum and the radial normal, in the
          // transverse plane; the polar part of the incidence is `invSinTheta`
          const float cosIncidence =
              std::cos(helix.phi0 - helix.charge * gamma - crossing.phi);
          crossing.pathLength =
              invSinTheta / std::max(std::abs(cosIncidence), 1e-3f);
          crossings.push_back(crossing);
        }
      }
      continue;
    }

    if (s == productionSurface || helix.cotTheta == 0.f) {
      continue;
    }
    // A disc is crossed once however far the track curls: z grows with the
    // turning angle at a fixed rate, so there is no return branch.
    crossing.surface = s;
    crossing.layer = kNoLayer;
    crossing.z = surface.refCoord;
    crossing.gamma = (crossing.z - helix.z0) / drift;
    if (crossing.gamma <= helix.minGamma || crossing.gamma >= maxGamma) {
      // the disc is behind the particle, or it curls up before reaching it
      continue;
    }
    crossing.r = helixRadiusAtGamma(helix, crossing.gamma);
    if (crossing.r <= 0.f) {
      continue;
    }
    if (!onSurface(surface, crossing.r) ||
        !acceptCrossing(layout, surface, crossing.r, crossing)) {
      continue;
    }
    crossing.phi = helixPhiAtRadius(helix, crossing.r);
    if (helixIsInbound(crossing.gamma)) {
      // on this half-turn the track is on its way back in, and its azimuth is
      // mirrored just as a cylinder crossing's is
      crossing.phi = helixPhiMirrored(helix, crossing.phi);
    }
    // the normal of a disc is the beam axis, so only the polar angle enters
    crossing.pathLength =
        invSinTheta / std::max(std::abs(helix.cotTheta), 1e-3f);
    crossings.push_back(crossing);
  }
}

}  // namespace ActsFatras
