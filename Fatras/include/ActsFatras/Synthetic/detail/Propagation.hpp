// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/detail/Helix.hpp"

#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <vector>

namespace ActsFatras::Synthetic::detail {

/// Stands for "no layer" where a layer index is optional.
constexpr std::uint32_t kNoLayer = std::numeric_limits<std::uint32_t>::max();

/// One crossing of a detector surface by a helix.
struct SurfaceCrossing {
  /// Radius
  float r{};
  /// Azimuth
  float phi{};
  /// Longitudinal position
  float z{};
  /// Turning angle, which gives the momentum direction there
  float gamma{};
  /// Path through the surface in units of its normal thickness, i.e. the
  /// reciprocal of the cosine of the incidence. Unbounded, so it is clamped
  /// where it is used.
  float pathLength{1.f};
  /// The band of the surface the crossing landed in, owned by the layout
  const Acts::MaterialSlab* material{};
  /// Index into `DetectorLayout::surfaces` of the crossed surface
  std::uint32_t surface{kNoSurface};
  /// Index into `DetectorLayout::layers` of the layer the crossing landed on;
  /// `kNoLayer` where it missed every one and is kept for its material alone
  std::uint32_t layer{kNoLayer};

  /// Whether the crossing measures the particle, which needs a layer to do
  constexpr bool sensitive() const { return layer != kNoLayer; }
};

/// The index into `layout.layers` of the layer of `surface` a position falls
/// into, or nothing if it is outside the bounds of every one of them.
std::optional<std::uint32_t> findLayer(const DetectorLayout& layout,
                                       const DetectorSurface& surface, float r,
                                       float z);

/// Propagate a helix through the detector, appending its crossings to
/// `crossings` in surface order. `productionSurface` indexes the surface the
/// track started on and therefore does not cross. Past half a turn of
/// `maxGamma` a track curls back through the surfaces it has already crossed.
void propagateHelix(const DetectorLayout& layout, const Helix& helix,
                    std::vector<SurfaceCrossing>& crossings,
                    std::uint32_t productionSurface = kNoSurface,
                    float maxGamma = std::numbers::pi_v<float>);

}  // namespace ActsFatras::Synthetic::detail
