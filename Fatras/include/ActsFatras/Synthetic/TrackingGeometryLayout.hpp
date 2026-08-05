// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Reduce a real detector to the RZ-symmetric skeleton the synthetic event
/// generator runs on, by collapsing each layer of an `Acts::TrackingGeometry`
/// onto a cylinder or a disc. This is the general way to get a layout.

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"

#include <cstdint>
#include <functional>
#include <optional>

namespace Acts {
class Surface;
class TrackingGeometry;
}  // namespace Acts

namespace ActsFatras::Synthetic {

/// Steering for `makeLayoutFromTrackingGeometry`.
struct TrackingGeometryLayoutOptions {
  /// Keep only surfaces for which this returns true; empty keeps everything.
  /// @note This applies to the material surfaces too, so a volume-based
  ///       selector also decides which passive material the layout gets. A
  ///       pixel-only selector excludes the beam pipe, which lives in a volume
  ///       of its own.
  std::function<bool(const Acts::Surface&)> surfaceSelector;

  /// Also add every non-sensitive surface that carries material as a passive
  /// surface. Off by default: a real geometry carries far more material
  /// surfaces than sensitive ones, each of them another crossing at which a
  /// secondary can be produced, so turning this on means the rates in
  /// `SecondaryConfig` have to be recalibrated.
  bool includeMaterialSurfaces = false;

  /// Fill `DetectorSurface::material` with the material the geometry carries.
  /// A mapped geometry puts a layer's material on its representing and approach
  /// surfaces and none on its sensors, so this is taken per layer and shared
  /// out over the surfaces that layer reduces to. A geometry without material
  /// is unaffected.
  bool materialFromGeometry = true;

  /// Bands each surface's material is sampled into along the coordinate it
  /// extends in, before the uniform ones among them are merged back together.
  /// One collapses a surface to its average, which puts a service cylinder's
  /// endcap material into the barrel it also passes through.
  std::uint32_t materialBands = 50;

  /// Relative difference below which two neighbouring bands are merged into
  /// one. Zero keeps every band, which makes for tables no one can read;
  /// large enough and a surface comes out uniform again.
  double materialBandTolerance = 0.2;

  /// Insert a passive cylinder at this radius, standing in for the beam pipe:
  /// the only material in front of the innermost layer, and hence the only
  /// source of secondaries there, but usually excluded by `surfaceSelector`.
  std::optional<float> passiveBeamPipeRadius;

  /// Number of logical layers each detected surface is split into, along its
  /// extended coordinate. For a disc this subdivides each of its rings rather
  /// than the disc as a whole.
  std::uint32_t modulesPerSurface = 1;

  /// Surfaces of one layer within this distance in z are one disc. Small,
  /// `maxRingOverlap` below being what guards against splitting too eagerly.
  double discZTolerance = 2 * Acts::UnitConstants::mm;

  /// Fraction of the narrower band by which two z planes of a layer may overlap
  /// radially and still count as separate rings; above it they are merged back
  /// into one disc.
  /// Two planes covering the same radii are one ring whose modules alternate
  /// in z, and splitting those would give a track two space points where the
  /// detector gives it one; genuinely different rings barely overlap in r.
  ///
  /// @note Every pair of z planes has to be tested, not only neighbours in z:
  ///       the generic detector's strip endcap interleaves its rings, so the
  ///       two halves of one ring can have another ring between them.
  double maxRingOverlap = 0.25;

  /// Radial gap between two rings of a disc below which they count as one.
  /// Modules of one ring overlap slightly in r and neighbouring rings are tens
  /// of millimetres apart, so anything in between is a matter of taste.
  double ringRTolerance = 5 * Acts::UnitConstants::mm;

  /// Radial extent below which a group of surfaces counts as a cylinder rather
  /// than a disc. The classification cannot use `Acts::Surface::type`, the
  /// sensitive surfaces of a real tracker being planes however they are
  /// arranged; what separates a barrel layer is that it is thin in r.
  double cylinderRTolerance = 30 * Acts::UnitConstants::mm;
};

/// Collapse a tracking geometry onto a synthetic detector layout.
///
/// Sensitive surfaces are grouped by the volume and layer of their geometry
/// identifier. A group thin in r becomes one cylinder; a group thin in z
/// becomes one disc per z plane its surfaces fall into, carrying one ring per
/// radial band of them.
///
/// @param trackingGeometry the detector to reduce
/// @param gctx the geometry context to evaluate surface positions in
/// @param options steering, see `TrackingGeometryLayoutOptions`
/// @return the layout
DetectorLayout makeLayoutFromTrackingGeometry(
    const Acts::TrackingGeometry& trackingGeometry,
    const Acts::GeometryContext& gctx,
    const TrackingGeometryLayoutOptions& options = {});

}  // namespace ActsFatras::Synthetic
