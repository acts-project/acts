// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// The detector the synthetic event generator runs on: cylinders at a fixed
/// radius and discs at a fixed z, unresolved in azimuth. ACTS native units.

#include "Acts/Material/MaterialSlab.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <vector>

namespace ActsFatras::Synthetic {

/// The two shapes a layout is made of.
enum class SurfaceShape {
  /// A cylinder at a fixed radius, extending along z
  Cylinder,
  /// A disc at a fixed z, extending in r
  Disc,
};

/// Where along the beam axis a surface sits, valued as the sign of z.
enum class SurfaceSide : std::int32_t {
  /// The endcap at negative z
  Negative = -1,
  /// The barrel, which straddles the interaction point rather than sitting on
  /// one side of it
  Barrel = 0,
  /// The endcap at positive z
  Positive = +1,
};

/// A slab from the numbers a geometry reports for it.
/// @param x0 radiation length
/// @param l0 nuclear interaction length
/// @param ar relative atomic mass
/// @param z nuclear charge
/// @param molarDensity molar density
/// @param thickness thickness of the slab
/// @return the slab
Acts::MaterialSlab materialSlab(float x0, float l0, float ar, float z,
                                float molarDensity, float thickness);

/// What the bands of one surface share: what they are made of and the
/// thickness they are quoted at.
///
/// @note A band's density follows its own radiation length, so the molar
///       density and the radiation length enter only as their product.
struct BandComposition {
  /// Relative atomic mass
  float ar{};
  /// Nuclear charge
  float z{};
  /// Molar density times radiation length, zero for a surface of nothing
  float molarDensityX0{};
  /// Thickness every band is quoted at
  float thickness{};
};

/// What a surface is made of, band by band along the coordinate it extends in:
/// r on a disc, |z| on a cylinder.
struct SurfaceMaterial {
  SurfaceMaterial() = default;

  /// A surface of one material all over, i.e. a single unbounded band.
  /// @param uniform the slab it carries everywhere
  explicit SurfaceMaterial(const Acts::MaterialSlab& uniform)
      : bounds{std::numeric_limits<float>::infinity()}, bands{uniform} {}

  /// A surface whose bands each say what they are made of.
  /// @param bandBounds upper bound of each band
  /// @param bandMaterials what each band is made of, one per bound
  SurfaceMaterial(std::vector<float> bandBounds,
                  std::vector<Acts::MaterialSlab> bandMaterials);

  /// The same, from the two lengths per band the generated tables carry.
  /// @param composition what the bands share, see `BandComposition`
  /// @param bandBounds upper bound of each band
  /// @param radiationLengths radiation length of each band, zero for a band
  ///        that holds nothing
  /// @param nuclearLengths nuclear interaction length of each band
  SurfaceMaterial(const BandComposition& composition,
                  std::vector<float> bandBounds,
                  const std::vector<float>& radiationLengths,
                  const std::vector<float>& nuclearLengths);

  /// Upper bound of each band, increasing; band `k` starts at `bounds[k - 1]`
  /// or at zero. Past the last there is nothing.
  std::vector<float> bounds;
  /// What each band is made of, one per entry of `bounds`
  std::vector<Acts::MaterialSlab> bands;

  /// Scale how much material the surface carries.
  /// @param scale what to multiply the thickness by
  void scaleThickness(const float scale) {
    for (Acts::MaterialSlab& band : bands) {
      band.scaleThickness(scale);
    }
  }

  /// What a crossing at a position along the surface meets.
  /// @param along r on a disc, |z| on a cylinder
  /// @return the material there, vacuum past the end of it
  const Acts::MaterialSlab& at(const float along) const {
    static const Acts::MaterialSlab nothing = Acts::MaterialSlab::Vacuum(0.f);
    if (bounds.empty() || along >= bounds.back()) {
      return nothing;
    }
    const auto band = std::ranges::upper_bound(bounds, along);
    return bands[band - bounds.begin()];
  }

  /// What the surface is made of on average; a crossing meets `at()`.
  /// @return the bands combined, taken back to one band's worth of thickness
  Acts::MaterialSlab average() const {
    if (bands.empty()) {
      return Acts::MaterialSlab::Vacuum(0.f);
    }
    // the bands sit side by side rather than stacked
    Acts::MaterialSlab combined = Acts::MaterialSlab::combineLayers(bands);
    combined.scaleThickness(1.f / static_cast<float>(bands.size()));
    return combined;
  }
};

/// One logical layer, the granularity a seeder reasons about: a cylinder split
/// into eta modules along z, a disc into rings in r.
struct DetectorLayer {
  /// Cylinder or disc
  SurfaceShape shape{};
  /// Radius for a cylinder, signed z position for a disc
  float refCoord{};
  /// Lower bound along the extended coordinate: z for a cylinder, r for a disc
  float minBound{};
  /// Upper bound along the extended coordinate: z for a cylinder, r for a disc
  float maxBound{};
  /// The endcap the layer belongs to, or the barrel
  SurfaceSide side{SurfaceSide::Barrel};
  /// Barrel layer index, or disc index counting outwards from the interaction
  /// point within one endcap
  std::uint32_t layer{};
  /// Eta module index for a cylinder, ring index for a disc
  std::uint32_t moduleIndex{};
};

/// The radial bounds of one ring of an endcap disc.
struct RingBounds {
  /// Inner radius
  float rMin{};
  /// Outer radius
  float rMax{};
};

/// One endcap disc: where it sits along the beam axis and the rings it carries.
struct DiscDescription {
  /// Absolute longitudinal position
  float absZ{};
  /// Radial bounds of the rings, increasing. Gaps are allowed and leave no
  /// space point.
  std::vector<RingBounds> rings;
  /// What a crossing meets, resolved in r.
  ///
  /// @note Last, so that tables initialising the position and the rings
  ///       positionally keep working.
  SurfaceMaterial material;
};

/// A passive surface the geometry carries away from any sensitive layer: a
/// support, a service disc, a cable run. Material without readout.
struct PassiveSurfaceDescription {
  /// Cylinder or disc
  SurfaceShape shape{};
  /// Radius for a cylinder, absolute z for a disc. Mirrored onto both sides.
  float refCoord{};
  /// Where its material starts along the coordinate it extends in: |z| for a
  /// cylinder, r for a disc
  float minBound{};
  /// Where it ends; zero leaves it unbounded, as a service surface spanning the
  /// detector is
  float maxBound{};
  /// What a crossing meets, read off the geometry
  SurfaceMaterial material;
};

/// A barrel of cylinders plus two mirrored endcaps of discs, optionally behind
/// a passive beam pipe. Nothing about it is uniform: every cylinder has its own
/// half-length and every disc its own z and rings.
struct BarrelEndcapDescription {
  /// Radius of the passive beam pipe, if there is one
  std::optional<float> beamPipeRadius;
  /// What a crossing of the beam pipe meets, read off the geometry
  SurfaceMaterial beamPipeMaterial;
  /// Passive surfaces the geometry carries away from any sensitive layer
  std::vector<PassiveSurfaceDescription> passiveSurfaces;
  /// Radii of the barrel cylinders, increasing
  std::vector<float> barrelRadii;
  /// Half-length along z of each barrel cylinder, one per entry of
  /// `barrelRadii`
  std::vector<float> barrelHalfLengthsZ;
  /// What a crossing of each barrel cylinder meets, read off the geometry and
  /// resolved in |z|. Empty leaves them all at vacuum.
  std::vector<SurfaceMaterial> barrelMaterials;
  /// Chance per barrel cylinder that a crossing also leaves a space point on
  /// the module next to it; see `DetectorSurface::overlapProbability`.
  std::vector<float> barrelOverlapProbabilities;
  /// The same for every endcap disc, one number because the measurement is flat
  /// across them
  float discOverlapProbability{0.f};
  /// How far an overlapping module sits from the one that was crossed; see
  /// `DetectorSurface::overlapOffset`
  float barrelOverlapOffset{0.f};
  /// @copydoc barrelOverlapOffset
  float discOverlapOffset{0.f};
  /// Extent of the tracker these surfaces sit inside; see
  /// `DetectorLayout::escapeRadius`
  float escapeRadius{1100.f};
  /// @copydoc escapeRadius
  float escapeHalfZ{3000.f};
  /// Eta modules each barrel cylinder is split into
  std::uint32_t barrelModules{1};
  /// The endcap discs, increasing in absolute z. Each is mirrored onto both
  /// sides.
  std::vector<DiscDescription> discs;
};

/// Stands for "no surface" where an index into `DetectorLayout::surfaces` is
/// optional.
constexpr std::uint32_t kNoSurface = std::numeric_limits<std::uint32_t>::max();

/// A physical detector element and the logical layers it is split into.
/// Propagation crosses surfaces and the position picks the layer. One without
/// layers is passive: material without readout.
struct DetectorSurface {
  /// Cylinder or disc
  SurfaceShape shape{};
  /// Radius for a cylinder, signed z position for a disc
  float refCoord{};
  /// Indices into `DetectorLayout::layers` of the logical layers of this
  /// surface, empty for a passive surface
  std::vector<std::uint32_t> layers;
  /// What a crossing meets, and the whole of the material model: `x/X0`
  /// scatters and makes electrons, `x/L0` makes nuclear products, and the slab
  /// takes the energy. Nothing is fitted on top.
  SurfaceMaterial material;
  /// Chance that a crossing also leaves a space point on the module next to it
  /// in phi, the modules of a stave or a ring overlapping there.
  float overlapProbability{0.f};
  /// How far that neighbour sits along the surface normal, adjacent modules
  /// alternating in depth rather than sitting side by side.
  float overlapOffset{0.f};
  /// Where the surface holds anything at all, in the coordinate it extends in:
  /// |z| on a cylinder, r on a disc. Outside it a crossing meets neither a
  /// layer nor material, which is what most crossings of its plane do.
  ///
  /// A passive surface is bounded by what it was declared with, a sensitive one
  /// by what its layers and its bands of material span together;
  /// `updateSurfaceExtents` derives the latter and leaves it wide open until
  /// then, which is correct and slower.
  float minBound{0.f};
  /// @copydoc minBound
  float maxBound{std::numeric_limits<float>::infinity()};
};

/// The full synthetic detector description.
struct DetectorLayout {
  /// All logical layers, in the index space used by the `layerId` column
  std::vector<DetectorLayer> layers;
  /// The physical surfaces the layers are grouped into
  std::vector<DetectorSurface> surfaces;

  /// @name Containment
  ///
  /// Where the *enclosing* tracker ends, not this layout: a track leaving the
  /// pixels curls back through the strips in the same field.
  ///
  /// @{

  /// Radius past which a track has left the tracker for good
  float escapeRadius{1100.f};
  /// Longitudinal half length of the same volume
  float escapeHalfZ{3000.f};

  /// @}
};

/// Rings tiling `[rMin, rMax]` without gaps.
/// @param rMin the inner radius of the first ring
/// @param rMax the outer radius of the last ring
/// @param numRings the number of rings
/// @return the ring bounds
std::vector<RingBounds> uniformRings(float rMin, float rMax,
                                     std::uint32_t numRings);

/// Split each ring into `parts` of equal radial width, keeping the structure.
/// @param rings the rings to subdivide
/// @param parts how many parts to split each ring into; one returns them
///        unchanged
/// @return the subdivided rings
std::vector<RingBounds> subdivideRings(std::span<const RingBounds> rings,
                                       std::uint32_t parts);

/// Assembles a `DetectorLayout` surface by surface, keeping each surface and
/// its layers cross-linked. Layers are numbered in call order.
class DetectorLayoutBuilder {
 public:
  /// Add a passive cylinder, i.e. material without readout.
  /// @param radius the radius
  /// @param maxAbsZ where it ends along z
  /// @param minAbsZ where it starts, for a band of a longer cylinder
  /// @return this builder, for chaining
  DetectorLayoutBuilder& addPassiveCylinder(
      float radius, float maxAbsZ = std::numeric_limits<float>::infinity(),
      float minAbsZ = 0.f);

  /// Add a passive disc on one side.
  /// @param side the endcap the disc sits on, which cannot be the barrel
  /// @param absZ the absolute longitudinal position
  /// @param rMin the inner radius of the material
  /// @param rMax the outer radius of the material
  /// @return this builder, for chaining
  DetectorLayoutBuilder& addPassiveDisc(
      SurfaceSide side, float absZ, float rMin = 0.f,
      float rMax = std::numeric_limits<float>::infinity());

  /// Add a sensitive cylinder, split into eta modules along z.
  /// @param radius the radius
  /// @param halfLengthZ the half-length along z
  /// @param numModules the number of eta modules to split it into
  /// @return this builder, for chaining
  DetectorLayoutBuilder& addCylinder(float radius, float halfLengthZ,
                                     std::uint32_t numModules);

  /// Add a sensitive disc on one side, with the rings as given so that they
  /// may leave radial gaps.
  ///
  /// @param side the endcap the disc sits on, which cannot be the barrel
  /// @param absZ the absolute longitudinal position
  /// @param rings the radial bounds of the rings, increasing and non-overlapping
  /// @return this builder, for chaining
  DetectorLayoutBuilder& addDisc(SurfaceSide side, float absZ,
                                 std::span<const RingBounds> rings);

  /// Add a sensitive disc on one side, split into rings tiling `[rMin, rMax]`.
  /// @param side the endcap the disc sits on, which cannot be the barrel
  /// @param absZ the absolute longitudinal position
  /// @param rMin the inner radius
  /// @param rMax the outer radius
  /// @param numRings the number of rings to split it into
  /// @return this builder, for chaining
  DetectorLayoutBuilder& addDisc(SurfaceSide side, float absZ, float rMin,
                                 float rMax, std::uint32_t numRings);

  /// Take the layout that has been assembled, leaving the builder empty.
  /// @return the layout
  DetectorLayout build();

 private:
  DetectorLayout m_layout;
  /// Number of sensitive cylinders added so far
  std::uint32_t m_numCylinders{};
  /// Number of sensitive discs added so far, per side
  std::array<std::uint32_t, 2> m_numDiscs{};

  /// @param side the endcap, which cannot be the barrel
  /// @return a reference to the disc counter of that side
  std::uint32_t& discCounter(SurfaceSide side);
};

/// Fill in `DetectorSurface::minBound` and `DetectorSurface::maxBound` of every
/// sensitive surface of a layout, a passive one carrying its own. Call it once
/// the layers and the material are in place, and again whenever either changes.
///
/// @param layout the layout to update in place
void updateSurfaceExtents(DetectorLayout& layout);

/// Build a barrel-endcap layout.
/// @param description the detector to build
/// @return the layout
DetectorLayout makeLayout(const BarrelEndcapDescription& description);

/// A coarse stand-in for the ACTS Generic detector's pixel system.
/// @return the description of the generic pixel layout
BarrelEndcapDescription genericDetectorPixelDescription();

/// Build a stand-in for the ACTS Generic detector's pixel system.
/// @return the generic pixel layout
DetectorLayout makeGenericDetectorPixelLayout();

}  // namespace ActsFatras::Synthetic
