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
#include <string>
#include <utility>
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

/// Which sides of the interaction point a described disc is built on. Positions
/// are quoted as absolute z whichever it is, so a mirrored endcap is written
/// once.
enum class EndcapPlacement {
  /// Negative z only
  Negative,
  /// Positive z only
  Positive,
  /// Both, from the one description
  Mirrored,
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

/// What the bands of one surface share, where they do: what they are made of
/// and the thickness they are quoted at. A convenience for stating a surface
/// whose composition does not vary along itself, which is most of them — the
/// bands themselves are free to differ in every respect.
///
/// @note A band's density follows its own radiation length, so the molar
///       density and the radiation length enter only as their product. `ar` and
///       `z` reach nothing but the energy loss, and vary by a tenth across the
///       materials a tracker is built of, which is why sharing them costs
///       little.
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
///
/// @note Unresolved in azimuth, as the whole layout is. A surface that is
///       copper on one side and carbon on the other is the average of the two.
struct SurfaceMaterial {
  SurfaceMaterial() = default;

  /// A surface of one material all over, i.e. a single unbounded band.
  /// @param uniform the slab it carries everywhere
  explicit SurfaceMaterial(const Acts::MaterialSlab& uniform)
      : bounds{0.f, std::numeric_limits<float>::infinity()}, bands{uniform} {}

  /// A surface whose bands each say what they are made of.
  /// @param bandEdges the band edges, increasing, one more than there are bands
  /// @param bandMaterials what each band is made of
  SurfaceMaterial(std::vector<float> bandEdges,
                  std::vector<Acts::MaterialSlab> bandMaterials);

  /// The same, from the two lengths per band the generated tables carry.
  /// @param composition what the bands share, see `BandComposition`
  /// @param bandEdges the band edges, increasing, one more than there are bands
  /// @param radiationLengths radiation length of each band, zero for a band
  ///        that holds nothing
  /// @param nuclearLengths nuclear interaction length of each band
  SurfaceMaterial(const BandComposition& composition,
                  std::vector<float> bandEdges,
                  const std::vector<float>& radiationLengths,
                  const std::vector<float>& nuclearLengths);

  /// The edges of the bands, increasing: band `k` runs from `bounds[k]` to
  /// `bounds[k + 1]`. Outside the first and the last there is nothing, so a
  /// surface whose material starts away from the beam axis needs no band of
  /// vacuum in front of it.
  std::vector<float> bounds;
  /// What each band is made of, one fewer than there are edges
  std::vector<Acts::MaterialSlab> bands;

  /// Scale how much material the surface carries.
  /// @param scale what to multiply the thickness by
  void scaleThickness(const float scale) {
    for (Acts::MaterialSlab& band : bands) {
      band.scaleThickness(scale);
    }
  }

  /// The edges of one band.
  /// @param band index of the band
  /// @return where it starts and where it ends
  std::pair<float, float> bandBounds(const std::size_t band) const {
    return {bounds[band], bounds[band + 1]};
  }

  /// What a crossing at a position along the surface meets.
  /// @param along r on a disc, |z| on a cylinder
  /// @return the material there, vacuum outside the bands
  const Acts::MaterialSlab& at(const float along) const {
    static const Acts::MaterialSlab nothing = Acts::MaterialSlab::Vacuum(0.f);
    if (bands.empty() || along < bounds.front() || along >= bounds.back()) {
      return nothing;
    }
    const auto edge = std::ranges::upper_bound(bounds, along);
    return bands[edge - bounds.begin() - 1];
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
  /// Index into `DetectorLayout::subsystems` of the subsystem it belongs to, so
  /// that a space point still says where it came from once several of them have
  /// been built together
  std::uint16_t subsystem{};
};

/// The radial bounds of one ring of an endcap disc.
struct RingBounds {
  /// Inner radius
  float rMin{};
  /// Outer radius
  float rMax{};
};

/// One described layer of a barrel: a cylinder and the eta modules along it.
///
/// @note Everything a cylinder needs it carries itself, so that no two lists
///       have to be kept aligned by index.
struct CylinderDescription {
  /// Radius
  float radius{};
  /// Half-length along z, the cylinder being centred on the interaction point
  float halfLengthZ{};
  /// Eta modules it is split into along z
  std::uint32_t modules{1};
  /// Chance that a crossing also leaves a space point on the module next to it;
  /// see `DetectorSurface::overlapProbability`
  float overlapProbability{};
  /// How far that neighbour sits from the one crossed; see
  /// `DetectorSurface::overlapOffset`
  float overlapOffset{};
  /// The layer index this is known by, which is what material keys onto.
  /// Defaults to its position in the barrel it belongs to.
  std::optional<std::uint32_t> layer;
  /// What a crossing meets, resolved in |z|
  SurfaceMaterial material;
};

/// One described layer of an endcap: a disc and the rings across it.
struct DiscDescription {
  /// Absolute longitudinal position, the side coming from the endcap
  float absZ{};
  /// Radial bounds of the rings, increasing. Gaps are allowed and leave no
  /// space point, being the support between two rings.
  std::vector<RingBounds> rings;
  /// @copydoc CylinderDescription::overlapProbability
  float overlapProbability{};
  /// @copydoc CylinderDescription::overlapOffset
  float overlapOffset{};
  /// The layer index this is known by, counting outwards from the interaction
  /// point. Defaults to its position in the endcap it belongs to.
  std::optional<std::uint32_t> layer;
  /// What a crossing meets, resolved in r
  SurfaceMaterial material;
};

/// A passive surface the detector carries away from any sensitive layer: a
/// support, a service disc, a cable run, a beam pipe. Material without readout.
struct PassiveSurfaceDescription {
  /// Cylinder or disc
  SurfaceShape shape{};
  /// Radius for a cylinder, absolute z for a disc
  float refCoord{};
  /// Which sides a disc is built on; a cylinder straddles the interaction point
  /// and has to leave this at its default
  EndcapPlacement placement{EndcapPlacement::Mirrored};
  /// Where its material starts along the coordinate it extends in: |z| for a
  /// cylinder, r for a disc
  float minBound{};
  /// Where it ends, which a service surface spanning the detector leaves at
  /// infinity
  float maxBound{std::numeric_limits<float>::infinity()};
  /// The index this is known by within the passives it is listed with, which is
  /// what material keys onto. Defaults to its position in that list.
  std::optional<std::uint32_t> layer;
  /// What a crossing meets, read off the geometry
  SurfaceMaterial material;
};

/// Cylinders around the beam axis. Nothing about them is uniform: every one has
/// its own half-length, its own material and its own module structure.
struct BarrelDescription {
  /// The cylinders, increasing in radius
  std::vector<CylinderDescription> cylinders;
};

/// Discs normal to the beam axis, on one side of the interaction point or
/// mirrored onto both.
struct EndcapDescription {
  /// Which sides it is built on
  EndcapPlacement placement{EndcapPlacement::Mirrored};
  /// The discs, increasing in absolute z
  std::vector<DiscDescription> discs;
};

/// A named part of the detector, which can be built on its own: the ITk pixels,
/// the ITk strips, one of the ODD's systems. It is what a layer belongs to, and
/// therefore what material and selection key onto.
struct SubsystemDescription {
  /// The name it is selected by, which has to be unique within a detector
  std::string name;
  /// Its barrels
  std::vector<BarrelDescription> barrels;
  /// Its endcaps
  std::vector<EndcapDescription> endcaps;
  /// Its supports and services
  std::vector<PassiveSurfaceDescription> passives;
};

/// A detector as a set of named subsystems, which is the form that is written
/// to and read from file, that material is decorated onto and that subsystems
/// are selected out of. `makeLayout` expands it into the `DetectorLayout` the
/// generator runs on.
struct DetectorDescription {
  /// Passives belonging to no subsystem, the beam pipe being the one that
  /// matters: the only material in front of the innermost layer, and hence the
  /// only source of secondaries there. Kept whichever subsystems are selected.
  std::vector<PassiveSurfaceDescription> passives;
  /// The subsystems, in the order their surfaces are built
  std::vector<SubsystemDescription> subsystems;

  /// @name Containment
  ///
  /// Where the *enclosing* tracker ends, not what this description holds: a
  /// track leaving the pixels curls back through the strips in the same field.
  /// Untouched by a selection, for that reason.
  ///
  /// @{

  /// Radius past which a track has left the tracker for good
  float escapeRadius{1100.f};
  /// Longitudinal half length of the same volume
  float escapeHalfZ{3000.f};

  /// @}
};

/// Which of a subsystem's lists a described layer is in, the other half of
/// naming it.
enum class LayerKind {
  /// A cylinder of one of its barrels
  Barrel,
  /// A disc of one of its endcaps
  Endcap,
  /// One of its supports or services
  Passive,
};

/// Names one described layer of one subsystem, which is what material is keyed
/// onto: it survives a layer being inserted in front of it and a subsystem
/// being selected away, neither of which a position in a list does.
struct LayerId {
  /// The subsystem it belongs to, empty for a passive belonging to the detector
  std::string subsystem;
  /// Which of that subsystem's lists it is in
  LayerKind kind{};
  /// Which endcap of the subsystem, i.e. the sides that endcap is built on.
  /// Only an endcap needs it: a subsystem's barrels share one numbering, and so
  /// do its passives, while a positive-only and a negative-only endcap number
  /// their discs independently.
  EndcapPlacement placement{EndcapPlacement::Mirrored};
  /// The index it answers to within that
  std::uint32_t layer{};

  /// @param other the identifier to compare against
  /// @return whether they name the same layer
  bool operator==(const LayerId& other) const = default;
};

/// What one layer of a detector is made of, away from the description of where
/// that layer is.
struct MaterialEntry {
  /// The layer it belongs on
  LayerId layer;
  /// What a crossing of it meets
  SurfaceMaterial material;
};

/// The material of a whole detector, keyed by the layers of its description, so
/// that where the detector is and what it is made of are written down,
/// generated and reviewed apart from each other.
using MaterialDecoration = std::vector<MaterialEntry>;

/// Write down every layer index a description leaves to the position of its
/// layer in a list, so that what material keys onto is stated rather than
/// derived twice.
///
/// Cylinders are numbered across all of a subsystem's barrels, discs per side
/// of it across all of its endcaps, and passives within the list they are in.
/// `makeLayout` does this itself; `decorate` and `extractMaterial` need it
/// done, and do it.
///
/// @param description the detector to number in place
/// @throws std::invalid_argument if two layers already claim one index
void assignLayerIndices(DetectorDescription& description);

/// Put material onto the layers a decoration names, replacing whatever they
/// carried.
///
/// @param description the detector to decorate in place
/// @param decoration what its layers are made of
/// @throws std::invalid_argument if an entry names a layer the description does
///         not have, which is what catches a material file left behind by a
///         description that has been renumbered
void decorate(DetectorDescription& description,
              const MaterialDecoration& decoration);

/// Take the material off a description, which is the inverse of `decorate` and
/// how a detector read off a geometry is split into the two files it ships as.
///
/// @param description the detector to read
/// @return the material of every layer of it that carries any
MaterialDecoration extractMaterial(const DetectorDescription& description);

/// Leave a description holding where its layers are and nothing about what they
/// are made of.
/// @param description the detector to strip in place
void stripMaterial(DetectorDescription& description);

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

/// The expanded detector the generator runs on, crossing surfaces and reading
/// out layers. `makeLayout` builds it from a `DetectorDescription`.
struct DetectorLayout {
  /// All logical layers, in the index space used by the `layerId` column
  std::vector<DetectorLayer> layers;
  /// The physical surfaces the layers are grouped into
  std::vector<DetectorSurface> surfaces;
  /// Names of the subsystems the layers belong to, indexed by
  /// `DetectorLayer::subsystem`
  std::vector<std::string> subsystems;

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
/// its layers cross-linked. Layers are numbered in call order within the
/// subsystem being built, unless the caller says otherwise.
class DetectorLayoutBuilder {
 public:
  /// @param subsystem the name of the subsystem to start in, empty for a layout
  ///        whose layers belong to nothing in particular
  explicit DetectorLayoutBuilder(std::string subsystem = {});

  /// Start a new subsystem; everything added after this belongs to it and its
  /// layers are numbered from zero again.
  /// @param name the name, which has to be one the layout does not have yet
  /// @return this builder, for chaining
  DetectorLayoutBuilder& beginSubsystem(std::string name);

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
  /// @param layer the layer index to give it, by default the number of
  ///        cylinders added to this subsystem so far
  /// @return this builder, for chaining
  DetectorLayoutBuilder& addCylinder(
      float radius, float halfLengthZ, std::uint32_t numModules,
      std::optional<std::uint32_t> layer = std::nullopt);

  /// Add a sensitive disc on one side, with the rings as given so that they
  /// may leave radial gaps.
  ///
  /// @param side the endcap the disc sits on, which cannot be the barrel
  /// @param absZ the absolute longitudinal position
  /// @param rings the radial bounds of the rings, increasing and non-overlapping
  /// @param layer the layer index to give it, by default the number of discs
  ///        added to this side of this subsystem so far
  /// @return this builder, for chaining
  DetectorLayoutBuilder& addDisc(
      SurfaceSide side, float absZ, std::span<const RingBounds> rings,
      std::optional<std::uint32_t> layer = std::nullopt);

  /// Add a sensitive disc on one side, split into rings tiling `[rMin, rMax]`.
  /// @param side the endcap the disc sits on, which cannot be the barrel
  /// @param absZ the absolute longitudinal position
  /// @param rMin the inner radius
  /// @param rMax the outer radius
  /// @param numRings the number of rings to split it into
  /// @param layer the layer index to give it, by default the number of discs
  ///        added to this side of this subsystem so far
  /// @return this builder, for chaining
  DetectorLayoutBuilder& addDisc(
      SurfaceSide side, float absZ, float rMin, float rMax,
      std::uint32_t numRings,
      std::optional<std::uint32_t> layer = std::nullopt);

  /// Take the layout that has been assembled, leaving the builder empty.
  /// @return the layout
  DetectorLayout build();

 private:
  DetectorLayout m_layout;
  /// Index into `m_layout.subsystems` of the subsystem being built
  std::uint16_t m_subsystem{};
  /// Number of sensitive cylinders added to this subsystem so far
  std::uint32_t m_numCylinders{};
  /// Number of sensitive discs added to this subsystem so far, per side
  std::array<std::uint32_t, 2> m_numDiscs{};
  /// The layer indices handed out already, so that two layers of one subsystem
  /// and side cannot answer to the same one and leave material ambiguous
  std::vector<std::uint64_t> m_layerKeys;

  /// @param side the endcap, which cannot be the barrel
  /// @return a reference to the disc counter of that side
  std::uint32_t& discCounter(SurfaceSide side);

  /// Claim a layer index, or the next free one, for the current subsystem.
  /// @param side the side the layer sits on
  /// @param requested the index the caller asked for, if any
  /// @param counter the running count of layers of this kind
  /// @return the index to use
  /// @throws std::invalid_argument if that index is taken already
  std::uint32_t claimLayer(SurfaceSide side,
                           std::optional<std::uint32_t> requested,
                           std::uint32_t counter);
};

/// Fill in `DetectorSurface::minBound` and `DetectorSurface::maxBound` of every
/// sensitive surface of a layout, a passive one carrying its own. Call it once
/// the layers and the material are in place, and again whenever either changes.
///
/// @param layout the layout to update in place
void updateSurfaceExtents(DetectorLayout& layout);

/// Expand a description into the layout the generator runs on.
///
/// Surfaces come out in the order they are described: the detector's own
/// passives first, then subsystem by subsystem its passives, its barrels and
/// its endcaps.
///
/// @param original the detector to build
/// @return the layout
DetectorLayout makeLayout(const DetectorDescription& original);

/// Keep only the named subsystems of a description, so that a detector held in
/// one piece can be built a system at a time: the ITk pixels alone, or the
/// pixels and the strips together.
///
/// The detector's own passives are kept whichever subsystems are named -- the
/// beam pipe is the only material in front of the innermost layer and dropping
/// it would leave a pixel-only run with no secondaries there -- and so are the
/// escape bounds, a track leaving the pixels still curling back through the
/// strips.
///
/// @param description the detector to narrow
/// @param names the subsystems to keep, in the order they are given
/// @return the narrowed description
/// @throws std::invalid_argument if a name is not one the description has, a
///         typo being far more likely than a detector with a hole in it
DetectorDescription selectSubsystems(const DetectorDescription& description,
                                     std::span<const std::string> names);

/// Build one detector out of several descriptions, appending their subsystems
/// and their detector-level passives in the order given.
///
/// @param descriptions the descriptions to combine
/// @return the combined description, with the widest escape bounds of any of
///         them
/// @throws std::invalid_argument if two of them name the same subsystem
DetectorDescription merge(std::span<const DetectorDescription> descriptions);

/// A coarse stand-in for the ACTS Generic detector's pixel system, and the one
/// detector spelled out in C++ rather than read from file, so that the tests of
/// the generator itself need no data.
/// @return the description of the generic pixel layout
DetectorDescription genericDetectorPixelDescription();

/// Build a stand-in for the ACTS Generic detector's pixel system.
/// @return the generic pixel layout
DetectorLayout makeGenericDetectorPixelLayout();

}  // namespace ActsFatras::Synthetic
