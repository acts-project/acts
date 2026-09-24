// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/DetectorLayout.hpp"

#include "Acts/Material/Material.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>

namespace ActsFatras::Synthetic {

namespace {

void requireSide(const SurfaceSide side) {
  if (side != SurfaceSide::Positive && side != SurfaceSide::Negative) {
    throw std::invalid_argument(
        "DetectorLayoutBuilder: a disc sits on one endcap or the other");
  }
}

/// @param side the side a surface sits on
/// @return the sign it gives the longitudinal coordinate
float sideSign(const SurfaceSide side) {
  return static_cast<float>(static_cast<std::int32_t>(side));
}

/// The sides a placement builds on, in the order they are built.
/// @param placement which sides a described endcap or passive disc sits on
/// @return those sides
std::span<const SurfaceSide> placementSides(const EndcapPlacement placement) {
  static constexpr std::array<SurfaceSide, 2> both{SurfaceSide::Positive,
                                                   SurfaceSide::Negative};
  if (placement == EndcapPlacement::Positive) {
    return {both.data(), 1};
  }
  if (placement == EndcapPlacement::Negative) {
    return {both.data() + 1, 1};
  }
  return both;
}

/// @param side the side a disc sits on, which cannot be the barrel
/// @return which of the two per-side counters is its own
std::size_t sideOrdinal(const SurfaceSide side) {
  return side == SurfaceSide::Positive ? 0u : 1u;
}

/// Visit every described layer of a detector with the identifier it answers to,
/// numbering the ones that leave their index to their position in a list.
///
/// The one traversal there is, so that the numbering material keys onto is
/// decided in a single place: cylinders across all of a subsystem's barrels,
/// discs per side of it across all of its endcaps, passives within their own
/// list.
///
/// @param description the detector to walk, const or not
/// @param visit called with each layer's identifier and the layer itself
template <typename Description, typename Visit>
void walkLayers(Description&& description, Visit visit) {
  const auto walkPassives = [&visit](const std::string& subsystem,
                                     auto& passives) {
    std::uint32_t ordinal = 0;
    for (auto& passive : passives) {
      visit(LayerId{subsystem, LayerKind::Passive, EndcapPlacement::Mirrored,
                    passive.layer.value_or(ordinal)},
            passive);
      ++ordinal;
    }
  };

  walkPassives(std::string{}, description.passives);
  for (auto& subsystem : description.subsystems) {
    walkPassives(subsystem.name, subsystem.passives);

    std::uint32_t cylinders = 0;
    for (auto& barrel : subsystem.barrels) {
      for (auto& cylinder : barrel.cylinders) {
        visit(LayerId{subsystem.name, LayerKind::Barrel,
                      EndcapPlacement::Mirrored,
                      cylinder.layer.value_or(cylinders)},
              cylinder);
        ++cylinders;
      }
    }

    std::array<std::uint32_t, 2> discs{};
    for (auto& endcap : subsystem.endcaps) {
      const std::span<const SurfaceSide> sides =
          placementSides(endcap.placement);
      for (auto& disc : endcap.discs) {
        const std::uint32_t ordinal = discs[sideOrdinal(sides.front())];
        for (const SurfaceSide side : sides) {
          // A mirrored disc is one description and has to answer to one index
          // on both sides, so the two counters cannot have drifted apart. They
          // can only have done so by mixing a one-sided endcap into a subsystem
          // that also has a mirrored one, where the author has to say what they
          // mean.
          if (discs[sideOrdinal(side)] != ordinal) {
            throw std::invalid_argument(
                "assignLayerIndices: subsystem '" + subsystem.name +
                "' has endcaps that number their discs differently on the two "
                "sides; give the discs explicit layer indices");
          }
          ++discs[sideOrdinal(side)];
        }
        visit(LayerId{subsystem.name, LayerKind::Endcap, endcap.placement,
                      disc.layer.value_or(ordinal)},
              disc);
      }
    }
  }
}

/// @param layer the layer to name
/// @return what to call it in a message
std::string describeLayer(const LayerId& layer) {
  const std::string kind = layer.kind == LayerKind::Barrel   ? "barrel"
                           : layer.kind == LayerKind::Endcap ? "endcap"
                                                             : "passive";
  const std::string placement =
      layer.kind != LayerKind::Endcap                ? ""
      : layer.placement == EndcapPlacement::Positive ? " positive"
      : layer.placement == EndcapPlacement::Negative ? " negative"
                                                     : " mirrored";
  return (layer.subsystem.empty() ? "the detector's"
                                  : "'" + layer.subsystem + "'")
      .append(placement)
      .append(" ")
      .append(kind)
      .append(" layer ")
      .append(std::to_string(layer.layer));
}

/// @param names the names given out so far
/// @param name the name to look for
/// @return whether it is one of them
template <typename Names>
bool holdsName(Names&& names, const std::string& name) {
  return std::ranges::find(names, name) != std::ranges::end(names);
}

/// @param subsystems the subsystems to read the names off
/// @return their names, for looking one up among them
auto subsystemNames(const std::vector<SubsystemDescription>& subsystems) {
  return subsystems | std::views::transform(&SubsystemDescription::name);
}

/// Pack a layer index and the side it sits on into one key, so that the indices
/// handed out within one subsystem can be checked for collisions.
/// @param subsystem the subsystem the layer belongs to
/// @param side the side it sits on
/// @param layer the index it answers to
/// @return the key
std::uint64_t layerKey(const std::uint16_t subsystem, const SurfaceSide side,
                       const std::uint32_t layer) {
  const auto sideBits =
      static_cast<std::uint64_t>(static_cast<std::int32_t>(side) + 1);
  return (static_cast<std::uint64_t>(subsystem) << 34) | (sideBits << 32) |
         layer;
}

/// Spell out bands that state only their two lengths.
/// @param composition what the bands share
/// @param radiationLengths radiation length of each band
/// @param nuclearLengths nuclear interaction length of each band
/// @return one slab per band
std::vector<Acts::MaterialSlab> bandsFromLengths(
    const BandComposition& composition,
    const std::vector<float>& radiationLengths,
    const std::vector<float>& nuclearLengths) {
  if (radiationLengths.size() != nuclearLengths.size()) {
    throw std::invalid_argument(
        "SurfaceMaterial: one radiation and one nuclear length per band");
  }
  std::vector<Acts::MaterialSlab> bands;
  bands.reserve(radiationLengths.size());
  for (std::size_t k = 0; k < radiationLengths.size(); ++k) {
    if (!(radiationLengths[k] > 0.f) || !(nuclearLengths[k] > 0.f)) {
      bands.emplace_back(Acts::MaterialSlab::Vacuum(composition.thickness));
      continue;
    }
    // The density follows the radiation length, a band with half the X0 of the
    // surface holding twice the material. Ar and Z say what the stuff is, and
    // that does not vary along a surface.
    const float density = composition.molarDensityX0 / radiationLengths[k];
    bands.emplace_back(Acts::Material::fromMolarDensity(
                           radiationLengths[k], nuclearLengths[k],
                           composition.ar, composition.z, density),
                       composition.thickness);
  }
  return bands;
}

}  // namespace

SurfaceMaterial::SurfaceMaterial(std::vector<float> bandEdges,
                                 std::vector<Acts::MaterialSlab> bandMaterials)
    : bounds(std::move(bandEdges)), bands(std::move(bandMaterials)) {
  if (bands.empty()) {
    if (!bounds.empty()) {
      throw std::invalid_argument(
          "SurfaceMaterial: edges without a band between them");
    }
    return;
  }
  if (bounds.size() != bands.size() + 1) {
    throw std::invalid_argument(
        "SurfaceMaterial: one more band edge than there are bands");
  }
  if (!std::ranges::is_sorted(bounds)) {
    throw std::invalid_argument("SurfaceMaterial: band edges have to increase");
  }
  if (bounds.front() < 0.f) {
    throw std::invalid_argument(
        "SurfaceMaterial: band edges have to be non-negative, a surface "
        "extending in |z| or r");
  }
}

SurfaceMaterial::SurfaceMaterial(const BandComposition& composition,
                                 std::vector<float> bandEdges,
                                 const std::vector<float>& radiationLengths,
                                 const std::vector<float>& nuclearLengths)
    : SurfaceMaterial(
          std::move(bandEdges),
          bandsFromLengths(composition, radiationLengths, nuclearLengths)) {}

DetectorLayoutBuilder::DetectorLayoutBuilder(std::string subsystem) {
  m_layout.subsystems.push_back(std::move(subsystem));
}

DetectorLayoutBuilder& DetectorLayoutBuilder::beginSubsystem(std::string name) {
  if (name.empty()) {
    throw std::invalid_argument(
        "DetectorLayoutBuilder: a subsystem needs a name to be selected by");
  }
  if (holdsName(m_layout.subsystems, name)) {
    throw std::invalid_argument("DetectorLayoutBuilder: subsystem '" + name +
                                "' is one the layout has already, and a name "
                                "has to say which is which");
  }
  // The subsystem a builder starts in is unnamed and holds nothing until
  // something has been added to it, so the first call names that one rather
  // than opening a second. Passives added before it make no layers and belong
  // to no subsystem either way.
  if (m_layout.subsystems.size() == 1 && m_layout.subsystems.front().empty() &&
      m_layout.layers.empty()) {
    m_layout.subsystems.front() = std::move(name);
    return *this;
  }
  m_layout.subsystems.push_back(std::move(name));
  m_subsystem = static_cast<std::uint16_t>(m_layout.subsystems.size() - 1);
  m_numCylinders = 0;
  m_numDiscs = {};
  return *this;
}

std::uint32_t& DetectorLayoutBuilder::discCounter(const SurfaceSide side) {
  return m_numDiscs[side == SurfaceSide::Positive ? 0 : 1];
}

std::uint32_t DetectorLayoutBuilder::claimLayer(
    const SurfaceSide side, const std::optional<std::uint32_t> requested,
    const std::uint32_t counter) {
  const std::uint32_t layer = requested.value_or(counter);
  const std::uint64_t key = layerKey(m_subsystem, side, layer);
  if (std::ranges::find(m_layerKeys, key) != m_layerKeys.end()) {
    throw std::invalid_argument(
        "DetectorLayoutBuilder: two layers of subsystem '" +
        m_layout.subsystems[m_subsystem] + "' answer to layer index " +
        std::to_string(layer) + ", which leaves their material ambiguous");
  }
  m_layerKeys.push_back(key);
  return layer;
}

DetectorLayoutBuilder& DetectorLayoutBuilder::addPassiveCylinder(
    const float radius, const float maxAbsZ, const float minAbsZ) {
  if (radius <= 0.f) {
    throw std::invalid_argument(
        "DetectorLayoutBuilder: cylinder radius has to be positive");
  }
  if (minAbsZ < 0.f || maxAbsZ <= minAbsZ) {
    throw std::invalid_argument(
        "DetectorLayoutBuilder: cylinder z bounds have to be increasing and "
        "non-negative");
  }
  DetectorSurface surface;
  surface.shape = SurfaceShape::Cylinder;
  surface.refCoord = radius;
  surface.minBound = minAbsZ;
  surface.maxBound = maxAbsZ;
  m_layout.surfaces.push_back(std::move(surface));
  return *this;
}

DetectorLayoutBuilder& DetectorLayoutBuilder::addPassiveDisc(
    const SurfaceSide side, const float absZ, const float rMin,
    const float rMax) {
  requireSide(side);
  if (absZ <= 0.f) {
    throw std::invalid_argument(
        "DetectorLayoutBuilder: disc position has to be positive");
  }
  if (rMin < 0.f || rMax <= rMin) {
    throw std::invalid_argument(
        "DetectorLayoutBuilder: disc radii have to be increasing and "
        "non-negative");
  }
  DetectorSurface surface;
  surface.shape = SurfaceShape::Disc;
  surface.refCoord = sideSign(side) * absZ;
  surface.minBound = rMin;
  surface.maxBound = rMax;
  m_layout.surfaces.push_back(std::move(surface));
  return *this;
}

DetectorLayoutBuilder& DetectorLayoutBuilder::addCylinder(
    const float radius, const float halfLengthZ, const std::uint32_t numModules,
    const std::optional<std::uint32_t> layerIndex) {
  if (radius <= 0.f) {
    throw std::invalid_argument(
        "DetectorLayoutBuilder: cylinder radius has to be positive");
  }
  if (halfLengthZ <= 0.f) {
    throw std::invalid_argument(
        "DetectorLayoutBuilder: cylinder half-length has to be positive");
  }
  if (numModules == 0) {
    throw std::invalid_argument(
        "DetectorLayoutBuilder: a cylinder needs at least one module");
  }

  const std::uint32_t index =
      claimLayer(SurfaceSide::Barrel, layerIndex, m_numCylinders);

  DetectorSurface surface;
  surface.shape = SurfaceShape::Cylinder;
  surface.refCoord = radius;
  for (std::uint32_t m = 0; m < numModules; ++m) {
    DetectorLayer layer;
    layer.shape = SurfaceShape::Cylinder;
    layer.refCoord = radius;
    // computed from the ends rather than accumulated, so that the last module
    // lands exactly on halfLengthZ
    layer.minBound = -halfLengthZ + 2.f * halfLengthZ * static_cast<float>(m) /
                                        static_cast<float>(numModules);
    layer.maxBound = -halfLengthZ + 2.f * halfLengthZ *
                                        static_cast<float>(m + 1) /
                                        static_cast<float>(numModules);
    layer.side = SurfaceSide::Barrel;
    layer.layer = index;
    layer.moduleIndex = m;
    layer.subsystem = m_subsystem;
    surface.layers.push_back(
        static_cast<std::uint32_t>(m_layout.layers.size()));
    m_layout.layers.push_back(layer);
  }
  m_layout.surfaces.push_back(std::move(surface));
  ++m_numCylinders;
  return *this;
}

DetectorLayoutBuilder& DetectorLayoutBuilder::addDisc(
    const SurfaceSide side, const float absZ,
    const std::span<const RingBounds> rings,
    const std::optional<std::uint32_t> layerIndex) {
  requireSide(side);
  if (absZ <= 0.f) {
    throw std::invalid_argument(
        "DetectorLayoutBuilder: disc position has to be positive");
  }
  if (rings.empty()) {
    throw std::invalid_argument(
        "DetectorLayoutBuilder: a disc needs at least one ring");
  }
  for (std::size_t m = 0; m < rings.size(); ++m) {
    if (!(rings[m].rMin < rings[m].rMax) || rings[m].rMin < 0.f) {
      throw std::invalid_argument(
          "DetectorLayoutBuilder: a ring needs 0 <= rMin < rMax");
    }
    // gaps between the rings are the point of this overload, overlaps would
    // make `findLayer` pick whichever comes first and are always a mistake
    if (m > 0 && rings[m].rMin < rings[m - 1].rMax) {
      throw std::invalid_argument(
          "DetectorLayoutBuilder: disc rings have to increase in r without "
          "overlapping");
    }
  }

  const std::uint32_t index = claimLayer(side, layerIndex, discCounter(side));

  const float z = sideSign(side) * absZ;
  DetectorSurface surface;
  surface.shape = SurfaceShape::Disc;
  surface.refCoord = z;
  for (std::size_t m = 0; m < rings.size(); ++m) {
    DetectorLayer layer;
    layer.shape = SurfaceShape::Disc;
    layer.refCoord = z;
    layer.minBound = rings[m].rMin;
    layer.maxBound = rings[m].rMax;
    layer.side = side;
    layer.layer = index;
    layer.moduleIndex = static_cast<std::uint32_t>(m);
    layer.subsystem = m_subsystem;
    surface.layers.push_back(
        static_cast<std::uint32_t>(m_layout.layers.size()));
    m_layout.layers.push_back(layer);
  }
  m_layout.surfaces.push_back(std::move(surface));
  ++discCounter(side);
  return *this;
}

DetectorLayoutBuilder& DetectorLayoutBuilder::addDisc(
    const SurfaceSide side, const float absZ, const float rMin,
    const float rMax, const std::uint32_t numRings,
    const std::optional<std::uint32_t> layerIndex) {
  return addDisc(side, absZ, uniformRings(rMin, rMax, numRings), layerIndex);
}

DetectorLayout DetectorLayoutBuilder::build() {
  DetectorLayout layout = std::move(m_layout);
  m_layout = DetectorLayout{};
  m_layout.subsystems.emplace_back();
  m_subsystem = 0;
  m_numCylinders = 0;
  m_numDiscs = {};
  m_layerKeys.clear();
  return layout;
}

}  // namespace ActsFatras::Synthetic

namespace ActsFatras {

Acts::MaterialSlab Synthetic::materialSlab(const float x0, const float l0,
                                           const float ar, const float z,
                                           const float molarDensity,
                                           const float thickness) {
  return Acts::MaterialSlab{
      Acts::Material::fromMolarDensity(x0, l0, ar, z, molarDensity), thickness};
}

std::vector<Synthetic::RingBounds> Synthetic::uniformRings(
    const float rMin, const float rMax, const std::uint32_t numRings) {
  if (!(rMin < rMax) || rMin < 0.f) {
    throw std::invalid_argument("uniformRings: needs 0 <= rMin < rMax");
  }
  if (numRings == 0) {
    throw std::invalid_argument("uniformRings: needs at least one ring");
  }

  std::vector<RingBounds> rings;
  rings.reserve(numRings);
  for (std::uint32_t m = 0; m < numRings; ++m) {
    // computed from the ends rather than accumulated so that the last ring
    // lands exactly on rMax
    rings.emplace_back(rMin + (rMax - rMin) * static_cast<float>(m) /
                                  static_cast<float>(numRings),
                       rMin + (rMax - rMin) * static_cast<float>(m + 1) /
                                  static_cast<float>(numRings));
  }
  return rings;
}

std::vector<Synthetic::RingBounds> Synthetic::subdivideRings(
    const std::span<const RingBounds> rings, const std::uint32_t parts) {
  if (parts <= 1) {
    return {rings.begin(), rings.end()};
  }
  std::vector<RingBounds> split;
  split.reserve(rings.size() * parts);
  for (const RingBounds& ring : rings) {
    const std::vector<RingBounds> pieces =
        uniformRings(ring.rMin, ring.rMax, parts);
    split.insert(split.end(), pieces.begin(), pieces.end());
  }
  return split;
}

void Synthetic::updateSurfaceExtents(DetectorLayout& layout) {
  constexpr float infinity = std::numeric_limits<float>::infinity();
  for (DetectorSurface& surface : layout.surfaces) {
    if (surface.layers.empty()) {
      // a passive surface was declared with its extent and carries no material
      // beyond it
      continue;
    }
    float minBound = infinity;
    float maxBound = -infinity;
    const bool cylinder = surface.shape == SurfaceShape::Cylinder;
    for (const std::uint32_t l : surface.layers) {
      const DetectorLayer& layer = layout.layers[l];
      float lo = layer.minBound;
      float hi = layer.maxBound;
      if (cylinder) {
        // a cylinder is bounded in |z|, so a layer across z = 0 starts there
        const bool spansOrigin = lo <= 0.f && hi >= 0.f;
        const float absLo = std::abs(lo);
        const float absHi = std::abs(hi);
        lo = spansOrigin ? 0.f : std::min(absLo, absHi);
        hi = std::max(absLo, absHi);
      }
      minBound = std::min(minBound, lo);
      maxBound = std::max(maxBound, hi);
    }
    // Material off every layer is kept, so the bands extend the surface: the
    // support between two rings is where a crossing meets no layer and stops
    // anyway.
    const SurfaceMaterial& material = surface.material;
    for (std::size_t k = 0; k < material.bands.size(); ++k) {
      if (material.bands[k].isVacuum()) {
        continue;
      }
      const auto [lo, hi] = material.bandBounds(k);
      minBound = std::min(minBound, lo);
      maxBound = std::max(maxBound, hi);
    }
    surface.minBound = minBound;
    surface.maxBound = maxBound;
  }
}

Synthetic::DetectorLayout Synthetic::makeLayout(
    const DetectorDescription& original) {
  // Numbered up front rather than left to the builder's counters, so that the
  // index a layer is built with is the one its material is keyed by even where
  // the two would count differently -- a mirrored endcap, whose one description
  // becomes a disc on either side.
  DetectorDescription description = original;
  assignLayerIndices(description);

  DetectorLayoutBuilder builder;

  /// What a surface carries, recorded as it is added. Matching afterwards by
  /// reference coordinate cannot tell two surfaces at the same one apart.
  struct Carried {
    const SurfaceMaterial* material;
    float overlapProbability;
    float overlapOffset;
  };
  std::vector<Carried> carried;
  const auto record = [&carried](const SurfaceMaterial& material,
                                 float overlapProbability = 0.f,
                                 float overlapOffset = 0.f) {
    carried.emplace_back(&material, overlapProbability, overlapOffset);
  };

  const auto addPassives =
      [&](const std::vector<PassiveSurfaceDescription>& passives) {
        for (const PassiveSurfaceDescription& passive : passives) {
          if (passive.shape == SurfaceShape::Cylinder) {
            if (passive.placement != EndcapPlacement::Mirrored) {
              throw std::invalid_argument(
                  "makeLayout: a passive cylinder straddles the interaction "
                  "point and cannot sit on one side of it");
            }
            builder.addPassiveCylinder(passive.refCoord, passive.maxBound,
                                       passive.minBound);
            record(passive.material);
            continue;
          }
          for (const SurfaceSide side : placementSides(passive.placement)) {
            builder.addPassiveDisc(side, passive.refCoord, passive.minBound,
                                   passive.maxBound);
            record(passive.material);
          }
        }
      };

  addPassives(description.passives);
  for (const SubsystemDescription& subsystem : description.subsystems) {
    builder.beginSubsystem(subsystem.name);
    addPassives(subsystem.passives);
    for (const BarrelDescription& barrel : subsystem.barrels) {
      for (const CylinderDescription& cylinder : barrel.cylinders) {
        builder.addCylinder(cylinder.radius, cylinder.halfLengthZ,
                            cylinder.modules, cylinder.layer);
        record(cylinder.material, cylinder.overlapProbability,
               cylinder.overlapOffset);
      }
    }
    for (const EndcapDescription& endcap : subsystem.endcaps) {
      // Side-major, so that the discs of one endcap come out in the order they
      // are crossed rather than interleaved between the two sides.
      for (const SurfaceSide side : placementSides(endcap.placement)) {
        for (const DiscDescription& disc : endcap.discs) {
          builder.addDisc(side, disc.absZ, disc.rings, disc.layer);
          record(disc.material, disc.overlapProbability, disc.overlapOffset);
        }
      }
    }
  }

  DetectorLayout layout = builder.build();
  layout.escapeRadius = description.escapeRadius;
  layout.escapeHalfZ = description.escapeHalfZ;

  if (carried.size() != layout.surfaces.size()) {
    throw std::logic_error(
        "makeLayout: a surface was built without recording what it carries");
  }
  for (std::size_t s = 0; s < layout.surfaces.size(); ++s) {
    layout.surfaces[s].material = *carried[s].material;
    layout.surfaces[s].overlapProbability = carried[s].overlapProbability;
    layout.surfaces[s].overlapOffset = carried[s].overlapOffset;
  }
  updateSurfaceExtents(layout);

  return layout;
}

void Synthetic::assignLayerIndices(DetectorDescription& description) {
  std::vector<LayerId> seen;
  walkLayers(description, [&seen](const LayerId& id, auto& layer) {
    if (std::ranges::find(seen, id) != seen.end()) {
      throw std::invalid_argument("assignLayerIndices: two layers answer to " +
                                  describeLayer(id) +
                                  ", which leaves their material ambiguous");
    }
    seen.push_back(id);
    layer.layer = id.layer;
  });
}

void Synthetic::decorate(DetectorDescription& description,
                         const MaterialDecoration& decoration) {
  assignLayerIndices(description);
  for (const MaterialEntry& entry : decoration) {
    bool found = false;
    walkLayers(description, [&](const LayerId& id, auto& layer) {
      if (id == entry.layer) {
        layer.material = entry.material;
        found = true;
      }
    });
    if (!found) {
      throw std::invalid_argument(
          "decorate: this detector has no " + describeLayer(entry.layer) +
          "; the material belongs to a description that has since been "
          "renumbered");
    }
  }
}

Synthetic::MaterialDecoration Synthetic::extractMaterial(
    const DetectorDescription& description) {
  DetectorDescription numbered = description;
  assignLayerIndices(numbered);

  MaterialDecoration decoration;
  walkLayers(std::as_const(numbered),
             [&decoration](const LayerId& id, const auto& layer) {
               if (layer.material.bands.empty()) {
                 return;
               }
               decoration.push_back(MaterialEntry{id, layer.material});
             });
  return decoration;
}

void Synthetic::stripMaterial(DetectorDescription& description) {
  walkLayers(description, [](const LayerId& /*id*/, auto& layer) {
    layer.material = SurfaceMaterial{};
  });
}

Synthetic::DetectorDescription Synthetic::selectSubsystems(
    const DetectorDescription& description,
    const std::span<const std::string> names) {
  DetectorDescription selected;
  // whatever is selected sits inside the same tracker and behind the same beam
  // pipe as the whole detector does
  selected.passives = description.passives;
  selected.escapeRadius = description.escapeRadius;
  selected.escapeHalfZ = description.escapeHalfZ;

  selected.subsystems.reserve(names.size());
  for (const std::string& name : names) {
    const auto subsystem = std::ranges::find(description.subsystems, name,
                                             &SubsystemDescription::name);
    if (subsystem == description.subsystems.end()) {
      std::string known;
      for (const SubsystemDescription& candidate : description.subsystems) {
        known += known.empty() ? "" : ", ";
        known += candidate.name;
      }
      throw std::invalid_argument("selectSubsystems: the description has no '" +
                                  name + "', only " +
                                  (known.empty() ? "nothing at all" : known));
    }
    if (holdsName(subsystemNames(selected.subsystems), name)) {
      throw std::invalid_argument("selectSubsystems: '" + name +
                                  "' was asked for twice");
    }
    selected.subsystems.push_back(*subsystem);
  }
  return selected;
}

Synthetic::DetectorDescription Synthetic::merge(
    const std::span<const DetectorDescription> descriptions) {
  if (descriptions.empty()) {
    return {};
  }

  DetectorDescription merged;
  // The bound is on the tracker these sit inside rather than on any one of
  // them, so it is taken from them rather than left at a default that could be
  // wider than every input.
  merged.escapeRadius = descriptions.front().escapeRadius;
  merged.escapeHalfZ = descriptions.front().escapeHalfZ;

  for (const DetectorDescription& description : descriptions) {
    for (const SubsystemDescription& subsystem : description.subsystems) {
      if (holdsName(subsystemNames(merged.subsystems), subsystem.name)) {
        throw std::invalid_argument(
            "merge: two of these detectors call a subsystem '" +
            subsystem.name + "', and a name has to say which is which");
      }
      merged.subsystems.push_back(subsystem);
    }
    merged.passives.insert(merged.passives.end(), description.passives.begin(),
                           description.passives.end());
    merged.escapeRadius =
        std::max(merged.escapeRadius, description.escapeRadius);
    merged.escapeHalfZ = std::max(merged.escapeHalfZ, description.escapeHalfZ);
  }
  return merged;
}

Synthetic::DetectorDescription Synthetic::genericDetectorPixelDescription() {
  // Read off the built `ActsExamples::GenericDetector`. Its sensors sit at the
  // nominal layer radii, so no offset is taken off them.
  DetectorDescription description;

  // The beam pipe belongs to the detector rather than to the pixels, being in a
  // volume of its own there as it is in a real one. `kBeamPipeRadius` in
  // GenericDetectorBuilder.hpp; it ends where the outermost pixel disc does.
  description.passives.push_back(PassiveSurfaceDescription{
      .shape = SurfaceShape::Cylinder,
      .refCoord = 19.f,
      .maxBound = 1500.f,
      .material = {BandComposition{9.012f, 4.f, 0.07235f, 0.8f},
                   {0.f, 3000.f},
                   {352.8f},
                   {407.f}}});

  SubsystemDescription pixels;
  pixels.name = "generic-pixel";

  // every barrel cylinder is the same module, and so is every disc
  const SurfaceMaterial barrelMaterial{
      BandComposition{28.03f, 14.f, 7.924e-6f, 1.557f},
      {0.f, 493.f},
      {95.7f},
      {465.2f}};
  BarrelDescription barrel;
  for (const float radius : {32.f, 72.f, 116.f, 172.f}) {
    barrel.cylinders.push_back(CylinderDescription{.radius = radius,
                                                   .halfLengthZ = 491.f,
                                                   .modules = 1,
                                                   .material = barrelMaterial});
  }
  pixels.barrels.push_back(std::move(barrel));

  const SurfaceMaterial discMaterial{
      BandComposition{28.03f, 14.f, 7.924e-6f, 1.255f},
      {29.f, 177.2f},
      {80.04f},
      {389.1f}};
  EndcapDescription endcap;
  for (const float absZ :
       {600.f, 700.f, 820.f, 960.f, 1100.f, 1300.f, 1500.f}) {
    endcap.discs.push_back(DiscDescription{
        .absZ = absZ, .rings = {{31.15f, 176.2f}}, .material = discMaterial});
  }
  pixels.endcaps.push_back(std::move(endcap));

  // Services, whose extent matters as much as their amount: a tube grazed at
  // cosh(eta) is worth ten times its thickness forward.
  const SurfaceMaterial serviceMaterial{
      BandComposition{28.03f, 14.f, 7.924e-6f, 1.5f},
      {29.f, 177.2f},
      {95.7f},
      {465.2f}};
  for (const float absZ :
       {570.6f, 650.f, 760.f, 890.f, 1030.f, 1200.f, 1400.f, 2252.f}) {
    pixels.passives.push_back(
        PassiveSurfaceDescription{.shape = SurfaceShape::Disc,
                                  .refCoord = absZ,
                                  .minBound = 29.f,
                                  .maxBound = 177.2f,
                                  .material = serviceMaterial});
  }

  description.subsystems.push_back(std::move(pixels));
  return description;
}

Synthetic::DetectorLayout Synthetic::makeGenericDetectorPixelLayout() {
  return makeLayout(genericDetectorPixelDescription());
}

}  // namespace ActsFatras
