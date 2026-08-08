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
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
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

SurfaceMaterial::SurfaceMaterial(std::vector<float> bandBounds,
                                 std::vector<Acts::MaterialSlab> bandMaterials)
    : bounds(std::move(bandBounds)), bands(std::move(bandMaterials)) {
  if (bounds.size() != bands.size()) {
    throw std::invalid_argument("SurfaceMaterial: one material per band");
  }
}

SurfaceMaterial::SurfaceMaterial(const BandComposition& composition,
                                 std::vector<float> bandBounds,
                                 const std::vector<float>& radiationLengths,
                                 const std::vector<float>& nuclearLengths)
    : SurfaceMaterial(
          std::move(bandBounds),
          bandsFromLengths(composition, radiationLengths, nuclearLengths)) {}

std::uint32_t& DetectorLayoutBuilder::discCounter(const SurfaceSide side) {
  return m_numDiscs[side == SurfaceSide::Positive ? 0 : 1];
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
    const float radius, const float halfLengthZ,
    const std::uint32_t numModules) {
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
    layer.layer = m_numCylinders;
    layer.moduleIndex = m;
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
    const std::span<const RingBounds> rings) {
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
    layer.layer = discCounter(side);
    layer.moduleIndex = static_cast<std::uint32_t>(m);
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
    const float rMax, const std::uint32_t numRings) {
  return addDisc(side, absZ, uniformRings(rMin, rMax, numRings));
}

DetectorLayout DetectorLayoutBuilder::build() {
  DetectorLayout layout = std::move(m_layout);
  m_layout = DetectorLayout{};
  m_numCylinders = 0;
  m_numDiscs[0] = 0;
  m_numDiscs[1] = 0;
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
    for (std::size_t k = 0; k < material.bounds.size(); ++k) {
      if (material.bands[k].isVacuum()) {
        continue;
      }
      minBound = std::min(minBound, k == 0 ? 0.f : material.bounds[k - 1]);
      maxBound = std::max(maxBound, material.bounds[k]);
    }
    surface.minBound = minBound;
    surface.maxBound = maxBound;
  }
}

Synthetic::DetectorLayout Synthetic::makeLayout(
    const BarrelEndcapDescription& description) {
  if (description.barrelRadii.size() != description.barrelHalfLengthsZ.size()) {
    throw std::invalid_argument(
        "makeLayout: one barrel half-length per barrel radius");
  }

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

  if (description.beamPipeRadius.has_value()) {
    // The beam pipe ends where the detector does; any finite bound would do.
    float envelopeZ = 0.f;
    for (const float halfLengthZ : description.barrelHalfLengthsZ) {
      envelopeZ = std::max(envelopeZ, halfLengthZ);
    }
    for (const DiscDescription& disc : description.discs) {
      envelopeZ = std::max(envelopeZ, disc.absZ);
    }
    builder.addPassiveCylinder(
        *description.beamPipeRadius,
        envelopeZ > 0.f ? envelopeZ : *description.beamPipeRadius);
    record(description.beamPipeMaterial);
  }
  for (const PassiveSurfaceDescription& passive : description.passiveSurfaces) {
    constexpr float unbounded = std::numeric_limits<float>::infinity();
    if (passive.shape == SurfaceShape::Cylinder) {
      builder.addPassiveCylinder(
          passive.refCoord,
          passive.maxBound > 0.f ? passive.maxBound : unbounded,
          passive.minBound);
      record(passive.material);
    } else {
      for (const SurfaceSide side :
           {SurfaceSide::Positive, SurfaceSide::Negative}) {
        builder.addPassiveDisc(
            side, passive.refCoord, passive.minBound,
            passive.maxBound > 0.f ? passive.maxBound : unbounded);
        record(passive.material);
      }
    }
  }
  static const SurfaceMaterial vacuum;
  for (std::size_t l = 0; l < description.barrelRadii.size(); ++l) {
    builder.addCylinder(description.barrelRadii[l],
                        description.barrelHalfLengthsZ[l],
                        description.barrelModules);
    record(l < description.barrelMaterials.size()
               ? description.barrelMaterials[l]
               : vacuum,
           l < description.barrelOverlapProbabilities.size()
               ? description.barrelOverlapProbabilities[l]
               : 0.f,
           description.barrelOverlapOffset);
  }
  for (const SurfaceSide side :
       {SurfaceSide::Positive, SurfaceSide::Negative}) {
    for (const DiscDescription& disc : description.discs) {
      builder.addDisc(side, disc.absZ, disc.rings);
      record(disc.material, description.discOverlapProbability,
             description.discOverlapOffset);
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

Synthetic::BarrelEndcapDescription
Synthetic::genericDetectorPixelDescription() {
  // Read off the built `ActsExamples::GenericDetector`. Its sensors sit at the
  // nominal layer radii, so no offset is taken off them.
  BarrelEndcapDescription description;
  // kBeamPipeRadius in GenericDetectorBuilder.hpp
  description.beamPipeRadius = 19.f;
  description.barrelRadii = {32.f, 72.f, 116.f, 172.f};
  description.barrelHalfLengthsZ = {491.f, 491.f, 491.f, 491.f};
  // every barrel cylinder is the same module, and so is every disc
  const SurfaceMaterial barrelMaterial{
      BandComposition{28.03f, 14.f, 7.924e-6f, 1.557f},
      {493.f},
      {95.7f},
      {465.2f}};
  description.barrelMaterials.assign(description.barrelRadii.size(),
                                     barrelMaterial);
  description.beamPipeMaterial = {BandComposition{9.012f, 4.f, 0.07235f, 0.8f},
                                  {3000.f},
                                  {352.8f},
                                  {407.f}};
  const SurfaceMaterial discMaterial{
      BandComposition{28.03f, 14.f, 7.924e-6f, 1.255f},
      {29.f, 177.2f},
      {0.f, 80.04f},
      {0.f, 389.1f}};
  for (const float absZ :
       {600.f, 700.f, 820.f, 960.f, 1100.f, 1300.f, 1500.f}) {
    description.discs.push_back({absZ, {{31.15f, 176.2f}}, discMaterial});
  }
  // Services, whose extent matters as much as their amount: a tube grazed at
  // cosh(eta) is worth ten times its thickness forward.
  const SurfaceMaterial serviceMaterial{
      BandComposition{28.03f, 14.f, 7.924e-6f, 1.5f},
      {29.f, 177.2f},
      {0.f, 95.7f},
      {0.f, 465.2f}};
  for (const float absZ :
       {570.6f, 650.f, 760.f, 890.f, 1030.f, 1200.f, 1400.f, 2252.f}) {
    description.passiveSurfaces.push_back(
        {SurfaceShape::Disc, absZ, 29.f, 177.2f, serviceMaterial});
  }
  description.barrelModules = 1;
  return description;
}

Synthetic::DetectorLayout Synthetic::makeGenericDetectorPixelLayout() {
  return makeLayout(genericDetectorPixelDescription());
}

}  // namespace ActsFatras
