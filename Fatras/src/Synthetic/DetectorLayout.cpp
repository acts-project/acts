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

Synthetic::BarrelEndcapDescription Synthetic::itkPixelDescription() {
  // Positions transcribed from the ATLAS ITKLayouts package, there being no
  // ACTS ITk geometry to reduce: every number is a constant in `*Defines.gmx`.
  BarrelEndcapDescription description;
  // Passive: the layout wants the beam pipe only as the material in front of
  // layer 0, which is the sole source of secondaries there.
  description.beamPipeRadius = 25.f;
  // Stave radii, with the half-length `NrOfModulesPerHalfStave` modules come
  // to. Beyond it the outer layers continue as the inclined rings below.
  description.barrelRadii = {34.f, 99.f, 160.f, 228.f, 291.f};
  description.barrelHalfLengthsZ = {244.2f, 245.4f, 374.6f, 374.6f, 374.6f};
  description.barrelMaterials = {
      {BandComposition{5.86f, 5.2f, 0.05465f, 1.433f},
       {4.9f, 9.8f, 205.3f, 215.1f, 220.f, 224.9f, 229.8f, 244.4f},
       {62.86f, 48.22f, 61.16f, 84.26f, 118.2f, 176.8f, 235.7f, 434.f},
       {144.9f, 117.8f, 142.3f, 209.4f, 320.5f, 542.8f, 711.7f, 923.4f}},
      {BandComposition{6.1f, 5.8f, 0.05216f, 0.9743f},
       {221.7f, 226.6f, 231.5f, 236.4f, 241.4f, 246.3f},
       {38.53f, 52.21f, 65.01f, 92.99f, 253.9f, 3790.f},
       {91.39f, 126.3f, 159.7f, 237.5f, 739.7f, 1.338e4f}},
      {BandComposition{9.107f, 5.1f, 0.02785f, 2.037f},
       {315.4f, 345.4f, 360.4f, 367.9f, 375.5f},
       {50.9f, 74.65f, 111.2f, 176.2f, 527.8f},
       {175.f, 239.5f, 333.2f, 472.f, 902.6f}},
      {BandComposition{12.95f, 6.9f, 0.01758f, 2.073f},
       {7.5f, 345.4f, 360.4f, 367.9f, 375.5f},
       {65.79f, 39.34f, 54.72f, 69.36f, 95.17f},
       {243.f, 162.f, 226.5f, 292.9f, 424.5f}},
      {BandComposition{13.17f, 7.f, 0.01745f, 2.216f},
       {7.5f, 352.9f, 367.9f, 375.5f},
       {68.17f, 40.78f, 61.43f, 92.96f},
       {248.7f, 166.2f, 258.2f, 378.5f}},
  };
  description.beamPipeMaterial = {
      BandComposition{10.6f, 7.6f, 0.04109f, 3.903f},
      {240.f, 1260.f, 1440.f, 2820.f, 3000.f},
      {984.3f, 617.2f, 771.9f, 603.f, 1272.f},
      {1426.f, 1151.f, 1282.f, 1136.f, 1622.f}};
  // Services, whose extent matters as much as their amount: a tube grazed at
  // cosh(eta) is worth ten times its thickness forward.
  description.passiveSurfaces = {
      {SurfaceShape::Cylinder,
       124.f,
       0.f,
       3001.f,
       {BandComposition{12.81f, 5.2f, 0.01398f, 1.826f},
        {120.f, 240.1f, 480.2f, 3001.f},
        {129.1f, 65.99f, 37.2f, 20.94f},
        {463.4f, 276.6f, 182.f, 121.7f}}},
      {SurfaceShape::Cylinder,
       205.1f,
       0.f,
       3001.f,
       {BandComposition{5.961f, 3.8f, 0.0483f, 0.9058f},
        {360.1f, 420.1f, 480.2f, 780.3f, 840.3f, 960.3f, 1020.f, 1080.f, 1140.f,
         1200.f, 1260.f, 3001.f},
        {0.f, 126.3f, 78.58f, 81.17f, 148.5f, 75.25f, 54.33f, 83.59f, 60.71f,
         67.79f, 48.22f, 280.4f},
        {0.f, 354.8f, 242.3f, 286.3f, 398.2f, 299.6f, 243.2f, 264.4f, 223.1f,
         280.1f, 224.f, 492.9f}}},
      {SurfaceShape::Cylinder,
       265.1f,
       0.f,
       3001.f,
       {BandComposition{5.989f, 3.8f, 0.04769f, 0.8821f},
        {360.1f, 480.2f, 540.2f, 600.2f, 720.2f, 840.3f, 900.3f, 960.3f, 1020.f,
         1140.f, 1200.f, 1260.f, 3001.f},
        {0.f, 134.6f, 95.98f, 206.8f, 118.8f, 94.85f, 243.1f, 140.4f, 106.3f,
         90.65f, 98.92f, 239.3f, 154.4f},
        {0.f, 402.f, 271.6f, 390.9f, 379.2f, 280.5f, 508.f, 409.7f, 358.9f,
         274.1f, 347.7f, 503.9f, 397.6f}}},
      {SurfaceShape::Disc,
       378.3f,
       124.f,
       319.1f,
       {BandComposition{6.547f, 4.f, 0.04388f, 1.027f},
        {120.f, 300.1f, 360.1f, 480.2f, 540.2f, 600.2f, 720.2f, 840.3f, 900.3f,
         960.3f, 1020.f, 1140.f, 1200.f, 1260.f, 3001.f},
        {0.f, 59.17f, 281.8f, 156.8f, 114.3f, 210.6f, 138.3f, 110.3f, 288.1f,
         171.3f, 123.8f, 103.6f, 117.5f, 259.9f, 169.5f},
        {0.f, 169.4f, 901.7f, 468.2f, 317.5f, 397.6f, 441.7f, 317.6f, 594.7f,
         486.9f, 418.f, 310.5f, 408.5f, 575.6f, 442.f}}},
  };
  // Coarser than the real detector, whose staves carry nine or twelve modules
  // along z. Only seeders that reason about eta modules see the difference.
  description.barrelModules = 1;
  // Measured on the reference off the module identifiers its clusters carry
  // (`measure_overlaps.py`), layer 0 being low because its staves do not
  // alternate in radius the way the outer ones do.
  description.barrelOverlapProbabilities = {0.02f, 0.15f, 0.16f, 0.13f, 0.13f};
  description.discOverlapProbability = 0.15f;
  // Adjacent staves are eight millimetres apart in radius and the modules of a
  // ring five in z; the track's own slope covers the rest.
  description.barrelOverlapOffset = 7.9f;
  description.discOverlapOffset = 5.f;
  // the strips around these pixels, which a track leaving them curls back
  // through; past this is the ITk's outer support and the calorimeter
  description.escapeRadius = 1000.f;
  description.escapeHalfZ = 3050.f;
  // Seventy-five discs per side carrying ninety-five rings. The radial gaps
  // matter as much as the rings, one being only a module deep.
  description.discs = {
      {263.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{13.6f, 6.8f, 0.01832f, 2.522f},
        {32.2f, 34.f, 50.4f, 52.2f, 54.f, 57.6f, 59.4f, 75.8f, 77.6f, 79.4f,
         81.2f, 117.6f, 119.4f, 121.2f, 123.f},
        {0.f, 129.f, 102.7f, 136.f, 304.5f, 419.7f, 962.1f, 1694.f, 924.7f,
         183.7f, 91.4f, 77.17f, 123.f, 218.7f, 62.52f},
        {0.f, 480.1f, 397.1f, 472.2f, 686.4f, 788.f, 1806.f, 3180.f, 1736.f,
         535.3f, 336.9f, 296.9f, 414.3f, 588.4f, 302.8f}}},
      {291.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{13.58f, 6.8f, 0.01833f, 2.527f},
        {32.2f, 34.f, 50.4f, 52.2f, 54.f, 57.6f, 59.4f, 75.8f, 77.6f, 79.4f,
         81.2f, 117.6f, 119.4f, 121.2f, 123.f},
        {0.f, 132.5f, 102.8f, 130.9f, 297.8f, 420.5f, 919.8f, 1696.f, 1008.f,
         202.1f, 92.38f, 77.21f, 118.8f, 210.8f, 57.92f},
        {0.f, 494.7f, 397.9f, 461.7f, 680.6f, 789.5f, 1727.f, 3184.f, 1893.f,
         565.4f, 340.1f, 297.2f, 404.2f, 560.6f, 276.9f}}},
      {322.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{13.49f, 6.8f, 0.01845f, 2.54f},
        {32.2f, 34.f, 52.2f, 54.f, 57.6f, 59.4f, 75.8f, 77.6f, 79.4f, 81.2f,
         117.6f, 119.4f, 121.2f, 123.f},
        {0.f, 125.4f, 104.9f, 294.8f, 422.6f, 885.1f, 1702.f, 1106.f, 224.4f,
         93.73f, 77.47f, 114.7f, 205.4f, 51.92f},
        {0.f, 420.7f, 403.4f, 679.5f, 793.5f, 1662.f, 3196.f, 2077.f, 599.8f,
         344.3f, 298.3f, 395.f, 541.f, 248.1f}}},
      {357.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{13.27f, 6.8f, 0.0188f, 2.562f},
        {32.2f, 34.f, 52.2f, 54.f, 57.6f, 59.4f, 75.8f, 77.6f, 79.4f, 81.2f,
         117.6f, 119.4f, 121.2f, 123.f},
        {0.f, 108.f, 105.7f, 290.f, 426.2f, 851.f, 1718.f, 1197.f, 247.5f,
         94.81f, 78.14f, 108.6f, 198.f, 46.78f},
        {0.f, 297.f, 406.7f, 678.3f, 800.3f, 1597.f, 3225.f, 2247.f, 633.3f,
         348.2f, 300.9f, 381.6f, 524.4f, 222.7f}}},
      {396.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{13.22f, 6.8f, 0.01886f, 2.569f},
        {32.2f, 34.f, 52.2f, 54.f, 57.6f, 59.4f, 75.8f, 77.6f, 79.4f, 81.2f,
         117.6f, 119.4f, 121.2f, 123.f},
        {0.f, 107.6f, 105.9f, 288.6f, 427.4f, 812.6f, 1724.f, 1293.f, 270.6f,
         95.38f, 78.36f, 103.2f, 198.2f, 42.56f},
        {0.f, 288.3f, 407.8f, 678.f, 802.7f, 1526.f, 3236.f, 2427.f, 663.1f,
         350.1f, 301.8f, 368.9f, 516.f, 201.6f}}},
      {397.4f,
       {{277.f, 317.7f}},
       {BandComposition{9.57f, 4.8f, 0.02596f, 0.0867f},
        {277.3f, 283.f, 283.8f, 298.5f, 299.3f, 300.2f, 304.2f, 305.1f, 310.8f,
         311.6f, 312.4f, 318.1f},
        {0.f, 3.98f, 2.2f, 4.3f, 4.97f, 6.64f, 7.99f, 9.19f, 10.4f, 19.82f,
         1497.f, 0.f},
        {0.f, 14.05f, 8.87f, 15.01f, 17.f, 22.f, 26.46f, 32.14f, 38.86f, 76.82f,
         4999.f, 0.f}}},
      {400.5f,
       {{214.7f, 255.4f}},
       {BandComposition{9.347f, 4.7f, 0.02679f, 0.1765f},
        {214.3f, 249.7f, 250.6f, 256.5f},
        {0.f, 8.22f, 36.06f, 78.27f},
        {0.f, 28.88f, 105.4f, 181.5f}}},
      {404.7f,
       {{148.1f, 188.8f}},
       {BandComposition{9.373f, 4.4f, 0.02624f, 0.344f},
        {146.3f, 182.6f, 183.5f, 184.4f, 189.f, 189.9f, 191.7f},
        {0.f, 15.8f, 12.96f, 7.89f, 5.66f, 10.87f, 47.35f},
        {0.f, 56.5f, 34.37f, 25.94f, 22.25f, 35.84f, 76.21f}}},
      {437.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{13.17f, 6.7f, 0.01889f, 2.578f},
        {32.2f, 34.f, 52.2f, 54.f, 57.6f, 59.4f, 77.6f, 79.4f, 81.2f, 117.6f,
         119.4f, 121.2f, 123.f},
        {0.f, 107.5f, 106.2f, 283.9f, 428.8f, 789.5f, 1687.f, 292.1f, 96.67f,
         78.61f, 99.33f, 192.6f, 36.01f},
        {0.f, 279.6f, 408.9f, 675.f, 805.4f, 1482.f, 3168.f, 689.5f, 354.1f,
         302.8f, 359.5f, 504.8f, 173.2f}}},
      {450.1f,
       {{277.f, 317.7f}},
       {BandComposition{9.566f, 4.9f, 0.02596f, 0.1174f},
        {277.3f, 297.7f, 304.2f, 310.8f, 311.6f, 312.4f, 313.2f, 318.1f},
        {0.f, 5.42f, 3.93f, 5.68f, 8.16f, 49.58f, 1402.f, 0.f},
        {0.f, 18.98f, 13.97f, 20.16f, 31.6f, 200.5f, 4180.f, 0.f}}},
      {458.6f,
       {{214.7f, 255.4f}},
       {BandComposition{9.05f, 4.5f, 0.02792f, 0.2129f},
        {214.3f, 236.2f, 237.1f, 243.8f, 249.7f, 250.6f, 251.4f, 253.1f, 254.f,
         256.5f},
        {0.f, 9.16f, 7.87f, 6.89f, 10.06f, 20.27f, 75.54f, 87.17f, 35.3f,
         28.44f},
        {0.f, 32.28f, 28.04f, 24.75f, 35.33f, 67.43f, 187.1f, 205.f, 60.18f,
         46.95f}}},
      {483.3f,
       {{148.1f, 188.8f}},
       {BandComposition{8.156f, 4.f, 0.03164f, 0.4868f},
        {146.3f, 157.2f, 158.1f, 182.6f, 183.5f, 184.4f, 189.9f, 190.8f,
         191.7f},
        {0.f, 15.26f, 18.99f, 22.12f, 17.15f, 9.56f, 6.68f, 12.21f, 42.54f},
        {0.f, 38.89f, 58.54f, 78.57f, 44.08f, 31.64f, 26.41f, 41.61f, 79.38f}}},
      {486.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{13.1f, 6.7f, 0.01897f, 2.588f},
        {32.2f, 34.f, 52.2f, 54.f, 57.6f, 59.4f, 77.6f, 79.4f, 81.2f, 119.4f,
         121.2f, 123.f},
        {0.f, 106.5f, 106.4f, 283.4f, 430.5f, 766.3f, 1710.f, 317.8f, 97.39f,
         79.58f, 190.1f, 32.09f},
        {0.f, 268.f, 410.f, 676.1f, 808.6f, 1439.f, 3209.f, 719.f, 356.4f,
         305.8f, 496.1f, 154.2f}}},
      {508.9f,
       {{277.f, 317.7f}},
       {BandComposition{9.469f, 4.8f, 0.02634f, 0.1414f},
        {277.3f, 300.2f, 306.7f, 311.6f, 314.9f, 318.1f},
        {0.f, 6.27f, 4.42f, 6.29f, 10.78f, 0.f},
        {0.f, 21.96f, 15.78f, 21.76f, 36.06f, 0.f}}},
      {525.2f,
       {{214.7f, 255.4f}},
       {BandComposition{8.979f, 4.5f, 0.02826f, 0.217f},
        {214.3f, 243.8f, 246.4f, 249.7f, 250.6f, 251.4f, 253.1f, 254.f, 256.5f},
        {0.f, 9.12f, 6.75f, 9.92f, 16.07f, 29.31f, 86.34f, 29.38f, 26.58f},
        {0.f, 32.13f, 24.18f, 34.81f, 54.22f, 92.64f, 202.4f, 48.82f, 43.59f}}},
      {543.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{13.04f, 6.7f, 0.01899f, 2.599f},
        {32.2f, 34.f, 52.2f, 54.f, 57.6f, 59.4f, 77.6f, 79.4f, 81.2f, 119.4f,
         121.2f, 123.f},
        {0.f, 105.3f, 106.6f, 278.9f, 432.1f, 756.4f, 1728.f, 351.4f, 98.84f,
         79.79f, 182.7f, 27.42f},
        {0.f, 257.3f, 411.2f, 672.6f, 811.7f, 1420.f, 3243.f, 754.9f, 360.7f,
         306.7f, 485.8f, 133.5f}}},
      {574.5f,
       {{277.f, 317.7f}},
       {BandComposition{9.416f, 4.8f, 0.02648f, 0.1437f},
        {277.3f, 304.2f, 309.1f, 311.6f, 314.f, 314.9f, 318.1f},
        {0.f, 6.14f, 4.43f, 6.18f, 10.25f, 53.74f, 0.f},
        {0.f, 21.5f, 15.82f, 21.17f, 34.21f, 179.5f, 0.f}}},
      {580.9f,
       {{148.1f, 188.8f}},
       {BandComposition{8.028f, 3.9f, 0.03225f, 0.584f},
        {146.3f, 159.9f, 173.5f, 176.2f, 177.2f, 181.7f, 182.6f, 183.5f, 184.4f,
         190.8f, 191.7f},
        {0.f, 17.25f, 25.21f, 14.36f, 18.68f, 27.04f, 23.56f, 19.07f, 10.18f,
         6.6f, 14.57f},
        {0.f, 42.35f, 89.86f, 52.66f, 67.64f, 96.47f, 70.16f, 47.45f, 33.38f,
         26.34f, 47.35f}}},
      {601.6f,
       {{214.7f, 255.4f}},
       {BandComposition{8.948f, 4.4f, 0.02842f, 0.2193f},
        {214.3f, 249.7f, 250.6f, 251.4f, 252.3f, 253.1f, 254.f, 254.8f, 255.6f,
         256.5f},
        {0.f, 9.16f, 15.24f, 20.13f, 48.23f, 87.49f, 29.56f, 23.51f, 16.12f,
         23.64f},
        {0.f, 32.27f, 51.53f, 66.33f, 135.9f, 202.9f, 49.31f, 38.1f, 33.37f,
         38.32f}}},
      {604.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{13.f, 6.7f, 0.01897f, 2.611f},
        {32.2f, 34.f, 52.2f, 54.f, 57.6f, 59.4f, 77.6f, 79.4f, 81.2f, 119.4f,
         121.2f, 123.f},
        {0.f, 104.1f, 107.f, 273.6f, 434.2f, 739.5f, 1743.f, 378.3f, 99.95f,
         80.11f, 183.1f, 22.17f},
        {0.f, 247.7f, 412.9f, 669.f, 815.6f, 1388.f, 3272.f, 782.4f, 364.4f,
         308.1f, 479.2f, 110.8f}}},
      {647.6f,
       {{277.f, 317.7f}},
       {BandComposition{9.385f, 4.7f, 0.02661f, 0.1383f},
        {277.3f, 309.1f, 310.8f, 311.6f, 313.2f, 314.f, 318.1f},
        {0.f, 5.79f, 4.48f, 5.91f, 8.85f, 21.86f, 0.f},
        {0.f, 20.28f, 15.97f, 20.53f, 28.71f, 73.27f, 0.f}}},
      {675.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{12.93f, 6.7f, 0.01897f, 2.627f},
        {32.2f, 34.f, 52.2f, 54.f, 57.6f, 59.4f, 77.6f, 79.4f, 81.2f, 119.4f,
         121.2f, 123.f},
        {0.f, 101.5f, 107.7f, 263.8f, 436.6f, 718.2f, 1761.f, 410.1f, 101.5f,
         80.48f, 178.4f, 18.07f},
        {0.f, 233.3f, 415.4f, 661.f, 820.3f, 1348.f, 3306.f, 814.4f, 368.9f,
         309.6f, 470.6f, 92.38f}}},
      {689.1f,
       {{214.7f, 255.4f}},
       {BandComposition{9.123f, 4.5f, 0.02758f, 0.231f},
        {214.3f, 249.7f, 250.6f, 251.4f, 252.3f, 253.1f, 254.f, 254.8f, 256.5f},
        {0.f, 9.4f, 15.53f, 17.97f, 28.6f, 81.14f, 28.03f, 22.18f, 5.3f},
        {0.f, 33.2f, 52.58f, 59.92f, 90.79f, 195.f, 46.27f, 35.62f, 17.93f}}},
      {702.1f,
       {{148.1f, 188.8f}},
       {BandComposition{7.937f, 3.9f, 0.03259f, 0.6497f},
        {146.3f, 159.9f, 181.7f, 182.6f, 183.5f, 184.4f, 191.7f},
        {0.f, 17.93f, 28.27f, 25.32f, 19.66f, 9.92f, 6.03f},
        {0.f, 42.23f, 100.8f, 73.84f, 47.12f, 32.35f, 24.27f}}},
      {729.3f,
       {{277.f, 317.7f}},
       {BandComposition{9.392f, 4.7f, 0.02662f, 0.1357f},
        {277.3f, 311.6f, 312.4f, 313.2f, 314.f, 318.1f},
        {0.f, 5.6f, 13.77f, 39.5f, 54.76f, 0.f},
        {0.f, 19.64f, 41.91f, 97.56f, 124.6f, 0.f}}},
      {749.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{12.84f, 6.6f, 0.01901f, 2.643f},
        {32.2f, 34.f, 52.2f, 54.f, 57.6f, 59.4f, 77.6f, 79.4f, 81.2f, 119.4f,
         121.2f, 123.f},
        {0.f, 97.97f, 108.4f, 257.4f, 439.5f, 710.1f, 1779.f, 441.5f, 102.5f,
         80.9f, 173.7f, 15.47f},
        {0.f, 217.1f, 418.1f, 655.5f, 825.6f, 1333.f, 3339.f, 845.3f, 372.6f,
         311.3f, 452.4f, 80.21f}}},
      {789.4f,
       {{214.7f, 255.4f}},
       {BandComposition{9.211f, 4.5f, 0.02715f, 0.2411f},
        {214.3f, 249.7f, 250.6f, 252.3f, 253.1f, 254.f, 254.8f, 255.6f, 256.5f},
        {0.f, 9.64f, 14.41f, 20.1f, 53.49f, 26.61f, 20.69f, 4.61f, 2.95f},
        {0.f, 34.07f, 49.21f, 66.46f, 148.5f, 43.54f, 32.96f, 15.9f, 11.56f}}},
      {820.4f,
       {{277.f, 317.7f}},
       {BandComposition{9.417f, 4.8f, 0.02652f, 0.1431f},
        {277.3f, 311.6f, 312.4f, 313.2f, 314.f, 314.9f, 318.1f},
        {0.f, 5.77f, 8.56f, 13.84f, 44.14f, 59.13f, 0.f},
        {0.f, 20.28f, 29.f, 43.09f, 104.3f, 134.6f, 0.f}}},
      {835.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{12.76f, 6.6f, 0.01909f, 2.658f},
        {32.2f, 34.f, 52.2f, 54.f, 57.6f, 59.4f, 77.6f, 79.4f, 81.2f, 119.4f,
         121.2f, 123.f},
        {0.f, 95.19f, 108.9f, 253.2f, 441.6f, 716.9f, 1789.f, 462.f, 104.5f,
         81.29f, 170.8f, 14.19f},
        {0.f, 205.7f, 420.2f, 652.6f, 829.8f, 1346.f, 3359.f, 869.f, 378.3f,
         312.9f, 450.3f, 73.21f}}},
      {852.4f,
       {{148.1f, 188.8f}},
       {BandComposition{7.558f, 3.7f, 0.03493f, 0.7226f},
        {146.3f, 159.9f, 181.7f, 182.6f, 185.3f, 186.2f, 191.7f},
        {0.f, 18.5f, 30.76f, 26.48f, 21.46f, 9.09f, 5.37f},
        {0.f, 41.7f, 109.6f, 75.06f, 48.38f, 31.73f, 22.08f}}},
      {904.4f,
       {{214.7f, 255.4f}},
       {BandComposition{9.215f, 4.5f, 0.02711f, 0.2508f},
        {214.3f, 249.7f, 250.6f, 252.3f, 253.1f, 254.f, 254.8f, 255.6f, 256.5f},
        {0.f, 9.79f, 14.25f, 19.77f, 34.94f, 24.42f, 19.22f, 4.08f, 2.64f},
        {0.f, 34.66f, 48.82f, 65.42f, 107.7f, 39.57f, 30.38f, 14.23f, 10.4f}}},
      {922.f,
       {{277.f, 317.7f}},
       {BandComposition{9.431f, 4.8f, 0.02645f, 0.1502f},
        {277.3f, 311.6f, 312.4f, 313.2f, 314.f, 314.9f, 315.7f, 318.1f},
        {0.f, 5.93f, 9.28f, 10.46f, 15.18f, 53.06f, 79.72f, 0.f},
        {0.f, 20.86f, 31.18f, 34.68f, 47.81f, 120.8f, 181.4f, 0.f}}},
      {925.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{12.68f, 6.6f, 0.01914f, 2.674f},
        {32.2f, 34.f, 52.2f, 54.f, 57.6f, 59.4f, 77.6f, 79.4f, 81.2f, 119.4f,
         121.2f, 123.f},
        {0.f, 92.54f, 109.6f, 248.5f, 444.5f, 696.8f, 1801.f, 469.5f, 107.7f,
         81.75f, 165.6f, 12.67f},
        {0.f, 194.6f, 422.8f, 649.5f, 835.1f, 1308.f, 3380.f, 882.4f, 387.2f,
         314.7f, 433.2f, 65.76f}}},
      {1026.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{12.6f, 6.6f, 0.01918f, 2.693f},
        {32.2f, 34.f, 52.2f, 54.f, 57.6f, 59.4f, 77.6f, 79.4f, 81.2f, 119.4f,
         121.2f, 123.f},
        {0.f, 89.28f, 110.3f, 244.3f, 447.3f, 697.9f, 1816.f, 480.2f, 109.4f,
         82.3f, 165.2f, 11.26f},
        {0.f, 183.6f, 425.6f, 646.6f, 840.6f, 1310.f, 3408.f, 901.8f, 392.5f,
         316.8f, 426.1f, 58.78f}}},
      {1035.f,
       {{277.f, 317.7f}},
       {BandComposition{9.452f, 4.8f, 0.02636f, 0.1575f},
        {277.3f, 311.6f, 312.4f, 314.f, 314.9f, 315.7f, 316.5f, 318.1f},
        {0.f, 6.08f, 8.58f, 11.25f, 32.45f, 54.45f, 169.2f, 0.f},
        {0.f, 21.44f, 29.1f, 37.25f, 87.81f, 123.9f, 385.f, 0.f}}},
      {1036.f,
       {{214.7f, 255.4f}},
       {BandComposition{9.986f, 5.1f, 0.02405f, 0.2445f},
        {214.3f, 231.2f, 237.1f, 249.7f, 250.6f, 252.3f, 253.1f, 254.8f, 255.6f,
         256.5f},
        {0.f, 9.18f, 6.78f, 9.69f, 13.59f, 19.17f, 22.12f, 82.07f, 4.37f,
         2.35f},
        {0.f, 32.96f, 27.32f, 34.06f, 46.54f, 63.48f, 72.42f, 196.6f, 21.27f,
         10.1f}}},
      {1039.f,
       {{148.1f, 188.8f}},
       {BandComposition{7.95f, 4.1f, 0.03233f, 0.5895f},
        {146.3f, 159.9f, 184.4f, 185.3f, 187.1f, 188.f, 191.7f},
        {0.f, 13.87f, 24.86f, 34.77f, 134.6f, 5.34f, 3.71f},
        {0.f, 29.87f, 88.54f, 122.8f, 396.4f, 26.57f, 19.4f}}},
      {1103.f,
       {{58.1f, 78.5f}},
       {BandComposition{18.57f, 7.7f, 0.0128f, 1.192f},
        {57.7f, 78.2f, 78.7f, 79.5f},
        {0.f, 38.56f, 49.8f, 69.48f},
        {0.f, 157.1f, 180.8f, 320.3f}}},
      {1142.f,
       {{32.6f, 52.9f}, {79.8f, 120.4f}},
       {BandComposition{13.56f, 6.9f, 0.01768f, 3.402f},
        {32.2f, 34.f, 35.8f, 41.3f, 43.1f, 50.4f, 52.2f, 54.f, 59.4f, 61.3f,
         77.6f, 79.4f, 81.2f, 83.1f, 92.1f, 113.9f, 119.4f, 121.2f, 123.f},
        {0.f, 62.67f, 63.31f, 121.4f, 86.06f, 91.53f, 65.69f, 89.31f, 168.7f,
         765.9f, 2296.f, 616.3f, 139.7f, 84.06f, 69.55f, 101.3f, 91.52f, 149.8f,
         12.5f},
        {0.f, 169.7f, 288.6f, 489.f, 372.4f, 389.5f, 296.1f, 364.4f, 608.1f,
         2252.f, 4311.f, 1157.f, 499.5f, 335.f, 283.1f, 374.8f, 299.4f, 354.4f,
         62.81f}}},
      {1146.f,
       {{154.2f, 194.9f}, {214.3f, 255.f}, {274.3f, 315.1f}},
       {BandComposition{8.94f, 5.4f, 0.03341f, 1.681f},
        {153.1f, 160.1f, 161.9f, 172.5f, 175.1f, 193.6f, 195.4f, 196.3f,
         197.1f, 213.1f, 214.9f, 229.8f, 253.5f, 255.2f, 256.1f, 257.f,
         273.1f, 274.f,  289.8f, 312.5f, 315.1f, 316.f,  316.9f},
        {0.f,    85.52f, 93.69f, 118.f,  98.24f, 126.5f, 143.5f, 190.6f,
         210.6f, 0.f,    160.2f, 109.4f, 90.6f,  135.5f, 185.2f, 218.9f,
         0.f,    181.6f, 125.9f, 88.5f,  137.4f, 177.f,  210.2f},
        {0.f,   166.6f, 253.7f, 310.4f, 336.2f, 309.8f, 424.8f, 694.7f,
         781.f, 0.f,    392.8f, 292.5f, 278.9f, 425.8f, 635.1f, 807.f,
         0.f,   434.f,  321.7f, 276.4f, 423.5f, 545.7f, 759.5f}}},
      {1229.f,
       {{58.1f, 78.5f}},
       {BandComposition{18.57f, 7.7f, 0.0128f, 1.192f},
        {57.7f, 78.2f, 78.7f, 79.5f},
        {0.f, 38.57f, 49.46f, 69.61f},
        {0.f, 157.2f, 178.6f, 321.7f}}},
      {1249.f,
       {{154.2f, 194.9f}},
       {BandComposition{9.235f, 5.6f, 0.03487f, 2.f},
        {153.1f, 157.5f, 160.1f, 167.2f, 172.5f, 175.1f, 194.5f, 195.4f, 196.3f,
         197.1f},
        {0.f, 36.04f, 26.14f, 32.21f, 51.81f, 42.81f, 57.28f, 58.55f, 67.58f,
         82.27f},
        {0.f, 69.6f, 52.11f, 72.96f, 160.9f, 176.9f, 157.9f, 234.5f, 375.2f,
         569.3f}}},
      {1272.f,
       {{79.8f, 120.4f}},
       {BandComposition{12.75f, 6.7f, 0.01821f, 2.122f},
        {78.8f, 79.6f, 80.5f, 82.3f, 91.1f, 119.5f, 120.4f, 121.3f, 122.1f,
         123.f},
        {0.f, 3622.f, 94.4f, 65.22f, 42.12f, 61.69f, 86.15f, 151.3f, 17.96f,
         3.99f},
        {0.f, 7241.f, 329.7f, 251.3f, 172.6f, 222.3f, 249.4f, 275.3f, 69.15f,
         22.73f}}},
      {1277.f,
       {{274.3f, 315.1f}},
       {BandComposition{11.56f, 6.2f, 0.02581f, 0.7772f},
        {273.1f, 290.6f, 295.f, 312.5f, 313.4f, 314.2f, 315.1f, 316.f, 316.9f},
        {0.f, 21.69f, 17.16f, 22.31f, 32.82f, 43.01f, 64.83f, 98.08f, 195.8f},
        {0.f, 60.22f, 66.88f, 60.93f, 90.96f, 123.5f, 146.1f, 208.3f, 882.8f}}},
      {1297.f,
       {{214.3f, 255.f}},
       {BandComposition{11.58f, 6.3f, 0.02515f, 0.9907f},
        {213.1f, 214.9f, 230.6f, 235.f, 253.5f, 254.3f, 255.2f, 256.1f, 257.f},
        {0.f, 34.21f, 25.99f, 21.f, 27.16f, 33.52f, 45.75f, 95.91f, 142.9f},
        {0.f, 83.95f, 74.22f, 83.96f, 76.f, 117.2f, 146.9f, 355.2f, 945.1f}}},
      {1359.f,
       {{58.1f, 78.5f}},
       {BandComposition{18.56f, 7.7f, 0.0128f, 1.193f},
        {57.7f, 78.2f, 78.7f, 79.5f},
        {0.f, 38.6f, 49.23f, 68.62f},
        {0.f, 157.3f, 177.f, 313.6f}}},
      {1365.f,
       {{154.2f, 194.9f}},
       {BandComposition{9.44f, 5.7f, 0.0337f, 1.771f},
        {153.1f, 158.3f, 160.1f, 167.2f, 172.5f, 174.2f, 193.6f, 194.5f, 195.4f,
         196.3f, 197.1f},
        {0.f, 30.22f, 23.1f, 30.97f, 45.21f, 34.9f, 50.34f, 37.65f, 42.39f,
         46.38f, 54.22f},
        {0.f, 58.01f, 46.12f, 73.19f, 137.f, 164.3f, 138.9f, 137.2f, 182.4f,
         261.8f, 372.4f}}},
      {1403.f,
       {{79.8f, 120.4f}},
       {BandComposition{12.71f, 6.6f, 0.01811f, 2.146f},
        {78.8f, 79.6f, 80.5f, 82.3f, 91.1f, 119.5f, 120.4f, 121.3f, 122.1f,
         123.f},
        {0.f, 4472.f, 98.32f, 63.18f, 42.69f, 62.37f, 84.3f, 145.4f, 16.61f,
         3.48f},
        {0.f, 8952.f, 340.2f, 245.2f, 174.9f, 224.3f, 245.7f, 260.2f, 64.11f,
         19.97f}}},
      {1427.f,
       {{274.3f, 315.1f}},
       {BandComposition{11.55f, 6.2f, 0.02586f, 0.784f},
        {273.1f, 290.6f, 295.f, 313.4f, 314.2f, 315.1f, 316.f, 316.9f},
        {0.f, 22.06f, 17.04f, 22.6f, 39.85f, 50.83f, 90.49f, 213.8f},
        {0.f, 60.71f, 67.8f, 61.66f, 113.1f, 134.1f, 186.7f, 775.6f}}},
      {1473.f,
       {{214.3f, 255.f}},
       {BandComposition{11.61f, 6.3f, 0.025f, 1.003f},
        {213.1f, 214.9f, 230.6f, 235.f, 254.3f, 255.2f, 256.1f, 257.f},
        {0.f, 35.95f, 26.5f, 20.92f, 27.19f, 31.2f, 80.07f, 144.7f},
        {0.f, 86.4f, 75.12f, 85.46f, 77.12f, 119.8f, 287.4f, 957.3f}}},
      {1492.f,
       {{154.2f, 194.9f}},
       {BandComposition{9.334f, 5.7f, 0.03391f, 1.791f},
        {153.1f, 158.3f, 160.1f, 167.2f, 172.5f, 175.1f, 193.6f, 194.5f, 195.4f,
         196.3f, 197.1f},
        {0.f, 29.65f, 22.79f, 31.35f, 45.97f, 37.01f, 50.97f, 31.21f, 35.25f,
         39.25f, 41.27f},
        {0.f, 56.07f, 45.06f, 74.01f, 137.3f, 161.9f, 139.7f, 118.9f, 160.4f,
         230.1f, 266.5f}}},
      {1503.f,
       {{58.1f, 78.5f}},
       {BandComposition{19.14f, 7.9f, 0.0122f, 1.275f},
        {57.7f, 63.8f, 78.2f, 78.7f, 79.5f},
        {0.f, 27.73f, 41.29f, 52.26f, 73.11f},
        {0.f, 121.8f, 168.6f, 187.6f, 332.9f}}},
      {1553.f,
       {{79.8f, 120.4f}},
       {BandComposition{12.67f, 6.6f, 0.01802f, 2.161f},
        {78.8f, 79.6f, 80.5f, 82.3f, 91.1f, 119.5f, 120.4f, 121.3f, 122.1f,
         123.f},
        {0.f, 5827.f, 101.1f, 64.76f, 42.87f, 62.83f, 83.37f, 136.5f, 15.21f,
         3.1f},
        {0.f, 1.154e4f, 347.7f, 250.5f, 175.6f, 226.1f, 243.8f, 244.3f, 58.95f,
         17.88f}}},
      {1597.f,
       {{274.3f, 315.1f}},
       {BandComposition{11.55f, 6.2f, 0.02587f, 0.7891f},
        {273.1f, 274.9f, 291.5f, 295.f, 313.4f, 314.2f, 315.1f, 316.f, 316.9f},
        {0.f, 27.71f, 21.63f, 16.79f, 22.59f, 32.07f, 44.93f, 76.66f, 262.9f},
        {0.f, 67.49f, 60.83f, 69.17f, 61.8f, 83.88f, 127.6f, 171.f, 723.5f}}},
      {1633.f,
       {{154.2f, 194.9f}},
       {BandComposition{9.23f, 5.7f, 0.03409f, 1.815f},
        {153.1f, 159.2f, 160.1f, 167.2f, 168.f, 171.6f, 175.1f, 192.7f, 193.6f,
         194.5f, 195.4f, 196.3f, 197.1f},
        {0.f, 28.24f, 22.31f, 31.77f, 38.91f, 50.48f, 37.56f, 52.11f, 37.65f,
         26.65f, 29.45f, 31.38f, 33.19f},
        {0.f, 52.98f, 44.01f, 74.81f, 99.97f, 141.2f, 164.3f, 141.9f, 125.f,
         105.7f, 138.8f, 187.f, 214.7f}}},
      {1665.f,
       {{58.1f, 78.5f}},
       {BandComposition{19.14f, 7.9f, 0.0122f, 1.276f},
        {57.7f, 64.3f, 78.2f, 78.7f, 79.5f},
        {0.f, 28.19f, 41.41f, 52.12f, 72.36f},
        {0.f, 123.6f, 169.f, 187.f, 326.5f}}},
      {1676.f,
       {{214.3f, 255.f}},
       {BandComposition{11.64f, 6.3f, 0.02476f, 1.018f},
        {213.1f, 214.9f, 231.5f, 235.f, 252.6f, 254.3f, 255.2f, 256.1f, 257.f},
        {0.f, 37.07f, 26.67f, 20.79f, 27.75f, 20.9f, 24.68f, 39.07f, 146.8f},
        {0.f, 88.46f, 76.45f, 88.08f, 77.56f, 72.28f, 103.8f, 171.5f, 970.5f}}},
      {1721.f,
       {{79.8f, 120.4f}},
       {BandComposition{12.66f, 6.6f, 0.01779f, 2.189f},
        {78.8f, 79.6f, 80.5f, 82.3f, 91.1f, 119.5f, 120.4f, 121.3f, 122.1f,
         123.f},
        {0.f, 9131.f, 104.5f, 66.27f, 43.37f, 63.64f, 81.82f, 132.3f, 12.21f,
         2.64f},
        {0.f, 1.809e4f, 357.2f, 256.f, 177.7f, 229.2f, 240.7f, 233.3f, 49.91f,
         15.34f}}},
      {1789.f,
       {{154.2f, 194.9f}, {274.3f, 315.1f}},
       {BandComposition{10.04f, 5.9f, 0.03018f, 1.002f},
        {153.1f, 159.2f, 160.1f, 168.f,  171.6f, 175.1f, 191.9f, 192.7f,
         193.6f, 194.5f, 195.4f, 196.3f, 197.1f, 273.1f, 274.9f, 291.5f,
         295.f,  313.4f, 314.2f, 315.1f, 316.f,  316.9f},
        {0.f,    29.15f, 25.95f, 34.33f, 53.15f, 39.58f, 53.57f, 41.2f,
         30.71f, 23.63f, 24.68f, 29.19f, 34.17f, 0.f,    72.28f, 55.19f,
         42.18f, 57.14f, 73.49f, 110.5f, 125.f,  437.4f},
        {0.f,    54.91f, 50.4f,  83.12f, 152.1f, 177.2f, 152.1f, 134.9f,
         116.9f, 104.3f, 116.6f, 167.1f, 231.7f, 0.f,    173.2f, 154.1f,
         178.3f, 156.7f, 187.9f, 307.3f, 344.3f, 1414.f}}},
      {1846.f,
       {{58.1f, 78.5f}},
       {BandComposition{19.55f, 8.1f, 0.0118f, 1.337f},
        {57.7f, 59.f, 63.8f, 64.3f, 78.2f, 78.7f, 79.5f},
        {0.f, 33.06f, 21.73f, 29.15f, 43.37f, 54.55f, 75.52f},
        {0.f, 137.1f, 101.1f, 124.2f, 177.f, 195.8f, 340.2f}}},
      {1909.f,
       {{79.8f, 120.4f}},
       {BandComposition{12.66f, 6.6f, 0.01752f, 2.22f},
        {78.8f, 79.6f, 80.5f, 82.3f, 91.1f, 119.5f, 120.4f, 121.3f, 122.1f,
         123.f},
        {0.f, 1.409e4f, 106.8f, 67.78f, 43.96f, 64.54f, 81.59f, 127.1f, 9.98f,
         2.24f},
        {0.f, 2.795e4f, 364.6f, 261.4f, 180.2f, 232.5f, 240.8f, 223.2f, 42.61f,
         13.06f}}},
      {1910.f,
       {{214.3f, 255.f}},
       {BandComposition{11.69f, 6.3f, 0.02443f, 1.031f},
        {213.1f, 214.9f, 231.5f, 232.4f, 235.f, 252.6f, 254.3f, 255.2f, 256.1f,
         257.f},
        {0.f, 38.13f, 27.04f, 23.43f, 20.66f, 28.1f, 17.69f, 19.83f, 23.41f,
         100.1f},
        {0.f, 90.39f, 77.69f, 83.62f, 90.38f, 78.57f, 64.7f, 87.75f, 116.5f,
         660.8f}}},
      {1961.f,
       {{154.2f, 194.9f}},
       {BandComposition{9.16f, 5.7f, 0.03297f, 1.99f},
        {153.1f, 160.1f, 168.f, 171.6f, 175.1f, 191.f, 192.7f, 193.6f, 195.4f,
         196.3f, 197.1f},
        {0.f, 27.53f, 33.96f, 52.98f, 39.07f, 53.82f, 33.69f, 24.47f, 19.58f,
         22.45f, 30.92f},
        {0.f, 51.31f, 82.04f, 150.4f, 176.2f, 151.7f, 121.1f, 101.4f, 93.67f,
         131.6f, 209.7f}}},
      {2007.f,
       {{274.3f, 315.1f}},
       {BandComposition{11.55f, 6.2f, 0.02585f, 0.7972f},
        {273.1f, 274.9f, 291.5f, 295.f, 314.2f, 315.1f, 316.f, 316.9f},
        {0.f, 29.35f, 22.05f, 16.55f, 22.81f, 38.29f, 46.59f, 160.f},
        {0.f, 69.58f, 61.2f, 71.39f, 62.5f, 101.7f, 135.1f, 459.8f}}},
      {2120.f,
       {{79.8f, 120.4f}},
       {BandComposition{12.37f, 6.5f, 0.01718f, 1.799f},
        {78.8f, 79.6f, 80.5f, 119.5f, 120.4f, 121.3f, 122.1f, 123.f},
        {0.f, 1.636e4f, 89.01f, 55.02f, 83.19f, 139.5f, 6.99f, 1.58f},
        {0.f, 3.279e4f, 301.6f, 211.8f, 293.3f, 233.f, 31.56f, 9.38f}}},
      {2151.f,
       {{154.2f, 194.9f}},
       {BandComposition{9.092f, 5.7f, 0.03273f, 2.039f},
        {153.1f, 160.1f, 168.f, 172.5f, 175.1f, 191.f, 192.7f, 193.6f, 195.4f,
         196.3f, 197.1f},
        {0.f, 27.18f, 34.66f, 51.99f, 38.75f, 55.15f, 27.61f, 18.71f, 16.97f,
         18.2f, 27.78f},
        {0.f, 49.89f, 83.56f, 156.8f, 185.5f, 155.4f, 109.2f, 86.89f, 82.93f,
         108.4f, 188.f}}},
      {2180.f,
       {{214.3f, 255.f}},
       {BandComposition{11.77f, 6.3f, 0.02396f, 1.048f},
        {213.1f, 214.9f, 231.5f, 235.f, 252.6f, 255.2f, 256.1f, 257.f},
        {0.f, 39.56f, 27.73f, 20.77f, 28.55f, 15.4f, 17.96f, 35.09f},
        {0.f, 92.95f, 78.49f, 92.16f, 79.82f, 62.91f, 92.08f, 231.5f}}},
      {2253.f,
       {{274.3f, 315.1f}},
       {BandComposition{11.56f, 6.2f, 0.02588f, 0.8005f},
        {273.1f, 274.9f, 291.5f, 295.9f, 314.2f, 315.1f, 316.f, 316.9f},
        {0.f, 29.77f, 22.27f, 17.27f, 22.96f, 33.88f, 47.83f, 94.02f},
        {0.f, 70.23f, 61.49f, 70.03f, 62.63f, 86.85f, 135.3f, 323.1f}}},
      {2357.f,
       {{79.8f, 120.4f}},
       {BandComposition{12.48f, 6.3f, 0.01891f, 3.501f},
        {78.8f, 79.6f, 80.5f, 81.4f, 97.3f, 108.9f, 119.5f, 120.4f, 121.3f,
         122.1f, 123.f},
        {0.f, 5.724e4f, 178.3f, 113.f, 80.31f, 65.17f, 83.46f, 110.f, 152.6f,
         12.37f, 2.73f},
        {0.f, 1.152e5f, 599.f, 431.2f, 252.8f, 207.1f, 293.3f, 308.7f, 266.7f,
         52.73f, 15.94f}}},
      {2361.f,
       {{154.2f, 194.9f}},
       {BandComposition{9.048f, 5.7f, 0.03234f, 2.092f},
        {153.1f, 160.1f, 168.f, 172.5f, 175.1f, 191.f, 192.7f, 196.3f, 197.1f},
        {0.f, 26.88f, 35.39f, 53.61f, 39.4f, 56.62f, 23.02f, 14.89f, 22.63f},
        {0.f, 48.61f, 85.18f, 160.f, 190.7f, 159.6f, 98.45f, 77.89f, 152.5f}}},
      {2491.f,
       {{214.3f, 255.f}},
       {BandComposition{11.86f, 6.3f, 0.0234f, 1.064f},
        {213.1f, 214.f, 214.9f, 231.5f, 232.4f, 235.f, 252.6f, 253.5f, 255.2f,
         256.1f, 257.f},
        {0.f, 43.61f, 33.81f, 28.28f, 23.94f, 20.64f, 29.f, 14.92f, 12.13f,
         14.38f, 18.39f},
        {0.f, 98.36f, 87.31f, 80.23f, 86.49f, 95.03f, 81.07f, 58.22f, 54.33f,
         75.96f, 121.1f}}},
      {2533.f,
       {{274.3f, 315.1f}},
       {BandComposition{11.66f, 6.2f, 0.02511f, 0.8035f},
        {273.1f, 274.9f, 292.4f, 295.f, 314.2f, 315.1f, 316.f, 316.9f},
        {0.f, 30.2f, 21.13f, 15.46f, 21.74f, 29.03f, 42.2f, 60.57f},
        {0.f, 70.96f, 60.67f, 72.3f, 61.55f, 75.62f, 125.9f, 233.1f}}},
      {2593.f,
       {{154.2f, 194.9f}},
       {BandComposition{9.001f, 5.7f, 0.03192f, 2.15f},
        {153.1f, 160.1f, 168.f, 172.5f, 175.1f, 191.f, 191.9f, 192.7f, 193.6f,
         196.3f, 197.1f},
        {0.f, 26.39f, 36.31f, 54.87f, 40.33f, 58.14f, 20.6f, 17.28f, 14.06f,
         12.55f, 18.09f},
        {0.f, 47.04f, 87.37f, 161.4f, 197.7f, 163.9f, 92.31f, 81.75f, 72.28f,
         68.52f, 121.4f}}},
      {2621.f,
       {{79.8f, 120.4f}},
       {BandComposition{12.32f, 6.4f, 0.01694f, 1.871f},
        {78.8f, 79.6f, 80.5f, 119.5f, 120.4f, 121.3f, 122.1f, 123.f},
        {0.f, 7.403e4f, 95.47f, 57.23f, 81.05f, 117.9f, 5.75f, 1.28f},
        {0.f, 1.609e5f, 321.6f, 220.3f, 287.4f, 194.f, 26.15f, 7.59f}}},
      {2850.f,
       {{154.2f, 194.9f}, {214.3f, 255.f}, {274.3f, 315.1f}},
       {BandComposition{10.22f, 6.f, 0.02795f, 1.005f},
        {153.1f, 160.1f, 168.f,  172.5f, 175.1f, 192.7f, 193.6f,
         196.3f, 197.1f, 213.1f, 214.f,  215.7f, 232.4f, 235.f,
         252.6f, 253.5f, 255.2f, 256.1f, 257.f,  273.1f, 274.9f,
         292.4f, 295.f,  295.9f, 314.2f, 315.1f, 316.f,  316.9f},
        {0.f,    35.28f, 50.86f, 76.29f, 56.08f, 81.46f, 19.96f,
         14.57f, 18.79f, 0.f,    137.3f, 88.54f, 79.55f, 57.4f,
         82.14f, 38.04f, 27.96f, 32.58f, 35.82f, 0.f,    115.6f,
         79.66f, 58.9f,  69.4f,  81.71f, 104.5f, 146.6f, 200.8f},
        {0.f,    61.95f, 122.4f, 221.6f, 278.3f, 229.7f, 99.39f,
         82.27f, 125.5f, 0.f,    293.f,  237.8f, 227.2f, 274.f,
         229.7f, 154.4f, 130.8f, 178.1f, 235.8f, 0.f,    269.1f,
         227.6f, 268.8f, 246.7f, 230.7f, 268.2f, 432.5f, 739.5f}}},
  };

  return description;
}

Synthetic::DetectorLayout Synthetic::makeItkPixelLayout() {
  return makeLayout(itkPixelDescription());
}

Synthetic::BarrelEndcapDescription
Synthetic::openDataDetectorPixelDescription() {
  // Read off the pixel volumes of the geometry ACTS builds for the ODD, and
  // kept here so that the preset works without DD4hep. Positions are of the
  // sensitive silicon rather than of the layers the ODD declares.
  BarrelEndcapDescription description;
  // Beryllium beam pipe at 23.6 to 24.4 mm. It sits in a volume of its own that
  // the pixel selector excludes, so the reduction is told to add it.
  description.beamPipeRadius = 24.f;
  // The carbon fibre cylinder each barrel layer is carried on, transcribed
  // from `OpenDataPixels.xml`: an ACTS geometry keeps material for its layers
  // alone. All four share a stave of 14 modules 72 mm long.
  description.barrelRadii = {32.21f, 68.21f, 114.2f, 170.2f};
  description.barrelHalfLengthsZ = {507.2f, 507.2f, 507.2f, 507.2f};
  description.barrelMaterials = {
      {BandComposition{14.59f, 6.5f, 0.02114f, 5.012f},
       {491.8f, 512.2f},
       {298.9f, 237.8f},
       {810.2f, 793.2f}},
      {BandComposition{13.64f, 6.3f, 0.02486f, 5.086f},
       {491.8f, 502.f, 512.2f},
       {368.6f, 411.5f, 326.4f},
       {868.7f, 1804.f, 1833.f}},
      {BandComposition{14.17f, 6.4f, 0.02198f, 4.679f},
       {276.6f, 481.5f, 491.8f, 502.f, 512.2f},
       {320.4f, 290.f, 367.8f, 256.4f, 183.4f},
       {806.f, 786.8f, 822.8f, 1100.f, 1346.f}},
      {BandComposition{14.74f, 6.5f, 0.01995f, 4.093f},
       {235.6f, 338.1f, 348.3f, 399.6f, 481.5f, 491.8f, 502.f, 512.2f},
       {294.2f, 261.5f, 208.9f, 247.3f, 199.6f, 136.3f, 102.4f, 114.5f},
       {719.8f, 697.3f, 667.2f, 690.9f, 649.8f, 570.7f, 640.4f, 921.f}},
  };
  description.beamPipeMaterial = {
      BandComposition{9.046f, 4.1f, 0.07136f, 0.1641f},
      {800.2f, 4001.f},
      {71.03f, 0.f},
      {80.26f, 0.f}};
  // Services, whose extent matters as much as their amount: a tube grazed at
  // cosh(eta) is worth ten times its thickness forward.
  description.passiveSurfaces = {
      {SurfaceShape::Disc,
       559.f,
       29.9f,
       190.f,
       {BandComposition{26.91f, 8.1f, 0.006668f, 1.739f},
        {29.9f, 33.1f, 36.3f, 42.7f, 74.7f, 77.9f, 81.1f, 122.8f, 126.f, 129.2f,
         132.4f, 183.6f, 186.8f, 190.f},
        {0.f, 0.f, 1181.f, 6.9f, 0.f, 30.94f, 25.52f, 0.f, 350.9f, 32.03f,
         42.11f, 0.f, 124.8f, 43.37f},
        {0.f, 0.f, 5737.f, 48.06f, 0.f, 137.3f, 107.4f, 0.f, 1777.f, 182.2f,
         211.5f, 0.f, 929.9f, 264.f}}},
      {SurfaceShape::Disc,
       1950.f,
       40.f,
       189.f,
       {BandComposition{22.96f, 8.3f, 0.007943f, 10.24f},
        {40.f, 69.8f, 135.4f, 138.3f, 141.3f, 144.3f, 162.2f, 165.2f, 177.1f,
         180.1f, 183.f, 189.f},
        {0.f, 0.f, 287.6f, 17.44f, 6.03f, 7.9f, 287.6f, 47.42f, 296.8f, 29.97f,
         7.76f, 12.99f},
        {0.f, 0.f, 521.8f, 139.8f, 58.04f, 73.59f, 521.8f, 273.8f, 538.6f,
         208.2f, 72.56f, 108.f}}},
  };
  // left unsplit as in the ITk layout above; the real detector has fourteen
  // modules along z
  description.barrelModules = 1;
  // Measured on the ColliderML reference through its surface identifiers: a
  // crossing above a GeV leaves 1.24 clusters, flat across the detector.
  description.barrelOverlapProbabilities = {0.24f, 0.24f, 0.24f, 0.24f};
  description.discOverlapProbability = 0.24f;
  description.barrelOverlapOffset = 1.65f;
  description.discOverlapOffset = 1.2f;
  // the ODD's short and long strips, out to the solenoid
  description.escapeRadius = 1100.f;
  description.escapeHalfZ = 3000.f;
  // Seven layers per side, each *two* rings rather than one disc: they leave
  // no radial gap, and the stagger in z is all the structure buys here.
  description.discs = {
      {616.2f,
       {{42.85f, 111.f}},
       {BandComposition{12.83f, 6.3f, 0.02915f, 1.608f},
        {32.f, 35.f, 47.2f, 50.2f, 68.4f, 71.5f, 86.6f, 89.7f, 92.7f, 120.f,
         123.1f, 144.3f, 147.4f, 168.6f, 171.6f, 174.7f, 177.7f, 183.8f},
        {0.f, 0.f, 207.5f, 90.89f, 193.1f, 113.5f, 195.3f, 97.54f, 123.f,
         186.2f, 112.6f, 190.3f, 109.1f, 191.2f, 93.16f, 187.2f, 783.5f, 0.f},
        {0.f, 0.f, 389.6f, 289.4f, 373.3f, 318.9f, 374.9f, 301.2f, 327.9f,
         364.3f, 315.9f, 369.7f, 312.1f, 371.4f, 294.3f, 390.9f, 1422.f, 0.f}}},
      {623.2f,
       {{105.5f, 173.8f}},
       {BandComposition{12.83f, 6.3f, 0.02915f, 1.608f},
        {32.f, 35.f, 47.2f, 50.2f, 68.4f, 71.5f, 86.6f, 89.7f, 92.7f, 120.f,
         123.1f, 144.3f, 147.4f, 168.6f, 171.6f, 174.7f, 177.7f, 183.8f},
        {0.f, 0.f, 207.5f, 90.89f, 193.1f, 113.5f, 195.3f, 97.54f, 123.f,
         186.2f, 112.6f, 190.3f, 109.1f, 191.2f, 93.16f, 187.2f, 783.5f, 0.f},
        {0.f, 0.f, 389.6f, 289.4f, 373.3f, 318.9f, 374.9f, 301.2f, 327.9f,
         364.3f, 315.9f, 369.7f, 312.1f, 371.4f, 294.3f, 390.9f, 1422.f, 0.f}}},
      {716.2f,
       {{42.85f, 111.f}},
       {BandComposition{12.81f, 6.3f, 0.02928f, 1.604f},
        {32.f, 35.f, 47.2f, 50.2f, 68.4f, 71.5f, 86.6f, 89.7f, 92.7f, 120.f,
         123.1f, 144.3f, 147.4f, 168.6f, 171.6f, 174.7f, 177.7f, 183.8f},
        {0.f, 0.f, 206.6f, 87.77f, 189.6f, 113.5f, 194.8f, 102.4f, 116.8f,
         185.5f, 110.4f, 190.3f, 111.3f, 192.1f, 107.f, 224.4f, 791.2f, 0.f},
        {0.f, 0.f, 388.1f, 284.8f, 370.9f, 318.3f, 374.f, 306.5f, 321.5f,
         363.2f, 312.8f, 370.f, 313.3f, 372.7f, 315.9f, 407.7f, 1436.f, 0.f}}},
      {723.2f,
       {{105.5f, 173.8f}},
       {BandComposition{12.81f, 6.3f, 0.02928f, 1.604f},
        {32.f, 35.f, 47.2f, 50.2f, 68.4f, 71.5f, 86.6f, 89.7f, 92.7f, 120.f,
         123.1f, 144.3f, 147.4f, 168.6f, 171.6f, 174.7f, 177.7f, 183.8f},
        {0.f, 0.f, 206.6f, 87.77f, 189.6f, 113.5f, 194.8f, 102.4f, 116.8f,
         185.5f, 110.4f, 190.3f, 111.3f, 192.1f, 107.f, 224.4f, 791.2f, 0.f},
        {0.f, 0.f, 388.1f, 284.8f, 370.9f, 318.3f, 374.f, 306.5f, 321.5f,
         363.2f, 312.8f, 370.f, 313.3f, 372.7f, 315.9f, 407.7f, 1436.f, 0.f}}},
      {836.2f,
       {{42.85f, 111.f}},
       {BandComposition{12.83f, 6.3f, 0.02919f, 1.607f},
        {32.f, 35.f, 47.2f, 50.2f, 68.4f, 71.5f, 86.6f, 92.7f, 120.f, 123.1f,
         144.3f, 147.4f, 168.6f, 171.6f, 174.7f, 177.7f, 183.8f},
        {0.f, 0.f, 206.2f, 86.92f, 187.7f, 114.3f, 195.2f, 106.2f, 186.7f,
         107.2f, 190.f, 113.6f, 192.8f, 106.3f, 224.f, 672.9f, 0.f},
        {0.f, 0.f, 387.7f, 284.1f, 370.5f, 319.5f, 374.7f, 311.1f, 365.4f,
         309.9f, 369.4f, 316.7f, 373.9f, 309.f, 407.7f, 1221.f, 0.f}}},
      {843.2f,
       {{105.5f, 173.8f}},
       {BandComposition{12.83f, 6.3f, 0.02919f, 1.607f},
        {32.f, 35.f, 47.2f, 50.2f, 68.4f, 71.5f, 86.6f, 92.7f, 120.f, 123.1f,
         144.3f, 147.4f, 168.6f, 171.6f, 174.7f, 177.7f, 183.8f},
        {0.f, 0.f, 206.2f, 86.92f, 187.7f, 114.3f, 195.2f, 106.2f, 186.7f,
         107.2f, 190.f, 113.6f, 192.8f, 106.3f, 224.f, 672.9f, 0.f},
        {0.f, 0.f, 387.7f, 284.1f, 370.5f, 319.5f, 374.7f, 311.1f, 365.4f,
         309.9f, 369.4f, 316.7f, 373.9f, 309.f, 407.7f, 1221.f, 0.f}}},
      {976.2f,
       {{42.85f, 111.f}},
       {BandComposition{12.83f, 6.3f, 0.02916f, 1.602f},
        {32.f, 35.f, 47.2f, 50.2f, 53.2f, 68.4f, 71.5f, 86.6f, 92.7f, 120.f,
         123.1f, 144.3f, 147.4f, 168.6f, 171.6f, 174.7f, 177.7f, 183.8f},
        {0.f, 0.f, 205.5f, 87.01f, 150.9f, 194.5f, 114.3f, 194.7f, 106.8f,
         185.4f, 101.9f, 191.8f, 113.5f, 194.1f, 109.4f, 228.4f, 680.8f, 0.f},
        {0.f, 0.f, 386.6f, 283.7f, 346.8f, 373.4f, 318.9f, 373.6f, 311.1f,
         363.f, 303.6f, 372.9f, 316.1f, 376.2f, 311.8f, 417.1f, 1235.f, 0.f}}},
      {983.2f,
       {{105.5f, 173.8f}},
       {BandComposition{12.83f, 6.3f, 0.02916f, 1.602f},
        {32.f, 35.f, 47.2f, 50.2f, 53.2f, 68.4f, 71.5f, 86.6f, 92.7f, 120.f,
         123.1f, 144.3f, 147.4f, 168.6f, 171.6f, 174.7f, 177.7f, 183.8f},
        {0.f, 0.f, 205.5f, 87.01f, 150.9f, 194.5f, 114.3f, 194.7f, 106.8f,
         185.4f, 101.9f, 191.8f, 113.5f, 194.1f, 109.4f, 228.4f, 680.8f, 0.f},
        {0.f, 0.f, 386.6f, 283.7f, 346.8f, 373.4f, 318.9f, 373.6f, 311.1f,
         363.f, 303.6f, 372.9f, 316.1f, 376.2f, 311.8f, 417.1f, 1235.f, 0.f}}},
      {1116.f,
       {{42.85f, 111.f}},
       {BandComposition{12.88f, 6.3f, 0.02887f, 1.479f},
        {32.f, 41.1f, 47.2f, 50.2f, 53.2f, 68.4f, 71.5f, 86.6f, 92.7f, 120.f,
         123.1f, 144.3f, 147.4f, 168.6f, 171.6f, 174.7f, 177.7f, 183.8f},
        {0.f, 0.f, 174.4f, 81.57f, 130.6f, 179.6f, 105.4f, 179.7f, 95.56f,
         171.5f, 92.2f, 174.5f, 108.f, 177.1f, 102.1f, 202.5f, 733.8f, 0.f},
        {0.f, 0.f, 338.5f, 263.9f, 313.6f, 344.9f, 294.3f, 345.f, 283.9f,
         335.4f, 278.2f, 339.5f, 300.9f, 343.6f, 289.6f, 371.2f, 1332.f, 0.f}}},
      {1123.f,
       {{105.5f, 173.8f}},
       {BandComposition{12.88f, 6.3f, 0.02887f, 1.479f},
        {32.f, 41.1f, 47.2f, 50.2f, 53.2f, 68.4f, 71.5f, 86.6f, 92.7f, 120.f,
         123.1f, 144.3f, 147.4f, 168.6f, 171.6f, 174.7f, 177.7f, 183.8f},
        {0.f, 0.f, 174.4f, 81.57f, 130.6f, 179.6f, 105.4f, 179.7f, 95.56f,
         171.5f, 92.2f, 174.5f, 108.f, 177.1f, 102.1f, 202.5f, 733.8f, 0.f},
        {0.f, 0.f, 338.5f, 263.9f, 313.6f, 344.9f, 294.3f, 345.f, 283.9f,
         335.4f, 278.2f, 339.5f, 300.9f, 343.6f, 289.6f, 371.2f, 1332.f, 0.f}}},
      {1316.f,
       {{42.85f, 111.f}},
       {BandComposition{12.92f, 6.4f, 0.02867f, 1.353f},
        {32.f, 47.2f, 50.2f, 53.2f, 68.4f, 71.5f, 86.6f, 92.7f, 120.f, 123.1f,
         144.3f, 147.4f, 168.6f, 171.6f, 174.7f, 177.7f, 183.8f},
        {0.f, 0.f, 74.11f, 116.6f, 164.3f, 96.68f, 164.4f, 91.73f, 155.4f,
         78.5f, 159.5f, 95.35f, 161.1f, 93.83f, 185.2f, 661.5f, 0.f},
        {0.f, 0.f, 240.5f, 285.2f, 315.5f, 269.5f, 315.6f, 264.4f, 306.1f,
         247.5f, 310.4f, 266.f, 312.7f, 264.6f, 339.5f, 1200.f, 0.f}}},
      {1323.f,
       {{105.5f, 173.8f}},
       {BandComposition{12.92f, 6.4f, 0.02867f, 1.353f},
        {32.f, 47.2f, 50.2f, 53.2f, 68.4f, 71.5f, 86.6f, 92.7f, 120.f, 123.1f,
         144.3f, 147.4f, 168.6f, 171.6f, 174.7f, 177.7f, 183.8f},
        {0.f, 0.f, 74.11f, 116.6f, 164.3f, 96.68f, 164.4f, 91.73f, 155.4f,
         78.5f, 159.5f, 95.35f, 161.1f, 93.83f, 185.2f, 661.5f, 0.f},
        {0.f, 0.f, 240.5f, 285.2f, 315.5f, 269.5f, 315.6f, 264.4f, 306.1f,
         247.5f, 310.4f, 266.f, 312.7f, 264.6f, 339.5f, 1200.f, 0.f}}},
      {1516.f,
       {{42.85f, 111.f}},
       {BandComposition{12.76f, 6.3f, 0.02962f, 1.294f},
        {32.f, 56.3f, 68.4f, 71.5f, 86.6f, 92.7f, 120.f, 123.1f, 144.3f, 147.4f,
         150.4f, 168.6f, 171.6f, 174.7f, 177.7f, 183.8f},
        {0.f, 0.f, 157.2f, 92.64f, 157.2f, 84.66f, 147.f, 72.21f, 154.9f, 90.9f,
         29.35f, 154.3f, 90.32f, 175.f, 696.5f, 0.f},
        {0.f, 0.f, 301.8f, 257.9f, 301.8f, 249.6f, 292.6f, 232.9f, 301.2f,
         254.4f, 53.95f, 299.4f, 253.9f, 322.4f, 1264.f, 0.f}}},
      {1523.f,
       {{105.5f, 173.8f}},
       {BandComposition{12.76f, 6.3f, 0.02962f, 1.294f},
        {32.f, 56.3f, 68.4f, 71.5f, 86.6f, 92.7f, 120.f, 123.1f, 144.3f, 147.4f,
         150.4f, 168.6f, 171.6f, 174.7f, 177.7f, 183.8f},
        {0.f, 0.f, 157.2f, 92.64f, 157.2f, 84.66f, 147.f, 72.21f, 154.9f, 90.9f,
         29.35f, 154.3f, 90.32f, 175.f, 696.5f, 0.f},
        {0.f, 0.f, 301.8f, 257.9f, 301.8f, 249.6f, 292.6f, 232.9f, 301.2f,
         254.4f, 53.95f, 299.4f, 253.9f, 322.4f, 1264.f, 0.f}}},
  };
  return description;
}

Synthetic::DetectorLayout Synthetic::makeOpenDataDetectorPixelLayout() {
  return makeLayout(openDataDetectorPixelDescription());
}

Synthetic::BarrelEndcapDescription
Synthetic::genericDetectorPixelDescription() {
  // Read off the built `ActsExamples::GenericDetector` as above. Its sensors
  // sit at the nominal layer radii, unlike the ODD's.
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
