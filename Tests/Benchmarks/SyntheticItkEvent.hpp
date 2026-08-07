// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// Synthetic ITk-like pixel event generator shared by the seeding benchmarks.
///
/// The layout is a coarse stand-in for the ATLAS ITk pixel detector: five
/// barrel cylinders subdivided into eta modules and nine endcap disks per side
/// subdivided into radial rings. Tracks are helices from a smeared luminous
/// region, with a soft-dominated 1/pT spectrum and a displaced component that
/// stands in for secondaries, loopers and combinatorial background.
///
/// The generator is deliberately seeder-agnostic: it produces a plain
/// `Acts::SpacePointContainer` and the layer index each space point belongs to,
/// so it can drive any seeding implementation. Seeder-specific geometry (such
/// as the GBTS layer-connection table) is built from `DetectorLayout` by the
/// individual benchmarks.

#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <cmath>
#include <cstdint>
#include <numbers>
#include <random>
#include <vector>

namespace ActsTests::SyntheticItk {

/// Whether a layer is a cylinder around the beam axis or a disk normal to it.
enum class LayerType { Barrel, Endcap };

/// One logical layer, i.e. the granularity a seeder reasons about. A barrel
/// cylinder is split into eta modules along z, an endcap disk into rings in r.
struct Layer {
  /// Barrel cylinder or endcap disk
  LayerType type{};
  /// Radius for a barrel layer, signed z position for an endcap disk
  float refCoord{};
  /// Lower bound along the extended coordinate: z for barrel, r for endcap
  float minBound{};
  /// Upper bound along the extended coordinate: z for barrel, r for endcap
  float maxBound{};
  /// +1 for the positive endcap, -1 for the negative endcap, 0 for the barrel
  int side{};
  /// Barrel layer index, or disk index counting outwards from the interaction
  /// point within one endcap
  int layer{};
  /// Eta module index for a barrel layer, ring index for an endcap disk
  int module{};
};

/// A physical detector element together with the logical layers it is split
/// into. Track propagation intersects surfaces; the resulting position then
/// selects one of the surface's modules.
struct Surface {
  /// Barrel cylinder or endcap disk
  LayerType type{};
  /// Radius for a barrel cylinder, signed z position for an endcap disk
  float refCoord{};
  /// Indices into DetectorLayout::layers of the modules of this surface
  std::vector<std::uint32_t> modules;
};

/// The full synthetic detector description.
struct DetectorLayout {
  /// All logical layers, in the index space used by the "layerId" column
  std::vector<Layer> layers;
  /// The physical surfaces the layers are grouped into
  std::vector<Surface> surfaces;
};

/// Steering for the synthetic event.
struct EventConfig {
  /// Number of generated tracks. Roughly 6 space points per track, so the
  /// default lands near the 200k space points of an ITk ttbar event at a
  /// pileup of 200.
  std::size_t numTracks = 33000;
  /// Lowest generated transverse momentum in GeV
  float minPt = 0.4f;
  /// Steepness of the 1/pT spectrum; closer to one means a harder tail
  float ptSpectrumSlope = 0.98f;
  /// Pseudorapidity range
  float minEta = -4.f;
  /// Upper end of the pseudorapidity range
  float maxEta = 4.f;
  /// Longitudinal spread of the luminous region in mm
  float beamspotSigmaZ = 45.f;
  /// Fraction of tracks that do not point back to the luminous region
  float displacedFraction = 0.4f;
  /// Half-width in mm of the flat z0 distribution of the displaced component
  float displacedRangeZ = 400.f;
  /// Gaussian position smearing per coordinate in mm
  float positionSmearing = 0.05f;
  /// Solenoid field along z in Tesla
  float bFieldZ = 2.f;
  /// Random seed, so that events are reproducible across runs and seeders
  unsigned int seed = 12345;
};

/// Build the ITk-like pixel layout.
/// @return the detector layout used by the seeding benchmarks
inline DetectorLayout makePixelLayout() {
  // Barrel cylinder radii in mm and their common half-length in z
  constexpr int numBarrel = 5;
  constexpr float barrelRadius[numBarrel] = {34.f, 99.f, 160.f, 228.f, 291.f};
  constexpr float barrelHalfZ = 250.f;
  // The barrel is left unsplit so that the logical layer ids come out as the
  // round 80000, 81000, ... that GBTS special-cases as its innermost layers.
  constexpr int barrelModules = 1;

  // Endcap disk positions in mm and their radial extent
  constexpr int numDisks = 9;
  constexpr float diskZ[numDisks] = {400.f,  600.f,  800.f,  1050.f, 1300.f,
                                     1600.f, 1900.f, 2300.f, 2800.f};
  constexpr float diskRMin = 30.f;
  constexpr float diskRMax = 350.f;
  constexpr int diskRings = 2;

  DetectorLayout layout;

  auto addSurface = [&](LayerType type, float refCoord, int side, int layerIdx,
                        float lo, float hi, int nSplit) {
    Surface surface;
    surface.type = type;
    surface.refCoord = refCoord;
    for (int m = 0; m < nSplit; ++m) {
      Layer layer;
      layer.type = type;
      layer.refCoord = refCoord;
      layer.minBound = lo + (hi - lo) * m / nSplit;
      layer.maxBound = lo + (hi - lo) * (m + 1) / nSplit;
      layer.side = side;
      layer.layer = layerIdx;
      layer.module = m;
      surface.modules.push_back(
          static_cast<std::uint32_t>(layout.layers.size()));
      layout.layers.push_back(layer);
    }
    layout.surfaces.push_back(std::move(surface));
  };

  for (int i = 0; i < numBarrel; ++i) {
    addSurface(LayerType::Barrel, barrelRadius[i], 0, i, -barrelHalfZ,
               barrelHalfZ, barrelModules);
  }
  for (int side : {+1, -1}) {
    for (int j = 0; j < numDisks; ++j) {
      addSurface(LayerType::Endcap, side * diskZ[j], side, j, diskRMin,
                 diskRMax, diskRings);
    }
  }

  return layout;
}

/// Generate a synthetic event.
///
/// The returned container carries the standard X/Y/Z/R/Phi columns plus three
/// dynamic columns: "layerId" with the index into `layout.layers`,
/// "clusterWidth" and "localPositionY" as expected by the GBTS seeder, and
/// "particleId" with the index of the generating track so that seeding
/// efficiency can be evaluated without a truth-matching stage.
///
/// @param layout the detector to generate hits on
/// @param cfg steering for the generated event
/// @return the generated space points
inline Acts::SpacePointContainer generateEvent(const DetectorLayout& layout,
                                               const EventConfig& cfg) {
  constexpr float pi = std::numbers::pi_v<float>;

  struct Hit {
    float x{};
    float y{};
    float z{};
    std::uint32_t layer{};
    std::uint32_t particle{};
  };

  std::mt19937 rng(cfg.seed);
  std::uniform_real_distribution<float> phiDist(-pi, pi);
  std::uniform_real_distribution<float> etaDist(cfg.minEta, cfg.maxEta);
  std::uniform_real_distribution<float> uniform(0.f, 1.f);
  std::normal_distribution<float> z0Dist(0.f, cfg.beamspotSigmaZ);
  std::uniform_real_distribution<float> z0Displaced(-cfg.displacedRangeZ,
                                                    cfg.displacedRangeZ);
  std::normal_distribution<float> smear(0.f, cfg.positionSmearing);

  std::vector<Hit> hits;
  hits.reserve(cfg.numTracks * 8);

  for (std::size_t t = 0; t < cfg.numTracks; ++t) {
    const float pt = cfg.minPt / (1.f - cfg.ptSpectrumSlope * uniform(rng));
    const float charge = uniform(rng) < 0.5f ? -1.f : +1.f;
    const float phi0 = phiDist(rng);
    const float cotTheta = std::sinh(etaDist(rng));
    const bool displaced = uniform(rng) < cfg.displacedFraction;
    const float z0 = displaced ? z0Displaced(rng) : z0Dist(rng);

    // Radius of curvature in mm for a momentum in GeV and a field in Tesla
    const float curvatureRadius = pt * 1000.f / (0.3f * cfg.bFieldZ);

    for (const Surface& surface : layout.surfaces) {
      float r = 0.f;
      float z = 0.f;

      if (surface.type == LayerType::Barrel) {
        r = surface.refCoord;
        const float arg = r / (2.f * curvatureRadius);
        if (arg >= 1.f) {
          // the track curls up before reaching this radius
          continue;
        }
        // arc length in the transverse plane up to this radius
        z = z0 + 2.f * curvatureRadius * std::asin(arg) * cotTheta;
      } else {
        z = surface.refCoord;
        if (cotTheta == 0.f) {
          continue;
        }
        const float arc = (z - z0) / cotTheta;
        if (arc <= 0.f) {
          // the disk is on the far side of the track's direction of flight
          continue;
        }
        const float arg = arc / (2.f * curvatureRadius);
        if (arg >= pi / 2.f) {
          continue;
        }
        r = 2.f * curvatureRadius * std::sin(arg);
      }

      // pick the module of this surface that the intersection falls into
      const float selector = surface.type == LayerType::Barrel ? z : r;
      std::uint32_t layerIndex = 0;
      bool found = false;
      for (const std::uint32_t m : surface.modules) {
        const Layer& layer = layout.layers[m];
        if (selector >= layer.minBound && selector < layer.maxBound) {
          layerIndex = m;
          found = true;
          break;
        }
      }
      if (!found) {
        continue;
      }

      const float arg = r / (2.f * curvatureRadius);
      if (arg >= 1.f) {
        continue;
      }
      const float phi = phi0 + charge * std::asin(arg);

      hits.push_back(Hit{r * std::cos(phi) + smear(rng),
                         r * std::sin(phi) + smear(rng), z + smear(rng),
                         layerIndex, static_cast<std::uint32_t>(t)});
    }
  }

  Acts::SpacePointContainer spacePoints(
      Acts::SpacePointColumns::CopiedFromIndex | Acts::SpacePointColumns::X |
      Acts::SpacePointColumns::Y | Acts::SpacePointColumns::Z |
      Acts::SpacePointColumns::R | Acts::SpacePointColumns::Phi);
  auto layerColumn = spacePoints.createColumn<std::uint32_t>("layerId");
  auto clusterWidthColumn = spacePoints.createColumn<float>("clusterWidth");
  auto localPositionColumn = spacePoints.createColumn<float>("localPositionY");
  auto particleColumn = spacePoints.createColumn<std::uint32_t>("particleId");
  spacePoints.reserve(hits.size());

  for (const Hit& hit : hits) {
    auto sp = spacePoints.createSpacePoint();
    sp.x() = hit.x;
    sp.y() = hit.y;
    sp.z() = hit.z;
    sp.r() = Acts::fastHypot(hit.x, hit.y);
    sp.phi() = std::atan2(hit.y, hit.x);
    sp.copiedFromIndex() = sp.index();
    sp.extra(layerColumn) = hit.layer;
    sp.extra(clusterWidthColumn) = 0.f;
    sp.extra(localPositionColumn) = 0.f;
    sp.extra(particleColumn) = hit.particle;
  }

  return spacePoints;
}

}  // namespace ActsTests::SyntheticItk
