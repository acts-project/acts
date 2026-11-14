// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/EventGenerator.hpp"
#include "ActsFatras/Synthetic/SyntheticEvent.hpp"
#include "ActsFatras/Synthetic/detail/Helix.hpp"
#include "ActsFatras/Synthetic/detail/Sampling.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

using namespace ActsFatras::Synthetic;
using ActsFatras::Synthetic::detail::helixRadius;
using ActsFatras::Synthetic::detail::samplePt;
using ActsFatras::Synthetic::detail::sampleRapidity;
using namespace Acts::UnitLiterals;

namespace ActsTests {

namespace {

DetectorLayout makeTestLayout() {
  // silicon at a percent and a half of a radiation length, which is a module
  const Acts::MaterialSlab sensor =
      materialSlab(93.7f, 465.2f, 28.0855f, 14.f, 8.2925e-5f, 0.015f * 93.7f);

  DetectorDescription description;
  description.passives = {
      PassiveSurfaceDescription{.shape = SurfaceShape::Cylinder,
                                .refCoord = 25.f,
                                .maxBound = 800.f,
                                .material = SurfaceMaterial{sensor}}};

  SubsystemDescription pixels;
  pixels.name = "pixel";
  BarrelDescription barrel;
  for (const float radius : {34.f, 99.f, 160.f, 228.f}) {
    barrel.cylinders.push_back(
        CylinderDescription{.radius = radius,
                            .halfLengthZ = 250.f,
                            .material = SurfaceMaterial{sensor}});
  }
  pixels.barrels.push_back(std::move(barrel));
  EndcapDescription endcap;
  for (const float absZ : {400.f, 600.f, 800.f}) {
    endcap.discs.push_back(
        DiscDescription{.absZ = absZ,
                        .rings = {{30.f, 200.f}},
                        .material = SurfaceMaterial{sensor}});
  }
  pixels.endcaps.push_back(std::move(endcap));
  description.subsystems.push_back(std::move(pixels));
  return makeLayout(description);
}

/// A handful of primaries and nothing else, so that the random sequence does
/// not depend on how much material the toggles under test leave behind.
EventConfig makeTestConfig() {
  EventConfig config;
  config.generation.pileup = 1;
  config.generation.chargedPerUnitRapidity = 25.f;
  config.generation.minRapidity = -1.f;
  config.generation.maxRapidity = 1.f;
  config.generation.minPt = 0.5_GeV;
  config.simulation.secondaries.electronRate = 0.f;
  config.simulation.secondaries.nuclearRate = 0.f;
  config.simulation.secondaries.stubRate = 0.f;
  config.simulation.secondaries.decayYield = 0.f;
  config.simulation.measurement.positionSmearing = 0.f;
  config.simulation.material.multipleScattering = false;
  config.simulation.material.energyLoss = false;
  config.simulation.propagation.maxTurns = 0.5f;
  return config;
}

/// The space points of one event, in the order the generator wrote them, which
/// per particle is the order the track met its surfaces.
struct Hit {
  float x{};
  float y{};
  float z{};
  /// Transverse momentum of the particle that left it
  float pt{};
  std::uint32_t particle{};
};

std::vector<Hit> generateHits(const DetectorLayout& layout,
                              const EventConfig& config) {
  const Event event = generateEvent(layout, config);
  auto particles = event.spacePoints.column<std::uint32_t>("particleId");
  std::vector<Hit> hits;
  hits.reserve(event.spacePoints.size());
  for (const auto sp : event.spacePoints) {
    const std::uint32_t particle = sp.extra(particles);
    hits.push_back(
        Hit{sp.x(), sp.y(), sp.z(), event.particles[particle].pt, particle});
  }
  return hits;
}

/// Radius of the circle through three points, which is what a three-point seed
/// reads a transverse momentum off.
float circleRadius(const Hit& a, const Hit& b, const Hit& c) {
  const float ab = std::hypot(b.x - a.x, b.y - a.y);
  const float bc = std::hypot(c.x - b.x, c.y - b.y);
  const float ca = std::hypot(a.x - c.x, a.y - c.y);
  const float cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  return ab * bc * ca / (2.f * std::abs(cross));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(SyntheticEventGeneratorSuite)

/// With the material effects switched off a hit sits at the nominal position of
/// its surface in the unsmeared directions.
BOOST_AUTO_TEST_CASE(NoMaterialEffectsIsAnExactHelix) {
  const DetectorLayout layout = makeTestLayout();
  const Event event = generateEvent(layout, makeTestConfig());
  BOOST_REQUIRE(!event.spacePoints.empty());

  auto layers = event.spacePoints.column<std::uint32_t>("layerId");
  for (const auto sp : event.spacePoints) {
    const DetectorLayer& layer = layout.layers[sp.extra(layers)];
    const float measured =
        layer.shape == SurfaceShape::Cylinder ? sp.r() : sp.z();
    BOOST_CHECK_SMALL(measured - layer.refCoord, 1e-3f);
  }
}

/// Energy loss leaves the first crossing of a track alone and moves the later
/// ones by more and more, the lever arm growing with the path since the loss.
BOOST_AUTO_TEST_CASE(EnergyLossGrowsAlongTheTrack) {
  const DetectorLayout layout = makeTestLayout();

  EventConfig config = makeTestConfig();
  const std::vector<Hit> exact = generateHits(layout, config);
  config.simulation.material.energyLoss = true;
  const std::vector<Hit> lossy = generateHits(layout, config);

  // no random number is drawn for the loss, so the two events hold the same
  // tracks and the same hits, only displaced
  BOOST_REQUIRE_EQUAL(exact.size(), lossy.size());
  BOOST_REQUIRE(!exact.empty());

  float largest = 0.f;
  std::uint32_t previous = std::numeric_limits<std::uint32_t>::max();
  float previousShift = 0.f;
  for (std::size_t h = 0; h < exact.size(); ++h) {
    BOOST_REQUIRE_EQUAL(exact[h].particle, lossy[h].particle);
    const float shift =
        std::hypot(lossy[h].x - exact[h].x, lossy[h].y - exact[h].y);
    if (exact[h].particle != previous) {
      // the first crossing has no material upstream of it
      BOOST_CHECK_SMALL(shift, 1e-4f);
      previous = exact[h].particle;
    } else {
      BOOST_CHECK_GE(shift, previousShift);
    }
    previousShift = shift;
    largest = std::max(largest, shift);
  }
  BOOST_CHECK_GT(largest, 1_um);
}

/// The displacement follows how much is lost: more material takes more energy,
/// and the most probable loss stays below the mean one, the two differing by
/// the tail of delta rays.
BOOST_AUTO_TEST_CASE(DisplacementFollowsWhatIsLost) {
  const DetectorLayout layout = makeTestLayout();

  EventConfig config = makeTestConfig();
  const std::vector<Hit> exact = generateHits(layout, config);
  config.simulation.material.energyLoss = true;

  // How far the material moved the hits of the whole event. No random number
  // is drawn for the loss, so every event here holds the same hits, displaced.
  const auto totalShift = [&](const EventConfig& lossyConfig) {
    const std::vector<Hit> lossy = generateHits(layout, lossyConfig);
    BOOST_REQUIRE_EQUAL(exact.size(), lossy.size());
    float shift = 0.f;
    for (std::size_t h = 0; h < exact.size(); ++h) {
      shift += std::hypot(lossy[h].x - exact[h].x, lossy[h].y - exact[h].y);
    }
    return shift;
  };

  float previous = 0.f;
  for (const float materialScale : {1.f, 2.f, 4.f}) {
    config.simulation.material.scale = materialScale;
    const float shift = totalShift(config);
    BOOST_CHECK_GT(shift, previous);
    previous = shift;
  }

  config.simulation.material.scale = 1.f;
  config.simulation.material.energyLossModel = EnergyLossModel::Mode;
  const float mode = totalShift(config);
  config.simulation.material.energyLossModel = EnergyLossModel::Mean;
  BOOST_CHECK_GT(mode, 0.f);
  BOOST_CHECK_LT(mode, totalShift(config));
}

/// Scattering turns a track either way, so half the triplets read a momentum
/// above the truth and half below; energy loss only ever takes momentum away.
BOOST_AUTO_TEST_CASE(ScatteringIsSymmetricInCurvatureAndEnergyLossIsNot) {
  const DetectorLayout layout = makeTestLayout();

  EventConfig config = makeTestConfig();
  config.generation.pileup = 20;

  for (const bool scattering : {true, false}) {
    EventConfig moved = config;
    moved.simulation.material.multipleScattering = scattering;
    moved.simulation.material.energyLoss = !scattering;
    const std::vector<Hit> hits = generateHits(layout, moved);

    std::size_t tighter = 0;
    std::size_t compared = 0;
    for (std::size_t h = 2; h < hits.size(); ++h) {
      if (hits[h].particle != hits[h - 2].particle) {
        continue;
      }
      ++compared;
      tighter += static_cast<std::size_t>(
          circleRadius(hits[h - 2], hits[h - 1], hits[h]) <
          helixRadius(hits[h].pt, moved.bFieldZ));
    }
    BOOST_REQUIRE_GT(compared, 100u);

    const float share =
        static_cast<float>(tighter) / static_cast<float>(compared);
    if (scattering) {
      BOOST_CHECK_GT(share, 0.4f);
      BOOST_CHECK_LT(share, 0.6f);
    } else {
      BOOST_CHECK_GT(share, 0.99f);
    }
  }
}

/// The same seed gives the same event, which is what makes a benchmark
/// comparable across runs and seeders.
BOOST_AUTO_TEST_CASE(SameSeedSameEvent) {
  const DetectorLayout layout = makeTestLayout();
  EventConfig config = makeTestConfig();
  config.simulation.measurement.positionSmearing = 15_um;
  config.simulation.secondaries.electronRate = 4.941f;
  config.simulation.secondaries.nuclearRate = 11.545f;
  config.simulation.material.multipleScattering = true;
  config.simulation.material.energyLoss = true;

  const std::vector<Hit> first = generateHits(layout, config);
  const std::vector<Hit> second = generateHits(layout, config);
  BOOST_REQUIRE(!first.empty());
  BOOST_REQUIRE_EQUAL(first.size(), second.size());
  for (std::size_t h = 0; h < first.size(); ++h) {
    BOOST_CHECK_EQUAL(first[h].x, second[h].x);
    BOOST_CHECK_EQUAL(first[h].y, second[h].y);
    BOOST_CHECK_EQUAL(first[h].z, second[h].z);
  }
}

/// A particle the core library has no mass for cannot be generated, and one
/// told to lose everything it has cannot be propagated.
BOOST_AUTO_TEST_CASE(RejectsAConfigItCannotHonour) {
  const DetectorLayout layout = makeTestLayout();

  EventConfig noMass = makeTestConfig();
  noMass.particlePdg = Acts::PdgParticle::eInvalid;
  BOOST_CHECK_THROW(EventGenerator(layout, noMass), std::invalid_argument);

  EventConfig everything = makeTestConfig();
  everything.simulation.material.maxEnergyLossFraction = 1.f;
  BOOST_CHECK_THROW(EventGenerator(layout, everything), std::invalid_argument);
}

/// A barrel whose modules always overlap measures every crossing twice, the
/// second space point a stagger further out in radius on the same layer.
BOOST_AUTO_TEST_CASE(ModuleOverlapStaggersASecondSpacePoint) {
  constexpr float offset = 8.f;
  const std::vector<float> radii = {34.f, 99.f, 160.f, 228.f};

  BarrelDescription barrel;
  for (const float radius : radii) {
    // long enough that no track leaves through an end: a stagger outboard in
    // radius is met further along z, and one that falls off has no partner
    barrel.cylinders.push_back(CylinderDescription{.radius = radius,
                                                   .halfLengthZ = 600.f,
                                                   .overlapProbability = 1.f,
                                                   .overlapOffset = offset});
  }
  DetectorDescription description;
  description.subsystems = {
      SubsystemDescription{.name = "pixel", .barrels = {std::move(barrel)}}};
  const DetectorLayout layout = makeLayout(description);

  EventConfig config = makeTestConfig();
  const std::vector<Hit> overlapping = generateHits(layout, config);
  config.simulation.measurement.overlapScale = 0.f;
  const std::vector<Hit> plain = generateHits(layout, config);

  BOOST_REQUIRE(!plain.empty());
  BOOST_CHECK_EQUAL(overlapping.size(), 2 * plain.size());

  // every extra one on a module a stagger outboard of a nominal radius
  std::size_t staggered = 0;
  for (const Hit& hit : overlapping) {
    const float r = std::hypot(hit.x, hit.y);
    for (const float radius : radii) {
      if (std::abs(r - radius - offset) < 1e-2f) {
        ++staggered;
      }
    }
  }
  BOOST_CHECK_EQUAL(staggered, plain.size());
}

/// A disc staggers along z instead, and away from the interaction point on both
/// sides of it rather than towards it on one.
BOOST_AUTO_TEST_CASE(ModuleOverlapStaggersADiscAwayFromTheOrigin) {
  constexpr float absZ = 400.f;
  constexpr float offset = 5.f;

  // discs alone, so that every space point below is one of theirs
  DetectorDescription description;
  description.subsystems = {SubsystemDescription{
      .name = "pixel",
      .endcaps = {EndcapDescription{
          .discs = {DiscDescription{.absZ = absZ,
                                    .rings = {{30.f, 200.f}},
                                    .overlapProbability = 1.f,
                                    .overlapOffset = offset}}}}}};
  const DetectorLayout layout = makeLayout(description);

  EventConfig config = makeTestConfig();
  config.generation.minRapidity = -3.f;
  config.generation.maxRapidity = 3.f;
  const std::vector<Hit> hits = generateHits(layout, config);

  std::size_t staggered = 0;
  for (const Hit& hit : hits) {
    if (std::abs(std::abs(hit.z) - absZ) < 1e-2f) {
      continue;
    }
    BOOST_CHECK_SMALL(std::abs(hit.z) - absZ - offset, 1e-2f);
    ++staggered;
  }
  BOOST_CHECK_GT(staggered, std::size_t{0});
}

/// `samplePt` inverts the survival function of the spectrum it documents, so
/// the survival at what it draws has to be the uniform it was given.
BOOST_AUTO_TEST_CASE(SamplePtInvertsItsSpectrum) {
  const std::array<std::array<float, 3>, 4> cases = {{{0.02f, 0.948f, 6.57f},
                                                      {0.1f, 1.404f, 5.00f},
                                                      {0.02f, 1.316f, 4.68f},
                                                      {0.05f, 2.0f, 9.00f}}};
  for (const auto& [minPt, scale, exponent] : cases) {
    const auto survival = [&](double pt) {
      const double v = 1. + pt / scale;
      return std::pow(v, 2. - exponent) / (exponent - 2.) -
             std::pow(v, 1. - exponent) / (exponent - 1.);
    };
    const double total = survival(minPt);

    float previous = 0.f;
    for (int i = 1; i < 2000; ++i) {
      const float u = static_cast<float>(i) / 2000.f;
      const float pt = samplePt(u, minPt, scale, exponent);
      // never below the threshold, and rising with the uniform
      BOOST_CHECK_GE(pt, minPt);
      BOOST_CHECK_GE(pt, previous);
      previous = pt;
      BOOST_CHECK_SMALL(survival(pt) / total - (1. - u), 1e-5);
    }
  }
}

/// A plateau with no edge is the flat draw it replaces, down to the random
/// stream: an event generated without one has to be unchanged.
BOOST_AUTO_TEST_CASE(SampleRapidityWithoutAnEdgeIsFlat) {
  std::mt19937 flat(4242);
  std::mt19937 edged(4242);
  for (int i = 0; i < 1000; ++i) {
    BOOST_CHECK_EQUAL(detail::sampleUniform(flat, -4.3f, 4.3f),
                      sampleRapidity(edged, -4.3f, 4.3f, 4.3f, 0.f));
  }
}

/// With one, the density it draws from is the Fermi factor: flat centrally and
/// half of that at the edge.
BOOST_AUTO_TEST_CASE(SampleRapidityFollowsItsEdge) {
  constexpr float edge = 3.f;
  constexpr float width = 0.4f;
  constexpr float span = 5.f;
  constexpr std::size_t bins = 20;
  constexpr int samples = 2000000;

  std::mt19937 rng(7);
  std::array<std::size_t, bins> counts = {};
  for (int i = 0; i < samples; ++i) {
    const float rapidity = sampleRapidity(rng, -span, span, edge, width);
    BOOST_REQUIRE_GE(rapidity, -span);
    BOOST_REQUIRE_LT(rapidity, span);
    ++counts[static_cast<std::size_t>(std::abs(rapidity) / span * bins)];
  }

  // the innermost bin is deep inside the plateau, so it normalises the rest
  const auto plateau = static_cast<float>(counts[0]);
  for (std::size_t b = 0; b < bins; ++b) {
    const float centre = (static_cast<float>(b) + 0.5f) * span / bins;
    const float expected = 1.f / (1.f + std::exp((centre - edge) / width));
    // a plateau bin holds a hundred thousand entries, so this is several
    // standard deviations either side of the value at the bin centre
    BOOST_CHECK_SMALL(static_cast<float>(counts[b]) / plateau - expected,
                      0.012f);
  }
}
BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
