// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsLayerConnection.hpp"
#include "Acts/Seeding/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"

#include <cmath>
#include <memory>
#include <numbers>
#include <sstream>
#include <vector>

// Regression harness for the GBTS seeding chain. GBTS has no other test
// coverage - it is not exercised by physmon either - so this pins the seed
// output of a fully synthetic setup: four pixel barrel layers and a set of
// straight tracks from the origin.
//
// The expected values below are recorded from the current implementation;
// they are not derived independently. Their purpose is to detect behaviour
// changes during refactoring, not to validate GBTS physics.

namespace Acts::Test {

namespace {

using namespace Acts::UnitLiterals;

// GBTS layer ids. The algorithm derives the subdetector from the id:
// `id / 10000 == 8` marks a pixel barrel layer, and id 80000 additionally
// selects the innermost-layer code path (z0 histogram bitmask).
constexpr std::int32_t kLayerIds[] = {80000, 81000, 82000, 83000};
constexpr float kLayerRadii[] = {40.f, 80.f, 120.f, 160.f};
constexpr std::size_t kNumLayers = 4;

// Half-length in z of every barrel layer.
constexpr float kHalfZ = 150.f;

/// Connector description consumed by GbtsLayerConnectionMap::fromStream.
///
/// Format: `nLinks etaBinWidth`, then one line per link with
/// `lIdx stage src dst height width nEntries` followed by `height * width`
/// integers. We use height = width = 0 so no bin table is read from the
/// stream - GbtsGeometry recomputes the tables from the layer geometry anyway.
///
/// `src` is the outer layer (n2), `dst` the inner one (n1).
std::string connectorText() {
  std::ostringstream os;
  os << "3 0.2\n";
  os << "0 0 81000 80000 0 0 0\n";
  os << "1 1 82000 81000 0 0 0\n";
  os << "2 2 83000 82000 0 0 0\n";
  return os.str();
}

std::shared_ptr<Acts::Experimental::GbtsGeometry> makeGeometry() {
  std::vector<Acts::Experimental::GbtsLayerDescription> layers;
  layers.reserve(kNumLayers);
  for (std::size_t i = 0; i < kNumLayers; ++i) {
    Acts::Experimental::GbtsLayerDescription layer;
    layer.id = kLayerIds[i];
    layer.type = Acts::Experimental::GbtsLayerType::Barrel;
    // For a barrel layer refCoord is the radius and the bounds are in z.
    layer.refCoord = kLayerRadii[i];
    layer.minBound = -kHalfZ;
    layer.maxBound = kHalfZ;
    layers.push_back(layer);
  }

  std::istringstream stream(connectorText());
  const Acts::Experimental::GbtsLayerConnectionMap connections =
      Acts::Experimental::GbtsLayerConnectionMap::fromStream(stream, false);

  return std::make_shared<Acts::Experimental::GbtsGeometry>(layers,
                                                            connections);
}

/// A straight track from the origin: constant phi, and z = r * tau with
/// tau = cot(theta). Such a track has zero curvature and a perfectly constant
/// tau ratio, so it passes the doublet and triplet cuts by construction.
struct Track {
  float phi{};
  float tau{};
};

/// Build the space point container GBTS consumes, one hit per layer per track.
/// The layout mirrors what GraphBasedSeedingAlgorithm::makeSpContainer
/// produces.
Acts::SpacePointContainer makeSpacePoints(const std::vector<Track>& tracks) {
  Acts::SpacePointContainer container(
      Acts::SpacePointColumns::CopiedFromIndex | Acts::SpacePointColumns::X |
      Acts::SpacePointColumns::Y | Acts::SpacePointColumns::Z |
      Acts::SpacePointColumns::R | Acts::SpacePointColumns::Phi);

  auto layerColumn = container.createColumn<std::uint32_t>("layerId");
  auto clusterWidthColumn = container.createColumn<float>("clusterWidth");
  auto localPositionColumn = container.createColumn<float>("localPositionY");

  container.reserve(tracks.size() * kNumLayers);

  for (const Track& track : tracks) {
    for (std::size_t layer = 0; layer < kNumLayers; ++layer) {
      const float r = kLayerRadii[layer];

      auto sp = container.createSpacePoint();
      sp.x() = r * std::cos(track.phi);
      sp.y() = r * std::sin(track.phi);
      sp.z() = r * track.tau;
      sp.r() = r;
      sp.phi() = track.phi;
      sp.copiedFromIndex() = sp.index();
      // The layer column carries the dense layer index, not the GBTS id.
      sp.extra(layerColumn) = static_cast<std::uint32_t>(layer);
      sp.extra(clusterWidthColumn) = 0.f;
      sp.extra(localPositionColumn) = 0.f;
    }
  }

  return container;
}

/// Evenly spaced tracks over [phiMin, phiMax] x [tauMin, tauMax]. The tau range
/// must stay within about +-0.9 so that every hit falls inside the layer z
/// bounds.
std::vector<Track> makeTracks(std::size_t nTracks, float phiMin, float phiMax,
                              float tauMin, float tauMax) {
  std::vector<Track> tracks;
  tracks.reserve(nTracks);
  for (std::size_t i = 0; i < nTracks; ++i) {
    const float frac = static_cast<float>(i) / nTracks;
    Track track;
    track.phi = phiMin + (phiMax - phiMin) * (frac + 0.5f / nTracks);
    track.tau = tauMin + (tauMax - tauMin) * frac;
    tracks.push_back(track);
  }
  return tracks;
}

/// Well separated tracks spanning the full phi range and a wide tau range.
/// Neighbouring tracks are far apart in both, so each track yields exactly one
/// clean seed.
std::vector<Track> makeSparseTracks() {
  return makeTracks(12, -std::numbers::pi_v<float>, std::numbers::pi_v<float>,
                    -0.6f, 0.6f);
}

/// Tracks packed closely in phi *and* tau. Neighbouring tracks fall inside each
/// other's phi sliding window and survive the tau-ratio cut, so cross-track
/// edges are created and the graph combinatorics, the tracking filter and clone
/// removal all do real work - none of which the sparse case reaches.
std::vector<Track> makeDenseTracks() {
  return makeTracks(40, -0.05f, 0.05f, -0.05f, 0.05f);
}

/// Render the seed collection as a stable, comparable string: one line per
/// seed, `quality:sp,sp,...`. This is what the expectations below compare.
std::string formatSeeds(const Acts::SeedContainer& seeds) {
  std::ostringstream os;
  for (const auto& seed : seeds) {
    os << seed.quality() << ":";
    for (const auto index : seed.spacePointIndices()) {
      os << index << ",";
    }
    os << "\n";
  }
  return os.str();
}

Acts::SeedContainer runSeeding(const Acts::SpacePointContainer& spacePoints) {
  auto geometry = makeGeometry();

  Acts::Experimental::GraphBasedTrackSeeder::Config config;
  config.minPt = 1_GeV;
  config.minZ0 = -150.f;
  config.maxZ0 = 150.f;
  config.maxOuterRadius = 200.f;
  // The toy geometry has no ML lookup table and no cluster widths.
  config.useMl = false;

  const Acts::Experimental::GraphBasedTrackSeeder::DerivedConfig derivedConfig(
      config);

  const Acts::Experimental::GraphBasedTrackSeeder seeder(
      derivedConfig, geometry,
      Acts::getDefaultLogger("GbtsTest", Acts::Logging::Level::WARNING));

  const Acts::Experimental::GbtsTrackingFilter filter(
      Acts::Experimental::GbtsTrackingFilter::Config{}, geometry);

  const Acts::Experimental::GbtsRoiDescriptor roi(
      0, -4.5, 4.5, 0, -std::numbers::pi, std::numbers::pi, 0, -150., 150.);

  const Acts::Experimental::GraphBasedTrackSeeder::Options options(2_T);

  // All layers in this setup are pixel layers.
  const std::vector<bool> isPixelLayer(kNumLayers, true);

  Acts::SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);

  seeder.createSeeds(spacePoints, roi, isPixelLayer, filter, options, seeds);

  return seeds;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(GbtsSeeding)

// Sanity check on the synthetic input rather than on GBTS itself: catches a
// broken fixture before it can be mistaken for a seeding regression.
BOOST_AUTO_TEST_CASE(SyntheticInputIsWellFormed) {
  const std::vector<Track> tracks = makeSparseTracks();
  const Acts::SpacePointContainer spacePoints = makeSpacePoints(tracks);

  BOOST_CHECK_EQUAL(spacePoints.size(), tracks.size() * kNumLayers);

  auto layerColumn = spacePoints.column<std::uint32_t>("layerId");
  for (const auto& sp : spacePoints) {
    BOOST_CHECK_LT(sp.extra(layerColumn), kNumLayers);
    BOOST_CHECK_LE(std::abs(sp.z()), kHalfZ);
  }
}

/// Structural invariants that must hold whatever the recorded values are: every
/// seed references in-range space points, ordered innermost first, all from a
/// single track. The container is filled track-major, so `index / kNumLayers`
/// identifies the track.
void checkSeedsAreWellFormed(const Acts::SeedContainer& seeds,
                             const Acts::SpacePointContainer& spacePoints) {
  auto layerColumn = spacePoints.column<std::uint32_t>("layerId");
  for (const auto& seed : seeds) {
    const auto indices = seed.spacePointIndices();
    BOOST_REQUIRE_GE(indices.size(), 3u);

    float previousR = -1.f;
    const std::uint32_t track = indices[0] / kNumLayers;
    for (const auto index : indices) {
      BOOST_REQUIRE_LT(index, spacePoints.size());
      const auto sp = spacePoints.at(index);
      BOOST_CHECK_GT(sp.r(), previousR);
      previousR = sp.r();
      BOOST_CHECK_EQUAL(index / kNumLayers, track);
      BOOST_CHECK_LT(sp.extra(layerColumn), kNumLayers);
    }
  }
}

// Well separated tracks: one seed per track, all four hits, in radial order.
// The expected values are recorded from the current implementation.
BOOST_AUTO_TEST_CASE(SeedsFromSeparatedTracks) {
  const std::vector<Track> tracks = makeSparseTracks();
  const Acts::SpacePointContainer spacePoints = makeSpacePoints(tracks);

  const Acts::SeedContainer seeds = runSeeding(spacePoints);

  BOOST_CHECK_EQUAL(seeds.size(), tracks.size());
  checkSeedsAreWellFormed(seeds, spacePoints);

  const std::string expected =
      "-10.5:0,1,2,3,\n"
      "-10.5:4,5,6,7,\n"
      "-10.5:8,9,10,11,\n"
      "-10.5:12,13,14,15,\n"
      "-10.5:16,17,18,19,\n"
      "-10.5:20,21,22,23,\n"
      "-10.5:24,25,26,27,\n"
      "-10.5:28,29,30,31,\n"
      "-10.5:32,33,34,35,\n"
      "-10.5:36,37,38,39,\n"
      "-10.5:40,41,42,43,\n"
      "-10.5:44,45,46,47,\n";
  BOOST_CHECK_EQUAL(formatSeeds(seeds), expected);
}

// Tracks packed closely enough in phi and tau that the sliding window sees
// several candidates per node, so the graph, the tracking filter and clone
// removal all do real work.
BOOST_AUTO_TEST_CASE(SeedsFromDenseTracks) {
  const std::vector<Track> tracks = makeDenseTracks();
  const Acts::SpacePointContainer spacePoints = makeSpacePoints(tracks);

  const Acts::SeedContainer seeds = runSeeding(spacePoints);

  BOOST_TEST_MESSAGE("dense scenario produced " << seeds.size() << " seeds");
  BOOST_TEST_MESSAGE("\n" << formatSeeds(seeds));

  // Deliberately not pinning the full seed list here. Under heavy branching
  // which candidates survive depends on float rounding, which differs between
  // compilers and architectures; a hard-coded expected string would fail on
  // other CI platforms for reasons unrelated to the code. What is stable is
  // that clone removal resolves the branching back to exactly one complete
  // seed per track. The message dump above lets a refactor be diffed by hand
  // on a fixed machine.
  BOOST_CHECK_EQUAL(seeds.size(), tracks.size());
  checkSeedsAreWellFormed(seeds, spacePoints);

  // Each track is seeded exactly once, with all four of its hits.
  std::vector<std::size_t> seedsPerTrack(tracks.size(), 0);
  for (const auto& seed : seeds) {
    const auto indices = seed.spacePointIndices();
    BOOST_CHECK_EQUAL(indices.size(), kNumLayers);
    seedsPerTrack.at(indices[0] / kNumLayers) += 1;
  }
  for (const std::size_t count : seedsPerTrack) {
    BOOST_CHECK_EQUAL(count, 1u);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
