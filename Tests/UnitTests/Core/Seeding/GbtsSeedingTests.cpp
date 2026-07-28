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

#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <sstream>
#include <vector>

// Regression harness for the GBTS seeding chain, which has no other test
// coverage. The expected values are recorded from the current implementation,
// not derived independently: they detect behaviour changes, not physics errors.

namespace Acts::Test {

namespace {

using namespace Acts::UnitLiterals;

constexpr std::size_t kNumLayers = 4;

// GBTS reads the subdetector off the layer id: `id / 10000 == 8` is a barrel
// layer, 80000 the innermost one (extra z0 cuts), ids 1000 apart are adjacent.
constexpr std::array<std::int32_t, kNumLayers> kLayerIds = {80000, 81000, 82000,
                                                            83000};
constexpr std::array<float, kNumLayers> kLayerRadii = {40.f, 80.f, 120.f,
                                                       160.f};

// Half-length in z of every layer, reused as the z0 and RoI z range.
constexpr float kHalfZ = 150.f;

/// Connector table for GbtsLayerConnectionMap::fromStream: `nLinks etaBinWidth`
/// then `lIdx stage src dst height width nEntries` per link, `src` outer and
/// `dst` inner. height = width = 0 leaves out the bin table, which GbtsGeometry
/// recomputes from the layer geometry anyway.
constexpr const char* kConnectorText =
    "3 0.2\n"
    "0 0 81000 80000 0 0 0\n"
    "1 1 82000 81000 0 0 0\n"
    "2 2 83000 82000 0 0 0\n";

std::shared_ptr<Experimental::GbtsGeometry> makeGeometry() {
  std::vector<Experimental::GbtsLayerDescription> layers;
  layers.reserve(kNumLayers);
  for (std::size_t i = 0; i < kNumLayers; ++i) {
    Experimental::GbtsLayerDescription layer;
    layer.id = kLayerIds[i];
    layer.type = Experimental::GbtsLayerType::Barrel;
    // For a barrel layer refCoord is the radius and the bounds are in z.
    layer.refCoord = kLayerRadii[i];
    layer.minBound = -kHalfZ;
    layer.maxBound = kHalfZ;
    layers.push_back(layer);
  }

  std::istringstream stream{kConnectorText};
  return std::make_shared<Experimental::GbtsGeometry>(
      layers, Experimental::GbtsLayerConnectionMap::fromStream(stream, false));
}

/// Straight track from the origin: fixed phi, z = r * tau, tau = cot(theta).
/// Zero curvature and constant tau ratio pass the doublet and triplet cuts.
struct Track {
  float phi{};
  float tau{};
};

/// Evenly spaced tracks over [phiMin, phiMax] x [tauMin, tauMax]. |tau| must
/// stay below about 0.9 to keep every hit inside the layer z bounds.
std::vector<Track> makeTracks(std::size_t nTracks, float phiMin, float phiMax,
                              float tauMin, float tauMax) {
  std::vector<Track> tracks;
  tracks.reserve(nTracks);
  for (std::size_t i = 0; i < nTracks; ++i) {
    const float frac = static_cast<float>(i) / nTracks;
    Track track;
    // half a step offset so a full phi range has no track on the wrap-around
    track.phi = phiMin + (phiMax - phiMin) * (frac + 0.5f / nTracks);
    track.tau = tauMin + (tauMax - tauMin) * frac;
    tracks.push_back(track);
  }
  return tracks;
}

/// Far apart in phi and tau: no cross-track edges, one clean seed per track.
std::vector<Track> makeSparseTracks() {
  return makeTracks(12, -std::numbers::pi_v<float>, std::numbers::pi_v<float>,
                    -0.6f, 0.6f);
}

/// Close enough in phi and tau for cross-track edges, so the graph
/// combinatorics, the tracking filter and clone removal all do real work.
std::vector<Track> makeDenseTracks() {
  return makeTracks(40, -0.05f, 0.05f, -0.05f, 0.05f);
}

/// One hit per layer per track, in the container layout
/// GraphBasedSeedingAlgorithm::makeSpContainer produces.
SpacePointContainer makeSpacePoints(const std::vector<Track>& tracks) {
  SpacePointContainer container(SpacePointColumns::CopiedFromIndex |
                                SpacePointColumns::X | SpacePointColumns::Y |
                                SpacePointColumns::Z | SpacePointColumns::R |
                                SpacePointColumns::Phi);

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
      // the dense layer index, not the GBTS layer id
      sp.extra(layerColumn) = static_cast<std::uint32_t>(layer);
      sp.extra(clusterWidthColumn) = 0.f;
      sp.extra(localPositionColumn) = 0.f;
    }
  }

  return container;
}

SeedContainer runSeeding(const SpacePointContainer& spacePoints) {
  auto geometry = makeGeometry();

  Experimental::GraphBasedTrackSeeder::Config config;
  config.minPt = 1_GeV;
  config.minZ0 = -kHalfZ;
  config.maxZ0 = kHalfZ;
  config.maxOuterRadius = 200.f;
  // the toy setup has no ML lookup table and no cluster widths
  config.useMl = false;

  const Experimental::GraphBasedTrackSeeder seeder(
      Experimental::GraphBasedTrackSeeder::DerivedConfig(config), geometry,
      getDefaultLogger("GbtsTest", Logging::Level::WARNING));

  const Experimental::GbtsTrackingFilter filter(
      Experimental::GbtsTrackingFilter::Config{}, geometry);

  const Experimental::GbtsRoiDescriptor roi(
      0, -4.5, 4.5, 0, -std::numbers::pi, std::numbers::pi, 0, -kHalfZ, kHalfZ);

  const Experimental::GraphBasedTrackSeeder::Options options(2_T);

  const std::vector<bool> isPixelLayer(kNumLayers, true);

  SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);

  seeder.createSeeds(spacePoints, roi, isPixelLayer, kNumLayers, filter,
                     options, seeds);

  return seeds;
}

/// One line per seed, `quality:sp,sp,...`.
std::string formatSeeds(const SeedContainer& seeds) {
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

/// Invariants that hold independently of the recorded values. The container is
/// filled track-major, so `index / kNumLayers` identifies the track.
void checkSeedsAreWellFormed(const SeedContainer& seeds,
                             const SpacePointContainer& spacePoints) {
  for (const auto& seed : seeds) {
    const auto indices = seed.spacePointIndices();
    BOOST_REQUIRE_GE(indices.size(), 3u);

    float previousR = -1.f;
    const std::uint32_t track = indices[0] / kNumLayers;
    for (const auto index : indices) {
      BOOST_REQUIRE_LT(index, spacePoints.size());
      BOOST_CHECK_EQUAL(index / kNumLayers, track);
      const auto sp = spacePoints.at(index);
      BOOST_CHECK_GT(sp.r(), previousR);
      previousR = sp.r();
    }
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(GbtsSeeding)

// Guards the fixture itself: a broken input must not look like a regression.
BOOST_AUTO_TEST_CASE(SyntheticInputIsWellFormed) {
  const std::vector<Track> tracks = makeSparseTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(tracks);

  BOOST_CHECK_EQUAL(spacePoints.size(), tracks.size() * kNumLayers);

  auto layerColumn = spacePoints.column<std::uint32_t>("layerId");
  for (const auto& sp : spacePoints) {
    BOOST_CHECK_LT(sp.extra(layerColumn), kNumLayers);
    BOOST_CHECK_LE(std::abs(sp.z()), kHalfZ);
  }
}

// One seed per track, all four hits, in radial order.
BOOST_AUTO_TEST_CASE(SeedsFromSeparatedTracks) {
  const std::vector<Track> tracks = makeSparseTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(tracks);

  const SeedContainer seeds = runSeeding(spacePoints);

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

// Dense tracks put several candidates per node into the sliding window.
BOOST_AUTO_TEST_CASE(SeedsFromDenseTracks) {
  const std::vector<Track> tracks = makeDenseTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(tracks);

  const SeedContainer seeds = runSeeding(spacePoints);

  // dumped rather than pinned: under heavy branching which candidates survive
  // depends on float rounding and differs between platforms
  BOOST_TEST_MESSAGE("dense seeds:\n" << formatSeeds(seeds));

  BOOST_CHECK_EQUAL(seeds.size(), tracks.size());
  checkSeedsAreWellFormed(seeds, spacePoints);

  // clone removal resolves the branching back to one complete seed per track
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
