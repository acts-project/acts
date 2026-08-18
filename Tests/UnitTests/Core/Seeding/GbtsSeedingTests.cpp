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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <numbers>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Regression harness for the GBTS seeding chain, which has no other test
// coverage. The expected values are recorded from the current implementation,
// not derived independently: they detect behaviour changes, not physics errors.

namespace Acts::Test {

namespace {

using namespace Acts::UnitLiterals;
using Experimental::GbtsLayerType;

/// Half-length in z of every barrel layer, reused as the z0 and RoI z range.
constexpr float kBarrelHalfZ = 150.f;

/// Radial extent of every endcap disc.
constexpr float kDiscMinR = 30.f;
constexpr float kDiscMaxR = 220.f;

/// Eta bin width declared by the connector table, as in createLinkingScheme.py.
constexpr float kEtaBinWidth = 0.2f;

/// One layer of the toy detector. GBTS reads the subdetector off the layer id:
/// 8xxxx is barrel, 9xxxx the positive and 7xxxx the negative endcap. 80000 is
/// the innermost barrel layer (extra z0 cuts), barrel ids 1000 apart are
/// adjacent.
struct LayerSpec {
  std::int32_t id{};
  GbtsLayerType type{};
  /// r for a barrel layer, z for an endcap disc.
  float refCoord{};
  /// z range for a barrel layer, r range for an endcap disc.
  float minBound{};
  float maxBound{};
};

constexpr LayerSpec barrelLayer(std::int32_t id, float radius) {
  return {id, GbtsLayerType::Barrel, radius, -kBarrelHalfZ, kBarrelHalfZ};
}

constexpr LayerSpec discLayer(std::int32_t id, float z) {
  return {id, GbtsLayerType::Endcap, z, kDiscMinR, kDiscMaxR};
}

/// Layers plus the connector links between them, outer (src) to inner (dst).
struct ToyDetector {
  std::vector<LayerSpec> layers;
  std::vector<std::pair<std::int32_t, std::int32_t>> links;
  /// Outer radius bound used by the doublet rz filter.
  float maxOuterRadius{};
};

ToyDetector barrelDetector() {
  return {{barrelLayer(80000, 40.f), barrelLayer(81000, 80.f),
           barrelLayer(82000, 120.f), barrelLayer(83000, 160.f)},
          {{81000, 80000}, {82000, 81000}, {83000, 82000}},
          200.f};
}

/// Barrel plus both endcaps. The discs are spaced so that every track in the
/// tau range of makeForwardTracks crosses at least four layers.
ToyDetector forwardDetector() {
  constexpr std::array<float, 4> discZ = {200.f, 260.f, 320.f, 380.f};

  ToyDetector detector = barrelDetector();
  detector.maxOuterRadius = 300.f;

  for (const std::int32_t idBase : {90000, 70000}) {
    const float sign = idBase == 90000 ? 1.f : -1.f;
    std::int32_t previousId = 0;

    for (std::size_t i = 0; i < discZ.size(); ++i) {
      const std::int32_t id = idBase + 1000 * static_cast<std::int32_t>(i);
      detector.layers.push_back(discLayer(id, sign * discZ[i]));

      if (i == 0) {
        // a track reaching the first disc still has hits in one of these
        // barrel layers, never in 83000
        for (const std::int32_t innerId : {80000, 81000, 82000}) {
          detector.links.emplace_back(id, innerId);
        }
      } else {
        detector.links.emplace_back(id, previousId);
      }
      previousId = id;
    }
  }

  return detector;
}

/// Connector table for GbtsLayerConnectionMap::fromStream: `nLinks etaBinWidth`
/// then `lIdx stage src dst height width nEntries` per link. height = width = 0
/// leaves out the bin table and the stage column does not fix the processing
/// order: GbtsGeometry rederives both from the layer geometry.
std::string makeConnectorText(const ToyDetector& detector) {
  std::ostringstream os;
  os << detector.links.size() << " " << kEtaBinWidth << "\n";
  for (std::size_t i = 0; i < detector.links.size(); ++i) {
    const auto& [src, dst] = detector.links[i];
    os << i << " " << i << " " << src << " " << dst << " 0 0 0\n";
  }
  return os.str();
}

std::shared_ptr<Experimental::GbtsGeometry> makeGeometry(
    const ToyDetector& detector) {
  std::vector<Experimental::GbtsLayerDescription> layers;
  layers.reserve(detector.layers.size());
  for (const LayerSpec& spec : detector.layers) {
    Experimental::GbtsLayerDescription layer;
    layer.id = spec.id;
    layer.type = spec.type;
    layer.refCoord = spec.refCoord;
    layer.minBound = spec.minBound;
    layer.maxBound = spec.maxBound;
    layers.push_back(layer);
  }

  const std::string connectorText = makeConnectorText(detector);
  std::istringstream stream{connectorText};
  return std::make_shared<Experimental::GbtsGeometry>(
      layers, Experimental::GbtsLayerConnectionMap::fromStream(stream, false));
}

/// Straight track from the origin: fixed phi, z = r * tau, tau = cot(theta).
/// Zero curvature and constant tau ratio pass the doublet and triplet cuts.
struct Track {
  float phi{};
  float tau{};
};

/// Radius and z where the track crosses the layer, empty if it misses it.
std::optional<std::pair<float, float>> intersect(const LayerSpec& layer,
                                                 const Track& track) {
  if (layer.type == GbtsLayerType::Barrel) {
    const float z = layer.refCoord * track.tau;
    if (z < layer.minBound || z > layer.maxBound) {
      return std::nullopt;
    }
    return std::make_pair(layer.refCoord, z);
  }

  // a disc is only reachable from the side it sits on
  if (track.tau == 0.f || (layer.refCoord > 0.f) != (track.tau > 0.f)) {
    return std::nullopt;
  }
  const float r = layer.refCoord / track.tau;
  if (r < layer.minBound || r > layer.maxBound) {
    return std::nullopt;
  }
  return std::make_pair(r, layer.refCoord);
}

/// Evenly spaced tracks over [phiMin, phiMax] x [tauMin, tauMax].
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
/// |tau| below 150/160 keeps every hit inside the barrel z bounds.
std::vector<Track> makeSparseTracks() {
  return makeTracks(12, -std::numbers::pi_v<float>, std::numbers::pi_v<float>,
                    -0.6f, 0.6f);
}

/// Close enough in phi and tau for cross-track edges, so the graph
/// combinatorics, the tracking filter and clone removal all do real work.
std::vector<Track> makeDenseTracks() {
  return makeTracks(40, -0.05f, 0.05f, -0.05f, 0.05f);
}

/// Well separated tracks over both endcaps. The tau range starts in the
/// barrel-endcap transition, where a seed mixes layer types, and ends where the
/// barrel is missed entirely.
std::vector<Track> makeForwardTracks() {
  constexpr float pi = std::numbers::pi_v<float>;
  std::vector<Track> tracks = makeTracks(8, -pi, 0.f, 1.05f, 6.3f);
  const std::vector<Track> negative = makeTracks(8, 0.f, pi, -1.05f, -6.3f);
  tracks.insert(tracks.end(), negative.begin(), negative.end());
  return tracks;
}

/// Dense forward tracks, all crossing the innermost barrel layer and all four
/// discs of the positive endcap.
std::vector<Track> makeDenseForwardTracks() {
  return makeTracks(20, -0.05f, 0.05f, 1.9f, 2.1f);
}

/// One hit per crossed layer per track, in the container layout
/// GraphBasedSeedingAlgorithm::makeSpContainer produces. `trackId` is test-only
/// bookkeeping: with endcaps the number of hits varies per track.
SpacePointContainer makeSpacePoints(const ToyDetector& detector,
                                    const std::vector<Track>& tracks) {
  SpacePointContainer container(SpacePointColumns::CopiedFromIndex |
                                SpacePointColumns::X | SpacePointColumns::Y |
                                SpacePointColumns::Z | SpacePointColumns::R |
                                SpacePointColumns::Phi);

  auto layerColumn = container.createColumn<std::uint32_t>("layerId");
  auto clusterWidthColumn = container.createColumn<float>("clusterWidth");
  auto localPositionColumn = container.createColumn<float>("localPositionY");
  auto trackColumn = container.createColumn<std::uint32_t>("trackId");

  container.reserve(tracks.size() * detector.layers.size());

  for (std::size_t track = 0; track < tracks.size(); ++track) {
    for (std::size_t layer = 0; layer < detector.layers.size(); ++layer) {
      const auto crossing = intersect(detector.layers[layer], tracks[track]);
      if (!crossing.has_value()) {
        continue;
      }
      const auto [r, z] = *crossing;

      auto sp = container.createSpacePoint();
      sp.x() = r * std::cos(tracks[track].phi);
      sp.y() = r * std::sin(tracks[track].phi);
      sp.z() = z;
      sp.r() = r;
      sp.phi() = tracks[track].phi;
      sp.copiedFromIndex() = sp.index();
      // the dense layer index, not the GBTS layer id
      sp.extra(layerColumn) = static_cast<std::uint32_t>(layer);
      sp.extra(clusterWidthColumn) = 0.f;
      sp.extra(localPositionColumn) = 0.f;
      sp.extra(trackColumn) = static_cast<std::uint32_t>(track);
    }
  }

  return container;
}

SeedContainer runSeeding(const ToyDetector& detector,
                         const SpacePointContainer& spacePoints) {
  auto geometry = makeGeometry(detector);

  const auto numLayers = static_cast<std::uint32_t>(detector.layers.size());

  Experimental::GraphBasedTrackSeeder::Config config;
  config.minPt = 1_GeV;
  config.minZ0 = -kBarrelHalfZ;
  config.maxZ0 = kBarrelHalfZ;
  config.maxOuterRadius = detector.maxOuterRadius;
  // the toy setup has no ML lookup table and no cluster widths
  config.useMl = false;

  const Experimental::GraphBasedTrackSeeder seeder(
      Experimental::GraphBasedTrackSeeder::DerivedConfig(config), geometry,
      getDefaultLogger("GbtsTest", Logging::Level::WARNING));

  const Experimental::GbtsTrackingFilter filter(
      Experimental::GbtsTrackingFilter::Config{}, geometry);

  const Experimental::GbtsRoiDescriptor roi(0, -4.5, 4.5, 0, -std::numbers::pi,
                                            std::numbers::pi, 0, -kBarrelHalfZ,
                                            kBarrelHalfZ);

  const Experimental::GraphBasedTrackSeeder::Options options(2_T);

  const std::vector<bool> isPixelLayer(numLayers, true);

  SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);

  seeder.createSeeds(spacePoints, roi, isPixelLayer, filter, options, seeds);

  return seeds;
}

/// One line per seed, `quality:sp,sp,...`, sorted by first space point. The
/// container order is not pinned: seeds are sorted by quality with an unstable
/// sort, so ties come out differently between standard libraries.
std::string formatSeeds(const SeedContainer& seeds) {
  std::vector<std::pair<std::uint32_t, std::string>> lines;
  for (const auto& seed : seeds) {
    const auto indices = seed.spacePointIndices();
    std::ostringstream line;
    line << seed.quality() << ":";
    for (const auto index : indices) {
      line << index << ",";
    }
    lines.emplace_back(indices[0], line.str());
  }
  std::ranges::sort(lines);

  std::ostringstream os;
  for (const auto& entry : lines) {
    os << entry.second << "\n";
  }
  return os.str();
}

/// Number of hits each track left in the container.
std::vector<std::size_t> hitsPerTrack(const SpacePointContainer& spacePoints,
                                      std::size_t nTracks) {
  auto trackColumn = spacePoints.column<std::uint32_t>("trackId");
  std::vector<std::size_t> counts(nTracks, 0);
  for (const auto& sp : spacePoints) {
    counts.at(sp.extra(trackColumn)) += 1;
  }
  return counts;
}

/// Invariants that hold independently of the recorded values: a seed collects
/// hits of a single track in radial order.
void checkSeedsAreWellFormed(const SeedContainer& seeds,
                             const SpacePointContainer& spacePoints) {
  auto trackColumn = spacePoints.column<std::uint32_t>("trackId");

  for (const auto& seed : seeds) {
    const auto indices = seed.spacePointIndices();
    BOOST_REQUIRE_GE(indices.size(), 3u);

    BOOST_REQUIRE_LT(indices[0], spacePoints.size());
    const std::uint32_t track = spacePoints.at(indices[0]).extra(trackColumn);

    float previousR = -1.f;
    for (const auto index : indices) {
      BOOST_REQUIRE_LT(index, spacePoints.size());
      const auto sp = spacePoints.at(index);
      BOOST_CHECK_EQUAL(sp.extra(trackColumn), track);
      BOOST_CHECK_GT(sp.r(), previousR);
      previousR = sp.r();
    }
  }
}

/// Exactly one seed per track, holding every hit of that track.
void checkOneSeedPerTrack(const SeedContainer& seeds,
                          const SpacePointContainer& spacePoints,
                          const std::vector<Track>& tracks) {
  BOOST_CHECK_EQUAL(seeds.size(), tracks.size());
  checkSeedsAreWellFormed(seeds, spacePoints);

  auto trackColumn = spacePoints.column<std::uint32_t>("trackId");
  const std::vector<std::size_t> hits =
      hitsPerTrack(spacePoints, tracks.size());

  std::vector<std::size_t> seedsPerTrack(tracks.size(), 0);
  for (const auto& seed : seeds) {
    const auto indices = seed.spacePointIndices();
    const std::uint32_t track = spacePoints.at(indices[0]).extra(trackColumn);
    BOOST_CHECK_EQUAL(indices.size(), hits.at(track));
    seedsPerTrack.at(track) += 1;
  }
  for (const std::size_t count : seedsPerTrack) {
    BOOST_CHECK_EQUAL(count, 1u);
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(GbtsSeeding)

// Guards the fixture itself: a broken input must not look like a regression.
BOOST_AUTO_TEST_CASE(BarrelInputIsWellFormed) {
  const ToyDetector detector = barrelDetector();
  const std::vector<Track> tracks = makeSparseTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  BOOST_CHECK_EQUAL(spacePoints.size(), tracks.size() * detector.layers.size());

  auto layerColumn = spacePoints.column<std::uint32_t>("layerId");
  for (const auto& sp : spacePoints) {
    BOOST_CHECK_LT(sp.extra(layerColumn), detector.layers.size());
    BOOST_CHECK_LE(std::abs(sp.z()), kBarrelHalfZ);
  }
}

// One seed per track, all four hits, in radial order.
BOOST_AUTO_TEST_CASE(SeedsFromSeparatedTracks) {
  const ToyDetector detector = barrelDetector();
  const std::vector<Track> tracks = makeSparseTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const SeedContainer seeds = runSeeding(detector, spacePoints);

  checkOneSeedPerTrack(seeds, spacePoints, tracks);

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
  const ToyDetector detector = barrelDetector();
  const std::vector<Track> tracks = makeDenseTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const SeedContainer seeds = runSeeding(detector, spacePoints);

  // dumped rather than pinned: under heavy branching which candidates survive
  // depends on float rounding and differs between platforms
  BOOST_TEST_MESSAGE("dense seeds:\n" << formatSeeds(seeds));

  // clone removal resolves the branching back to one complete seed per track
  checkOneSeedPerTrack(seeds, spacePoints, tracks);
}

// Guards the forward fixture: every track must reach the four layers a seed
// needs, and the tau range must straddle the barrel-endcap transition.
BOOST_AUTO_TEST_CASE(ForwardInputIsWellFormed) {
  const ToyDetector detector = forwardDetector();
  const std::vector<Track> tracks = makeForwardTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  for (const std::size_t count : hitsPerTrack(spacePoints, tracks.size())) {
    BOOST_CHECK_GE(count, 4u);
  }

  auto layerColumn = spacePoints.column<std::uint32_t>("layerId");
  auto trackColumn = spacePoints.column<std::uint32_t>("trackId");

  std::vector<bool> hasBarrel(tracks.size(), false);
  std::vector<bool> hasEndcap(tracks.size(), false);
  for (const auto& sp : spacePoints) {
    const LayerSpec& layer = detector.layers.at(sp.extra(layerColumn));
    const std::uint32_t track = sp.extra(trackColumn);
    if (layer.type == GbtsLayerType::Barrel) {
      hasBarrel[track] = true;
      BOOST_CHECK_LE(std::abs(sp.z()), kBarrelHalfZ);
    } else {
      hasEndcap[track] = true;
      BOOST_CHECK_GE(sp.r(), kDiscMinR);
      BOOST_CHECK_LE(sp.r(), kDiscMaxR);
    }
  }

  auto isSet = [](bool value) { return value; };
  BOOST_CHECK(std::ranges::all_of(hasEndcap, isSet));
  BOOST_CHECK(std::ranges::any_of(hasBarrel, isSet));
  BOOST_CHECK(!std::ranges::all_of(hasBarrel, isSet));
}

// Covers the endcap eta binning and the barrel-endcap and endcap-endcap bin
// compatibility. The layer type branches in the tracking filter and in the
// triplet cuts also run, but ideal tracks pass every cut by the same margin, so
// a change there does not move these values.
BOOST_AUTO_TEST_CASE(SeedsFromForwardTracks) {
  const ToyDetector detector = forwardDetector();
  const std::vector<Track> tracks = makeForwardTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const SeedContainer seeds = runSeeding(detector, spacePoints);

  checkOneSeedPerTrack(seeds, spacePoints, tracks);

  // the five hit seeds are the tracks that clear the barrel early enough to
  // cross all four discs
  const std::string expected =
      "-10.5:0,1,2,3,\n"
      "-11.2:4,5,6,7,8,\n"
      "-11.2:9,10,11,12,13,\n"
      "-11.2:14,15,16,17,18,\n"
      "-11.2:19,20,21,22,23,\n"
      "-10.5:24,25,26,27,\n"
      "-10.5:28,29,30,31,\n"
      "-10.5:32,33,34,35,\n"
      "-10.5:36,37,38,39,\n"
      "-11.2:40,41,42,43,44,\n"
      "-11.2:45,46,47,48,49,\n"
      "-11.2:50,51,52,53,54,\n"
      "-11.2:55,56,57,58,59,\n"
      "-10.5:60,61,62,63,\n"
      "-10.5:64,65,66,67,\n"
      "-10.5:68,69,70,71,\n";
  BOOST_CHECK_EQUAL(formatSeeds(seeds), expected);
}

// Dense forward tracks add the endcap combinatorics on top.
BOOST_AUTO_TEST_CASE(SeedsFromDenseForwardTracks) {
  const ToyDetector detector = forwardDetector();
  const std::vector<Track> tracks = makeDenseForwardTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const SeedContainer seeds = runSeeding(detector, spacePoints);

  BOOST_TEST_MESSAGE("dense forward seeds:\n" << formatSeeds(seeds));

  checkOneSeedPerTrack(seeds, spacePoints, tracks);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
