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
#include "Acts/Seeding/GbtsGraphBuilder.hpp"
#include "Acts/Seeding/GbtsLayerConnection.hpp"
#include "Acts/Seeding/GbtsProfiler.hpp"
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
#include <ranges>
#include <set>
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
using Experimental::GbtsLayerTechnology;
using Experimental::GbtsLayerType;

/// Half-length in z of every barrel layer, reused as the z0 and RoI z range.
constexpr float kBarrelHalfZ = 150.f;

/// Radial extent of every endcap disc.
constexpr float kDiscMinR = 30.f;
constexpr float kDiscMaxR = 220.f;

/// Eta bin width the layers are split into, as in ATLAS'
/// createLinkingScheme.py.
constexpr float kEtaBinWidth = 0.2f;

/// One layer of the toy detector. The seeder does not read the ids: it orders
/// the pixel barrel layers by radius.
struct LayerSpec {
  Experimental::GbtsExperimentLayerId id{};
  GbtsLayerType type{};
  GbtsLayerTechnology technology{GbtsLayerTechnology::Pixel};
  /// r for a barrel layer, z for an endcap disc.
  float refCoord{};
  /// z range for a barrel layer, r range for an endcap disc.
  float minBound{};
  float maxBound{};
};

constexpr LayerSpec barrelLayer(Experimental::GbtsExperimentLayerId id,
                                float radius) {
  return {.id = id,
          .type = GbtsLayerType::Barrel,
          .refCoord = radius,
          .minBound = -kBarrelHalfZ,
          .maxBound = kBarrelHalfZ};
}

constexpr LayerSpec discLayer(Experimental::GbtsExperimentLayerId id, float z) {
  return {.id = id,
          .type = GbtsLayerType::Endcap,
          .refCoord = z,
          .minBound = kDiscMinR,
          .maxBound = kDiscMaxR};
}

/// Layers plus the connector links between them, outer (src) to inner (dst).
struct ToyDetector {
  std::vector<LayerSpec> layers;
  std::vector<std::pair<Experimental::GbtsExperimentLayerId,
                        Experimental::GbtsExperimentLayerId>>
      links;
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

  for (const Experimental::GbtsExperimentLayerId idBase : {90000u, 70000u}) {
    const float sign = idBase == 90000 ? 1.f : -1.f;
    Experimental::GbtsExperimentLayerId previousId = 0;

    for (std::size_t i = 0; i < discZ.size(); ++i) {
      const Experimental::GbtsExperimentLayerId id =
          idBase + 1000 * static_cast<Experimental::GbtsExperimentLayerId>(i);
      detector.layers.push_back(discLayer(id, sign * discZ[i]));

      if (i == 0) {
        // a track reaching the first disc still has hits in one of these
        // barrel layers, never in 83000
        for (const Experimental::GbtsExperimentLayerId innerId :
             {80000u, 81000u, 82000u}) {
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

std::shared_ptr<Experimental::GbtsGeometry> makeGeometry(
    const ToyDetector& detector) {
  std::vector<Experimental::GbtsLayerDescription> layers;
  layers.reserve(detector.layers.size());
  for (const LayerSpec& spec : detector.layers) {
    Experimental::GbtsLayerDescription layer;
    layer.id = spec.id;
    layer.type = spec.type;
    layer.technology = spec.technology;
    layer.refCoord = spec.refCoord;
    layer.minBound = spec.minBound;
    layer.maxBound = spec.maxBound;
    layers.push_back(layer);
  }

  std::vector<Experimental::GbtsLayerConnection> connections;
  connections.reserve(detector.links.size());
  for (const auto& [src, dst] : detector.links) {
    connections.push_back({.src = src, .dst = dst});
  }
  Experimental::assignAdaptiveTauRatioCuts(connections, layers);

  return std::make_shared<Experimental::GbtsGeometry>(layers, connections,
                                                      kEtaBinWidth);
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

  auto layerColumn =
      container.createColumn<Experimental::GbtsLayerIndex>("gbtsLayerIndex");
  auto clusterWidthColumn = container.createColumn<float>("clusterWidth");
  auto localPositionColumn = container.createColumn<float>("localPositionY");
  auto trackColumn = container.createColumn<std::uint32_t>("trackId");

  container.reserve(
      static_cast<std::uint32_t>(tracks.size() * detector.layers.size()));

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
      sp.extra(layerColumn) = static_cast<Experimental::GbtsLayerIndex>(layer);
      sp.extra(clusterWidthColumn) = 0.f;
      sp.extra(localPositionColumn) = 0.f;
      sp.extra(trackColumn) = static_cast<std::uint32_t>(track);
    }
  }

  return container;
}

/// A strip module of the toy barrel, roughly the ITk strip barrel.
constexpr float kStereoAngle = 26e-3f;
constexpr float kModuleGap = 2.f;
constexpr float kStripHalfLength = 24.f;

/// The same hits, with every strip layer of the detector displaced along its
/// strip by `walk` and carrying the stereo pair it was formed from.
///
/// @param detector the toy detector
/// @param tracks the tracks crossing it
/// @param walk how far along the strip the beam spot puts the point
/// @return the space points
SpacePointContainer makeStripSpacePoints(const ToyDetector& detector,
                                         const std::vector<Track>& tracks,
                                         float walk) {
  SpacePointContainer container(
      SpacePointColumns::CopiedFromIndex | SpacePointColumns::X |
      SpacePointColumns::Y | SpacePointColumns::Z | SpacePointColumns::R |
      SpacePointColumns::Phi | SpacePointColumns::StripCalibrationDetails);

  auto layerColumn =
      container.createColumn<Experimental::GbtsLayerIndex>("gbtsLayerIndex");
  auto clusterWidthColumn = container.createColumn<float>("clusterWidth");
  auto localPositionColumn = container.createColumn<float>("localPositionY");
  auto trackColumn = container.createColumn<std::uint32_t>("trackId");

  container.reserve(
      static_cast<std::uint32_t>(tracks.size() * detector.layers.size()));

  for (std::size_t track = 0; track < tracks.size(); ++track) {
    const float cosPhi = std::cos(tracks[track].phi);
    const float sinPhi = std::sin(tracks[track].phi);
    // the toy tracks are straight, so this is the direction everywhere on one
    const std::array<float, 3> direction{cosPhi, sinPhi, tracks[track].tau};

    for (std::size_t layer = 0; layer < detector.layers.size(); ++layer) {
      const auto crossing = intersect(detector.layers[layer], tracks[track]);
      if (!crossing.has_value()) {
        continue;
      }
      const auto [r, z] = *crossing;
      const std::array<float, 3> point{r * cosPhi, r * sinPhi, z};
      const bool strip =
          detector.layers[layer].technology == GbtsLayerTechnology::Strip;

      auto sp = container.createSpacePoint();
      sp.x() = point[0];
      sp.y() = point[1];
      sp.z() = strip ? z + walk : z;
      sp.r() = r;
      sp.phi() = tracks[track].phi;
      sp.copiedFromIndex() = sp.index();
      // the dense layer index, not the GBTS layer id
      sp.extra(layerColumn) = static_cast<Experimental::GbtsLayerIndex>(layer);
      sp.extra(clusterWidthColumn) = 0.f;
      sp.extra(localPositionColumn) = 0.f;
      sp.extra(trackColumn) = static_cast<std::uint32_t>(track);
      if (!strip) {
        continue;
      }

      // strips along z, the two sensors half the gap either side of the
      // crossing and each rotated half the stereo angle about the normal
      const std::array<float, 3> across{sinPhi, -cosPhi, 0.f};
      const float half = 0.5f * kStereoAngle;
      OuterStripSpacePointCalibrationDetails details{};
      for (std::size_t side = 0; side < 2; ++side) {
        const float sign = side == 0 ? -1.f : 1.f;
        const std::array<float, 3> axis{sign * std::sin(half) * across[0],
                                        sign * std::sin(half) * across[1],
                                        std::cos(half)};
        std::array<float, 3>& halfVector =
            side == 0 ? details.innerHalfVector : details.outerHalfVector;
        for (std::size_t i = 0; i < 3; ++i) {
          halfVector[i] = axis[i] * kStripHalfLength;
        }
        if (side == 1) {
          for (std::size_t i = 0; i < 3; ++i) {
            // the crossing sits at the centre of its strip
            details.outerCenter[i] =
                point[i] + 0.5f * kModuleGap * direction[i];
            details.innerToOuterSeparation[i] = kModuleGap * direction[i];
          }
        }
      }
      sp.outerStripCalibrationDetails() = details;
    }
  }

  return container;
}

/// The seeder and everything it needs, for one toy detector.
struct SeederSetup {
  Experimental::GraphBasedTrackSeeder seeder;
  Experimental::GbtsGraphBuilder graphBuilder;
  Experimental::GbtsTrackingFilter filter;
  Experimental::GbtsRoiDescriptor roi;
  Experimental::GraphBasedTrackSeeder::Options options;
};

/// @param detector the toy detector
/// @param calibrateStrips whether to resolve a doublet's strip ends
/// @param quiet mute the loggers, for a case that drives the graph empty on
///              purpose: the log failure threshold makes the warning that says
///              so an error, and a level above it cannot be asked for either
/// @param profiler handed to the seeder, if not null
SeederSetup makeSeeder(
    const ToyDetector& detector, const bool calibrateStrips = true,
    const bool quiet = false,
    std::shared_ptr<Experimental::GbtsProfiler> profiler = nullptr) {
  auto geometry = makeGeometry(detector);

  const auto makeLogger = [quiet]() -> std::unique_ptr<const Logger> {
    return quiet ? getDummyLogger().clone()
                 : getDefaultLogger("GbtsTest", Logging::Level::WARNING);
  };

  Experimental::GraphBasedTrackSeeder::Config config;
  // the toy setup has no tau lookup table and no cluster widths
  config.useClusterWidthCuts = false;

  Experimental::GbtsGraphBuilder::Config graphConfig;
  graphConfig.minPt = 1_GeV;
  graphConfig.minZ0 = -kBarrelHalfZ;
  graphConfig.maxZ0 = kBarrelHalfZ;
  graphConfig.maxOuterRadius = detector.maxOuterRadius;
  graphConfig.calibrateStrips = calibrateStrips;

  return SeederSetup{
      .seeder = Experimental::GraphBasedTrackSeeder(
          Experimental::GraphBasedTrackSeeder::DerivedConfig(config), geometry,
          makeLogger(), std::move(profiler)),
      .graphBuilder =
          Experimental::GbtsGraphBuilder(graphConfig, geometry, makeLogger()),
      .filter = Experimental::GbtsTrackingFilter(
          Experimental::GbtsTrackingFilter::Config{}, geometry, makeLogger()),
      .roi = Experimental::GbtsRoiDescriptor(-4.5, 4.5, -kBarrelHalfZ,
                                             kBarrelHalfZ),
      .options = Experimental::GraphBasedTrackSeeder::Options{.bFieldInZ = 2_T},
  };
}

SeedContainer runSeeding(const ToyDetector& detector,
                         const SpacePointContainer& spacePoints) {
  const SeederSetup setup = makeSeeder(detector);

  SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);

  setup.seeder.createSeeds(spacePoints, setup.roi, setup.graphBuilder,
                           setup.filter, setup.options, seeds);

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

  auto layerColumn =
      spacePoints.column<Experimental::GbtsLayerIndex>("gbtsLayerIndex");
  for (const auto& sp : spacePoints) {
    BOOST_CHECK_LT(sp.extra(layerColumn), detector.layers.size());
    BOOST_CHECK_LE(std::abs(sp.z()), kBarrelHalfZ);
  }
}

// The binning the geometry hands out is what an external consumer flattens
// into a device array, so it has to tile the bin numbering exactly.
BOOST_AUTO_TEST_CASE(LayerBinningTilesTheBins) {
  const ToyDetector detector = barrelDetector();
  const auto geometry = makeGeometry(detector);

  BOOST_REQUIRE_EQUAL(geometry->numLayers(), detector.layers.size());

  std::uint32_t nextBin = 0;
  for (Experimental::GbtsLayerIndex i = 0; i < geometry->numLayers(); ++i) {
    const Experimental::GbtsLayerBinning& binning = geometry->layerBinning(i);

    BOOST_CHECK_GE(binning.numBins, 1u);
    BOOST_CHECK_GT(binning.etaBinWidth, 0.f);
    // contiguous, in layer order, with no gap or overlap
    BOOST_CHECK_EQUAL(binning.firstBin, nextBin);
    nextBin += binning.numBins;

    BOOST_CHECK_EQUAL(geometry->layerDescription(i).id, detector.layers[i].id);
  }
  BOOST_CHECK_EQUAL(nextBin, geometry->numBins());

  BOOST_CHECK(!geometry->binGroups().empty());
  for (const Experimental::GbtsBinGroup& group : geometry->binGroups()) {
    BOOST_CHECK_LT(group.bin, geometry->numBins());
    for (const Experimental::GbtsBinLink& link : group.links) {
      BOOST_CHECK_LT(link.bin, geometry->numBins());
    }
  }
}

// Filling the node storage from a caller's own EDM has to give the same seeds
// as the overload that fills it from a space point container.
BOOST_AUTO_TEST_CASE(SeedsFromCallerFilledNodeStorage) {
  const ToyDetector detector = barrelDetector();
  const std::vector<Track> tracks = makeSparseTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const SeederSetup setup = makeSeeder(detector);

  Experimental::GbtsNodeStorage storage = setup.seeder.makeNodeStorage();

  // r and phi as stored, so both paths see the same numbers - deriving them
  // from x and y instead costs the last bits and the seed quality with it
  auto layerColumn =
      spacePoints.column<Experimental::GbtsLayerIndex>("gbtsLayerIndex");
  for (const auto& sp : spacePoints) {
    const std::optional<std::uint32_t> bin =
        storage.insert(sp.index(), sp.x(), sp.y(), sp.z(), sp.r(), sp.phi(),
                       sp.extra(layerColumn));
    BOOST_CHECK(bin.has_value());
  }
  storage.finalize();

  BOOST_REQUIRE_EQUAL(storage.numberOfNodes(), spacePoints.size());

  // finalize reorders by (eta bin, phi), so this has to be a permutation
  std::vector<SpacePointIndex> mapped;
  for (std::uint32_t node = 0; node < storage.numberOfNodes(); ++node) {
    mapped.push_back(storage.spacePointIndex(node));
  }
  std::ranges::sort(mapped);
  BOOST_CHECK(std::ranges::equal(
      mapped, std::views::iota(SpacePointIndex{0},
                               static_cast<SpacePointIndex>(mapped.size()))));

  SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);
  setup.seeder.createSeeds(storage, setup.roi, setup.graphBuilder, setup.filter,
                           setup.options, seeds);

  BOOST_CHECK_EQUAL(formatSeeds(seeds),
                    formatSeeds(runSeeding(detector, spacePoints)));
}

// Nodes with identical phi are ordered by position, so the node order does not
// depend on the order in which the caller inserts them.
BOOST_AUTO_TEST_CASE(NodeOrderIgnoresInsertionOrderAtEqualPhi) {
  const ToyDetector detector = barrelDetector();
  const SeederSetup setup = makeSeeder(detector);

  const float r = detector.layers.front().refCoord;
  const float phi = 0.25f;
  const std::array<float, 4> zs = {1.f, -2.f, 3.f, 0.5f};

  const auto nodeOrder = [&](bool reversed) {
    Experimental::GbtsNodeStorage storage = setup.seeder.makeNodeStorage();
    for (std::size_t k = 0; k < zs.size(); ++k) {
      const std::size_t i = reversed ? zs.size() - 1 - k : k;
      storage.insert(static_cast<SpacePointIndex>(i), r * std::cos(phi),
                     r * std::sin(phi), zs[i], r, phi, 0);
    }
    storage.finalize();

    std::vector<SpacePointIndex> order;
    for (std::uint32_t node = 0; node < storage.numberOfNodes(); ++node) {
      order.push_back(storage.spacePointIndex(node));
    }
    return order;
  };

  const std::vector<SpacePointIndex> forward = nodeOrder(false);
  BOOST_CHECK_EQUAL(forward.size(), zs.size());
  BOOST_CHECK(forward == nodeOrder(true));
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

  auto layerColumn =
      spacePoints.column<Experimental::GbtsLayerIndex>("gbtsLayerIndex");
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

// The walk throws off the tau of every doublet through the point, and the
// linking stage rejects on it long before the tracking filter does. 3 mm is
// clear of both ends: the triplet tau ratio gives up past 1.9 mm, the filter's
// chi2 against the nominal positions past 4.4 mm.
BOOST_AUTO_TEST_CASE(StripLayersNeedTheirPairResolved) {
  const std::vector<Track> tracks = makeSparseTracks();
  constexpr float kWalk = 3.f;
  constexpr std::size_t kFirstStripLayer = 2;

  ToyDetector detector = barrelDetector();
  for (std::size_t layer = kFirstStripLayer; layer < detector.layers.size();
       ++layer) {
    detector.layers[layer].technology = GbtsLayerTechnology::Strip;
  }

  const auto seed = [&](const float walk, const bool calibrate) {
    const SpacePointContainer spacePoints =
        makeStripSpacePoints(detector, tracks, walk);
    // a graph with no connections is the point of the unresolved case, so the
    // warning that says so is not one here
    const SeederSetup setup = makeSeeder(detector, calibrate, /*quiet=*/true);
    SeedContainer seeds;
    seeds.assignSpacePointContainer(spacePoints);
    setup.seeder.createSeeds(spacePoints, setup.roi, setup.graphBuilder,
                             setup.filter, setup.options, seeds);
    return seeds.size();
  };

  // a pair that has not walked is seeded either way
  BOOST_CHECK_EQUAL(seed(0.f, false), tracks.size());
  BOOST_CHECK_EQUAL(seed(0.f, true), tracks.size());

  // unresolved, the walk loses every seed; resolved, it loses none
  BOOST_CHECK_EQUAL(seed(kWalk, false), 0u);
  BOOST_CHECK_EQUAL(seed(kWalk, true), tracks.size());
}

BOOST_AUTO_TEST_CASE(AdaptiveTauRatioCuts) {
  using Experimental::GbtsLayerConnection;

  const auto layer = [](Experimental::GbtsExperimentLayerId id,
                        GbtsLayerType type, GbtsLayerTechnology technology,
                        float refCoord) {
    Experimental::GbtsLayerDescription description;
    description.id = id;
    description.type = type;
    description.technology = technology;
    description.refCoord = refCoord;
    return description;
  };
  // listed out of radial order, so the barrel order has to be derived
  const std::vector<Experimental::GbtsLayerDescription> layers = {
      layer(2, GbtsLayerType::Barrel, GbtsLayerTechnology::Pixel, 100.f),
      layer(1, GbtsLayerType::Barrel, GbtsLayerTechnology::Pixel, 50.f),
      layer(3, GbtsLayerType::Barrel, GbtsLayerTechnology::Pixel, 150.f),
      layer(4, GbtsLayerType::Endcap, GbtsLayerTechnology::Pixel, 400.f),
      layer(5, GbtsLayerType::Barrel, GbtsLayerTechnology::Strip, 400.f),
  };

  std::vector<GbtsLayerConnection> connections = {
      {.src = 2, .dst = 1},  // adjacent
      {.src = 3, .dst = 1},  // skips a layer
      {.src = 4, .dst = 3},  // leaves the barrel
      {.src = 4, .dst = 4},  // off the barrel
      {.src = 5, .dst = 3},  // leaves the barrel onto a strip
  };
  const Experimental::GbtsAdaptiveTauRatioCuts cuts;
  Experimental::assignAdaptiveTauRatioCuts(connections, layers, cuts);

  const float base = cuts.tauRatioCut;
  const float corr = cuts.tauRatioCorr;
  const float strip = cuts.tauRatioCorrStrip;
  const std::vector<std::pair<float, float>> expected = {
      {base, base}, {base + corr, base + corr},          {base, base + corr},
      {base, base}, {base + strip, base + corr + strip},
  };
  for (std::size_t i = 0; i < connections.size(); ++i) {
    BOOST_CHECK_CLOSE(connections[i].tauRatioCutInner, expected[i].first, 1e-4);
    BOOST_CHECK_CLOSE(connections[i].tauRatioCutOuter, expected[i].second,
                      1e-4);
  }
}

// Profiling records what the cuts do and must not change what they do. Dense
// tracks so that every stage rejects something.
BOOST_AUTO_TEST_CASE(ProfilerLeavesTheSeedsAlone) {
  const ToyDetector detector = barrelDetector();
  const std::vector<Track> tracks = makeDenseTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  auto profiler = std::make_shared<Experimental::GbtsProfiler>();
  SeederSetup setup = makeSeeder(detector, true, false, profiler);
  setup.options.eventNumber = 42;

  SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);
  setup.seeder.createSeeds(spacePoints, setup.roi, setup.graphBuilder,
                           setup.filter, setup.options, seeds);

  BOOST_CHECK_EQUAL(formatSeeds(seeds),
                    formatSeeds(runSeeding(detector, spacePoints)));

  BOOST_CHECK(!profiler->take(41).has_value());
  std::optional<Experimental::GbtsProfiler::Event> event = profiler->take(42);
  BOOST_REQUIRE(event.has_value());
  BOOST_CHECK(profiler->takeAll().empty());
  BOOST_CHECK_EQUAL(event->eventNumber, 42u);

  BOOST_CHECK(event->rejectedNodes.empty());
  BOOST_CHECK(!event->rejectedDoublets.empty());
  BOOST_CHECK(!event->rejectedTriplets.empty());
  BOOST_CHECK(!event->links.empty());

  // a doublet is either an edge or rejected, never both
  std::set<std::pair<SpacePointIndex, SpacePointIndex>> edges;
  for (const auto& edge : event->edges) {
    BOOST_CHECK(edges.emplace(edge.inner, edge.outer).second);
  }
  for (const auto& doublet : event->rejectedDoublets) {
    BOOST_CHECK(!edges.contains({doublet.inner, doublet.outer}));
  }
  // a link joins two edges
  for (const auto& link : event->links) {
    BOOST_CHECK(edges.contains({link.inner, link.middle}));
    BOOST_CHECK(edges.contains({link.middle, link.outer}));
  }
  for (const auto& triplet : event->rejectedTriplets) {
    BOOST_CHECK(edges.contains({triplet.inner, triplet.middle}));
    BOOST_CHECK(edges.contains({triplet.middle, triplet.outer}));
  }

  // every edge is the head of a path once
  BOOST_REQUIRE_EQUAL(event->heads.size(), event->edges.size());
  std::size_t nCandidateHeads = 0;
  for (const auto& head : event->heads) {
    BOOST_CHECK(edges.contains({head.inner, head.outer}));
    nCandidateHeads +=
        head.outcome == Experimental::GbtsProfiler::HeadOutcome::Candidate;
  }

  // a filter rejection's chain runs from the failed edge out to its head, edge
  // by edge
  for (const auto& rejection : event->filterRejections) {
    const auto& chain = rejection.spacePoints;
    BOOST_REQUIRE_EQUAL(chain.size(), rejection.pathLength + 2u);
    BOOST_CHECK_EQUAL(chain[0], rejection.inner);
    BOOST_CHECK_EQUAL(chain[1], rejection.outer);
    BOOST_CHECK_EQUAL(chain[chain.size() - 2], rejection.headInner);
    BOOST_CHECK_EQUAL(chain.back(), rejection.headOuter);
    for (std::size_t k = 0; k + 1 < chain.size(); ++k) {
      BOOST_CHECK(edges.contains({chain[k], chain[k + 1]}));
    }
  }

  // the candidates are what the clone removal was given, and the output what
  // it kept
  BOOST_CHECK_EQUAL(event->candidates.size(), nCandidateHeads);
  std::size_t nOutput = 0;
  std::size_t nClones = 0;
  for (const auto& candidate : event->candidates) {
    using enum Experimental::GbtsProfiler::CandidateFate;
    nOutput += candidate.fate == Seed ? 1 : candidate.fate == Split ? 2 : 0;
    nClones += candidate.fate == Clone;
  }
  BOOST_CHECK_EQUAL(nOutput, seeds.size());
  BOOST_CHECK_GT(nClones, 0u);

  // a clone shares more than the threshold with better candidates, the one
  // it is a clone of among them; a kept candidate shares no more than that
  for (const auto& candidate : event->candidates) {
    using enum Experimental::GbtsProfiler::CandidateFate;
    if (candidate.fate == TooFewSpacePoints) {
      continue;
    }
    const float limit = candidate.shareThreshold * candidate.spacePoints.size();
    if (candidate.fate == Clone) {
      BOOST_CHECK_GT(candidate.nShared, limit);
      BOOST_REQUIRE_GE(candidate.sharedWith, 0);
      const auto& owner =
          event->candidates.at(static_cast<std::size_t>(candidate.sharedWith));
      BOOST_CHECK(owner.fate != TooFewSpacePoints);
      BOOST_CHECK_LE(owner.quality, candidate.quality);
    } else {
      BOOST_CHECK_LE(candidate.nShared, limit);
    }
    BOOST_CHECK_EQUAL(candidate.nShared == 0, candidate.sharedWith < 0);
  }
}

// The unresolved strip walk of `StripLayersNeedTheirPairResolved` loses every
// seed to the tau ratio of the links, so that is what the profiler must name.
BOOST_AUTO_TEST_CASE(ProfilerRecordsTheCutThatRejects) {
  const std::vector<Track> tracks = makeSparseTracks();
  constexpr float kWalk = 3.f;
  constexpr std::size_t kFirstStripLayer = 2;

  ToyDetector detector = barrelDetector();
  for (std::size_t layer = kFirstStripLayer; layer < detector.layers.size();
       ++layer) {
    detector.layers[layer].technology = GbtsLayerTechnology::Strip;
  }

  const SpacePointContainer spacePoints =
      makeStripSpacePoints(detector, tracks, kWalk);
  auto profiler = std::make_shared<Experimental::GbtsProfiler>();
  const SeederSetup setup =
      makeSeeder(detector, /*calibrateStrips=*/false, /*quiet=*/true, profiler);
  SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);
  setup.seeder.createSeeds(spacePoints, setup.roi, setup.graphBuilder,
                           setup.filter, setup.options, seeds);
  BOOST_CHECK_EQUAL(seeds.size(), 0u);

  std::vector<Experimental::GbtsProfiler::Event> events = profiler->takeAll();
  BOOST_REQUIRE_EQUAL(events.size(), 1u);
  const Experimental::GbtsProfiler::Event& event = events.front();

  auto trackColumn = spacePoints.column<std::uint32_t>("trackId");
  const auto track = [&](const SpacePointIndex index) {
    return spacePoints.at(index).extra(trackColumn);
  };

  // the outer triplet of a track is formed and rejected on its tau ratio,
  // which leaves the second layer isolated, so the innermost doublet is never
  // made: the first layer cuts on the z0 histogram of the second
  BOOST_CHECK(event.links.empty());
  std::size_t nTrueTriplets = 0;
  for (const auto& triplet : event.rejectedTriplets) {
    if (track(triplet.inner) != track(triplet.middle) ||
        track(triplet.middle) != track(triplet.outer)) {
      continue;
    }
    ++nTrueTriplets;
    BOOST_CHECK(triplet.cut ==
                Experimental::GbtsProfiler::TripletCut::TauRatio);
    BOOST_CHECK_GT(triplet.value, triplet.threshold);
  }
  BOOST_CHECK_EQUAL(nTrueTriplets, tracks.size());

  std::size_t nTrueDoublets = 0;
  for (const auto& doublet : event.rejectedDoublets) {
    if (track(doublet.inner) != track(doublet.outer)) {
      continue;
    }
    ++nTrueDoublets;
    BOOST_CHECK(doublet.cut ==
                Experimental::GbtsProfiler::DoubletCut::IsolatedOuterNode);
  }
  BOOST_CHECK_EQUAL(nTrueDoublets, tracks.size());
  BOOST_CHECK_EQUAL(event.edges.size(), 2 * tracks.size());

  // so no edge has a path long enough to seed
  BOOST_CHECK(!event.heads.empty());
  for (const auto& head : event.heads) {
    BOOST_CHECK(head.outcome ==
                Experimental::GbtsProfiler::HeadOutcome::ChainTooShort);
  }
  BOOST_CHECK(event.candidates.empty());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
