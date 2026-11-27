// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Throughput benchmark for the orthogonal triplet seeder on a synthetic event.
///
/// The event comes from `ActsFatras::Synthetic` through
/// SyntheticEventOptions.hpp and the triplet cuts from
/// TripletSeedingConfig.hpp, both shared with the grid triplet benchmark, so
/// the two measure the same seeder fed by two different ways of finding
/// candidates:
///
///  - the grid bins space points up front and offers the doublet finder whole
///    neighbouring bins,
///  - this one puts them in a three-dimensional k-d tree and queries it per
///    middle space point for the box the cuts allow.
///
/// The driving loop follows `ActsExamples::OrthogonalTripletSeedingAlgorithm`,
/// which is how ACTS runs this seeder, including its two passes per middle
/// space point: one assuming the track's z increases with radius and one
/// assuming it decreases.
///
/// What the tree offers the seeder is reported alongside the timing, in the
/// same terms the grid benchmark uses: `bottomPairs` and `topPairs` are the
/// space point pairs handed to the doublet finder, and `ranges` counts how
/// often it is entered, which for this seeder is twice per middle space point.
///
/// The time is not in the doublet finder, which the tree hands a comparable
/// number of pairs to and enters far less often than either grid does; it is in
/// the queries themselves, one box query per middle space point and per z
/// ordering over a tree of 220k points, against a grid lookup that is an array
/// index.
///
/// The longitudinal doublet cut is on by default here for the same reason as in
/// the spherical grid: without z bins nothing else bounds it. `--delta-z-max`
/// with a large value reproduces the ATLAS pixel setting of no cut.

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding/CylindricalSpacePointKDTree.hpp"
#include "Acts/Seeding/DoubletSeedFinder.hpp"
#include "Acts/Seeding/TripletSeedFinder.hpp"
#include "Acts/Seeding/TripletSeeder.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/EventCsvWriter.hpp"
#include "ActsFatras/Synthetic/EventGenerator.hpp"
#include "ActsFatras/Synthetic/SeedingTruth.hpp"
#include "ActsTests/CommonHelpers/BenchmarkTools.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "SyntheticEventOptions.hpp"
#include "TripletSeedingConfig.hpp"

namespace po = boost::program_options;
namespace Exp = Acts::Experimental;

using namespace Acts::UnitLiterals;
using namespace ActsTests;
using namespace ActsFatras::Synthetic;
using namespace ActsTests::TripletSeedingConfig;

namespace {

/// The cuts the k-d tree queries with, which have no counterpart in the grid
/// configuration: the grid states its acceptance through bin edges and
/// neighbour windows, this seeder states it as the box it asks the tree for.
struct OrthogonalConfig {
  /// Radial acceptance in mm
  float rMax = 320_mm;
  /// Longitudinal acceptance in mm
  float zMin = -3000_mm;
  float zMax = 3000_mm;
  /// Radial range middle candidates are taken from in mm, the counterpart of
  /// `rRangeMiddleSP` in the grid configuration
  float rMinMiddle = 40_mm;
  float rMaxMiddle = 260_mm;
  /// Longitudinal reach of middle candidates in mm. The grid gets this for free
  /// by not looping over its outermost z bins.
  float zMinMiddle = -2700_mm;
  float zMaxMiddle = 2700_mm;
  /// Azimuthal half-width of the query box. The grid's counterpart is its phi
  /// bin width times the number of neighbours, which for the ITk configuration
  /// comes out at a comparable value.
  float deltaPhiMax = 0.085f;
};

Exp::CylindricalSpacePointKDTree::Options makeTreeOptions(
    const ItkPixelConfig& cfg, const OrthogonalConfig& ortho) {
  Exp::CylindricalSpacePointKDTree::Options options;
  options.rMax = ortho.rMax;
  options.zMin = ortho.zMin;
  options.zMax = ortho.zMax;
  options.phiMin = cfg.phiMin;
  options.phiMax = cfg.phiMax;
  options.collisionRegionMin = cfg.collisionRegionMin;
  options.collisionRegionMax = cfg.collisionRegionMax;
  options.cotThetaMax = cfg.cotThetaMax;
  options.deltaPhiMax = ortho.deltaPhiMax;
  return options;
}

/// Everything that is reused between benchmark iterations, so that the timed
/// region allocates as little as the experiment's own does.
struct SeedingCache {
  Acts::BroadTripletSeedFilter::Config filterConfig;
  std::shared_ptr<const Acts::DoubletSeedFinder> bottomDoubletFinder;
  std::shared_ptr<const Acts::DoubletSeedFinder> topDoubletFinder;
  std::shared_ptr<const Acts::TripletSeedFinder> tripletFinder;
  Acts::TripletSeeder seeder;
  Acts::TripletSeeder::Cache seederCache;
  Exp::CylindricalSpacePointKDTree::Candidates candidates;
  std::unique_ptr<const Acts::Logger> logger;
};

/// Repack an event into the layout the seeder consumes and index it in a k-d
/// tree.
/// @param spacePoints the space points of the event
/// @param coreSpacePoints receives the repacked space points
/// @param rRange receives the radial extent of the event
/// @param logger the logger to build the tree with
/// @return the tree, indexing `coreSpacePoints`
Exp::CylindricalSpacePointKDTree buildTree(
    const Acts::SpacePointContainer& spacePoints,
    Acts::SpacePointContainer& coreSpacePoints, Acts::Extent& rRange,
    const Acts::Logger& logger) {
  Exp::CylindricalSpacePointKDTreeBuilder builder(logger.clone());
  builder.reserve(spacePoints.size());
  coreSpacePoints.reserve(spacePoints.size());

  for (std::size_t i = 0; i < spacePoints.size(); ++i) {
    const auto sp = spacePoints[i];
    const Acts::SpacePointIndex newIndex = coreSpacePoints.size();
    auto newSp = coreSpacePoints.createSpacePoint();
    newSp.copiedFromIndex() = sp.index();
    newSp.xy() = std::array<float, 2>{sp.x(), sp.y()};
    newSp.zr() = std::array<float, 2>{sp.z(), sp.r()};
    newSp.phi() = sp.phi();
    newSp.varianceZ() = sp.varianceZ();
    newSp.varianceR() = sp.varianceR();

    builder.insert(newIndex, newSp.phi(), newSp.zr()[1], newSp.zr()[0]);
    rRange.extend({newSp.xy()[0], newSp.xy()[1], newSp.zr()[0]});
  }

  return builder.build();
}

/// Whether a space point may be the middle one of a seed, and if so how many
/// top candidates the confirmation asks for.
/// @param cfg the seeding configuration
/// @param ortho the k-d tree configuration
/// @param spM the candidate middle space point
/// @param variableMiddleRange the radial range derived from the event extent
/// @param nTopSeedConf receives the required number of top candidates
/// @return whether the space point is an acceptable middle one
bool acceptMiddle(const ItkPixelConfig& cfg, const OrthogonalConfig& ortho,
                  const Acts::ConstSpacePointProxy& spM,
                  const Acts::Range1D<float>& variableMiddleRange,
                  std::size_t& nTopSeedConf) {
  const float rM = spM.zr()[1];
  if (cfg.useVariableMiddleSPRange) {
    if (rM < variableMiddleRange.min() || rM > variableMiddleRange.max()) {
      return false;
    }
  } else if (rM < ortho.rMinMiddle || rM > ortho.rMaxMiddle) {
    return false;
  }

  const float zM = spM.zr()[0];
  if (zM < ortho.zMinMiddle || zM > ortho.zMaxMiddle) {
    return false;
  }
  if (spM.phi() < cfg.phiMin || spM.phi() > cfg.phiMax) {
    return false;
  }

  nTopSeedConf = 0;
  if (cfg.seedConfirmation) {
    // the confirmation asks for more top candidates in the central region than
    // in the forward one, and for fewer at large radius than at small
    const Acts::BroadTripletSeedFilter::Config filterCfg =
        makeFilterConfig(cfg);
    const Acts::SeedConfirmationRangeConfig& range =
        (zM > filterCfg.centralSeedConfirmationRange.zMaxSeedConf ||
         zM < filterCfg.centralSeedConfirmationRange.zMinSeedConf)
            ? filterCfg.forwardSeedConfirmationRange
            : filterCfg.centralSeedConfirmationRange;
    nTopSeedConf =
        rM > range.rMaxSeedConf ? range.nTopForLargeR : range.nTopForSmallR;
  }
  return true;
}

/// Run one full seeding pass over an event, the way the ACTS algorithm does.
/// @param cfg the seeding configuration
/// @param ortho the k-d tree configuration
/// @param cache the reusable finders and their state
/// @param spacePoints the space points of the event
/// @param seeds receives the seeds
void createSeeds(const ItkPixelConfig& cfg, const OrthogonalConfig& ortho,
                 SeedingCache& cache,
                 const Acts::SpacePointContainer& spacePoints,
                 Acts::SeedContainer& seeds) {
  Acts::SpacePointContainer coreSpacePoints(
      Acts::SpacePointColumns::CopiedFromIndex |
      Acts::SpacePointColumns::PackedXY | Acts::SpacePointColumns::PackedZR |
      Acts::SpacePointColumns::Phi | Acts::SpacePointColumns::VarianceZ |
      Acts::SpacePointColumns::VarianceR);
  Acts::Extent rRangeExtent;
  const Exp::CylindricalSpacePointKDTree tree =
      buildTree(spacePoints, coreSpacePoints, rRangeExtent, *cache.logger);

  const Exp::CylindricalSpacePointKDTree::Options lhOptions = [&] {
    Exp::CylindricalSpacePointKDTree::Options options =
        makeTreeOptions(cfg, ortho);
    options.deltaRMin = cfg.deltaRMinBottomSP;
    options.deltaRMax = cfg.deltaRMaxBottomSP;
    return options;
  }();
  const Exp::CylindricalSpacePointKDTree::Options hlOptions = [&] {
    Exp::CylindricalSpacePointKDTree::Options options =
        makeTreeOptions(cfg, ortho);
    options.deltaRMin = cfg.deltaRMinTopSP;
    options.deltaRMax = cfg.deltaRMaxTopSP;
    return options;
  }();

  const Acts::Range1D<float> variableMiddleRange(
      std::floor(rRangeExtent.min(Acts::AxisDirection::AxisR) / 2) * 2 +
          cfg.deltaRMiddleMinSPRange,
      std::floor(rRangeExtent.max(Acts::AxisDirection::AxisR) / 2) * 2 -
          cfg.deltaRMiddleMaxSPRange);

  Acts::BroadTripletSeedFilter::State filterState;
  Acts::BroadTripletSeedFilter::Cache filterCache;
  Acts::BroadTripletSeedFilter filter(cache.filterConfig, filterState,
                                      filterCache, *cache.logger);

  Acts::SeedContainer candidateSeeds;
  candidateSeeds.assignSpacePointContainer(coreSpacePoints);

  for (const auto& middle : tree) {
    const auto spM = coreSpacePoints.at(middle.second).asConst();

    std::size_t nTopSeedConf = 0;
    if (!acceptMiddle(cfg, ortho, spM, variableMiddleRange, nTopSeedConf)) {
      continue;
    }

    cache.candidates.clear();
    tree.validTuples(lhOptions, hlOptions, spM, nTopSeedConf, cache.candidates);

    // two passes: one for a track whose z grows with radius and one for the
    // opposite, which is how this seeder replaces the grid's z binning
    auto bottomSps =
        coreSpacePoints.subset(cache.candidates.bottom_lh_v).asConst();
    auto topSps = coreSpacePoints.subset(cache.candidates.top_lh_v).asConst();
    cache.seeder.createSeedsFromGroup(
        cache.seederCache, *cache.bottomDoubletFinder, *cache.topDoubletFinder,
        *cache.tripletFinder, filter, coreSpacePoints, bottomSps, spM, topSps,
        candidateSeeds);

    bottomSps = coreSpacePoints.subset(cache.candidates.bottom_hl_v).asConst();
    topSps = coreSpacePoints.subset(cache.candidates.top_hl_v).asConst();
    cache.seeder.createSeedsFromGroup(
        cache.seederCache, *cache.bottomDoubletFinder, *cache.topDoubletFinder,
        *cache.tripletFinder, filter, coreSpacePoints, bottomSps, spM, topSps,
        candidateSeeds);
  }

  // map the seeds back onto the original container, and apply the same quality
  // selection the grid benchmark does so that the two are comparable
  seeds.assignSpacePointContainer(spacePoints);
  seeds.reserve(candidateSeeds.size());
  for (Acts::MutableSeedProxy candidate : candidateSeeds) {
    const auto& indices = candidate.spacePointIndices();
    if (cfg.seedQualitySelection) {
      const float quality = candidate.quality();
      const bool best =
          std::ranges::any_of(indices, [&](Acts::SpacePointIndex index) {
            return filterState.bestSeedQualityMap.at(index) <= quality;
          });
      if (!best) {
        continue;
      }
    }
    const std::array<Acts::SpacePointIndex, 3> original{
        coreSpacePoints.at(indices[0]).copiedFromIndex(),
        coreSpacePoints.at(indices[1]).copiedFromIndex(),
        coreSpacePoints.at(indices[2]).copiedFromIndex()};
    auto seed = seeds.createSeed();
    seed.assignSpacePointIndices(original);
    seed.quality() = candidate.quality();
    seed.vertexZ() = candidate.vertexZ();
  }
}

/// Count the space point pairs the tree offers the doublet finders, in the same
/// terms the grid benchmark reports.
/// @param cfg the seeding configuration
/// @param ortho the k-d tree configuration
/// @param spacePoints the space points of the event
/// @param logger the logger to build the tree with
void printTreeStatistics(const ItkPixelConfig& cfg,
                         const OrthogonalConfig& ortho,
                         const Acts::SpacePointContainer& spacePoints,
                         const Acts::Logger& logger) {
  Acts::SpacePointContainer coreSpacePoints(
      Acts::SpacePointColumns::CopiedFromIndex |
      Acts::SpacePointColumns::PackedXY | Acts::SpacePointColumns::PackedZR |
      Acts::SpacePointColumns::Phi | Acts::SpacePointColumns::VarianceZ |
      Acts::SpacePointColumns::VarianceR);
  Acts::Extent rRangeExtent;
  const Exp::CylindricalSpacePointKDTree tree =
      buildTree(spacePoints, coreSpacePoints, rRangeExtent, logger);

  Exp::CylindricalSpacePointKDTree::Options lhOptions =
      makeTreeOptions(cfg, ortho);
  lhOptions.deltaRMin = cfg.deltaRMinBottomSP;
  lhOptions.deltaRMax = cfg.deltaRMaxBottomSP;
  Exp::CylindricalSpacePointKDTree::Options hlOptions =
      makeTreeOptions(cfg, ortho);
  hlOptions.deltaRMin = cfg.deltaRMinTopSP;
  hlOptions.deltaRMax = cfg.deltaRMaxTopSP;

  const Acts::Range1D<float> variableMiddleRange(
      std::floor(rRangeExtent.min(Acts::AxisDirection::AxisR) / 2) * 2 +
          cfg.deltaRMiddleMinSPRange,
      std::floor(rRangeExtent.max(Acts::AxisDirection::AxisR) / 2) * 2 -
          cfg.deltaRMiddleMaxSPRange);

  Exp::CylindricalSpacePointKDTree::Candidates candidates;
  std::size_t middles = 0;
  std::size_t bottomPairs = 0;
  std::size_t topPairs = 0;
  std::size_t ranges = 0;
  for (const auto& middle : tree) {
    const auto spM = coreSpacePoints.at(middle.second).asConst();
    std::size_t nTopSeedConf = 0;
    if (!acceptMiddle(cfg, ortho, spM, variableMiddleRange, nTopSeedConf)) {
      continue;
    }
    ++middles;
    candidates.clear();
    tree.validTuples(lhOptions, hlOptions, spM, nTopSeedConf, candidates);
    bottomPairs +=
        candidates.bottom_lh_v.size() + candidates.bottom_hl_v.size();
    topPairs += candidates.top_lh_v.size() + candidates.top_hl_v.size();
    // the doublet finder is entered twice per middle space point, once per z
    // ordering
    ranges += 2;
  }

  std::cout << "treeSize=" << tree.size() << " middles=" << middles
            << " bottomPairs=" << bottomPairs << " topPairs=" << topPairs
            << " ranges=" << ranges << std::endl;
}

}  // namespace

int main(int argc, char* argv[]) {
  SyntheticEventOptions::Values shared;
  std::size_t numRuns = 10;
  std::size_t numEvents = 1;
  float minPt = ItkPixelConfig{}.minPt / 1_MeV;
  // Efficiency is counted over a harder threshold than the seeder is cut at, so
  // that the turn-on stays out of it.
  float truthPt = 1000.f;
  // The grid gets its longitudinal doublet bound from its z bins, so the ATLAS
  // pixel configuration switches the explicit cut off. A k-d tree has no z bins
  // either, so as with the spherical grid the cut has to be made explicitly.
  float deltaZMax = 900.f;
  std::string dumpPrefix;
  bool verbose = false;

  OrthogonalConfig ortho;

  try {
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message");
    SyntheticEventOptions::add(desc, shared);
    desc.add_options()("runs",
                       po::value<std::size_t>(&numRuns)->default_value(numRuns),
                       "number of benchmark runs")(
        "events", po::value<std::size_t>(&numEvents)->default_value(numEvents),
        "number of distinct events to generate and cycle through")(
        "min-pt", po::value<float>(&minPt)->default_value(minPt),
        "seed momentum threshold in MeV")(
        "truth-pt", po::value<float>(&truthPt)->default_value(truthPt),
        "momentum in MeV a primary has to reach to be counted in the "
        "efficiency")(
        "delta-z-max", po::value<float>(&deltaZMax)->default_value(deltaZMax),
        "longitudinal reach of a doublet in mm; the k-d tree has no z bins to "
        "bound it implicitly")(
        "delta-phi-max",
        po::value<float>(&ortho.deltaPhiMax)->default_value(ortho.deltaPhiMax),
        "azimuthal half-width of the k-d tree query box")(
        "middle-r",
        po::value<float>(&ortho.rMaxMiddle)->default_value(ortho.rMaxMiddle),
        "outer radius middle candidates are taken from in mm")(
        "dump", po::value<std::string>(&dumpPrefix),
        "write the first event to <arg>_spacepoints.csv and "
        "<arg>_particles.csv")("verbose", po::bool_switch(&verbose),
                               "log the seeder's own statistics");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help") > 0) {
      std::cout << desc << std::endl;
      return 0;
    }
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }

  if (numEvents == 0) {
    std::cerr << "error: --events needs at least one event" << std::endl;
    return 1;
  }

  DetectorLayout layout;
  EventConfig eventConfig;
  try {
    layout = SyntheticEventOptions::makeLayout(shared);
    eventConfig = SyntheticEventOptions::makeConfig(shared);
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }

  std::vector<Event> events;
  events.reserve(numEvents);
  {
    EventGenerator generator(layout, eventConfig);
    for (std::size_t i = 0; i < numEvents; ++i) {
      events.push_back(generator.generate());
    }
  }
  if (!dumpPrefix.empty()) {
    writeEventCsv(events.front(), layout, dumpPrefix);
  }

  ItkPixelConfig cfg;
  cfg.bFieldInZ = eventConfig.bFieldZ;
  cfg.minPt = minPt * 1_MeV;
  cfg.deltaZMax = deltaZMax;

  SeedingCache cache;
  cache.logger = Acts::getDefaultLogger(
      "Orthogonal",
      verbose ? Acts::Logging::Level::DEBUG : Acts::Logging::Level::FATAL);
  cache.filterConfig = makeFilterConfig(cfg);

  Acts::DoubletSeedFinder::Config bottomConfig = makeBottomDoubletConfig(cfg);
  // the tree hands over candidates in tree order, not sorted by radius
  bottomConfig.spacePointsSortedByRadius = false;
  Acts::DoubletSeedFinder::Config topConfig = bottomConfig;
  topConfig.candidateDirection = Acts::Direction::Forward();
  topConfig.deltaRMin = cfg.deltaRMinTopSP;
  topConfig.deltaRMax = cfg.deltaRMaxTopSP;
  cache.bottomDoubletFinder = Acts::DoubletSeedFinder::create(
      Acts::DoubletSeedFinder::DerivedConfig(bottomConfig, cfg.bFieldInZ));
  cache.topDoubletFinder = Acts::DoubletSeedFinder::create(
      Acts::DoubletSeedFinder::DerivedConfig(topConfig, cfg.bFieldInZ));
  cache.tripletFinder =
      Acts::TripletSeedFinder::create(Acts::TripletSeedFinder::DerivedConfig(
          makeTripletConfig(cfg), cfg.bFieldInZ));
  cache.seeder = Acts::TripletSeeder(cache.logger->clone());

  printTreeStatistics(cfg, ortho, events.front().spacePoints, *cache.logger);

  // Truth and efficiency over every event, outside the timed region so that the
  // truth lookup does not show up in the measurement.
  EventSummary summary;
  SeedingSummary seedSummary;
  for (const Event& event : events) {
    const EventSummary one = summarize(event, truthPt * 1_MeV / 1_GeV);
    summary.spacePoints += one.spacePoints;
    summary.primaries += one.primaries;
    summary.secondaries += one.secondaries;
    summary.seedablePrimaries += one.seedablePrimaries;
    summary.primaryHits += one.primaryHits;
    summary.secondaryHits += one.secondaryHits;

    Acts::SeedContainer reference;
    createSeeds(cfg, ortho, cache, event.spacePoints, reference);
    const SeedingSummary oneSeeding =
        evaluateSeeds(event, reference, truthPt * 1_MeV / 1_GeV);
    seedSummary.seeds += oneSeeding.seeds;
    seedSummary.trueSeeds += oneSeeding.trueSeeds;
    seedSummary.matchedPrimaries += oneSeeding.matchedPrimaries;
  }

  const auto perEvent = [&](std::size_t value) { return value / numEvents; };
  std::cout << "seeder=orthogonal events=" << numEvents
            << "\nmean/event: spacePoints=" << perEvent(summary.spacePoints)
            << " primaryHits=" << perEvent(summary.primaryHits)
            << " secondaryHits=" << perEvent(summary.secondaryHits)
            << " primaries=" << perEvent(summary.primaries)
            << " secondaries=" << perEvent(summary.secondaries)
            << " seedable=" << perEvent(summary.seedablePrimaries) << std::endl;

  std::size_t next = 0;
  const auto result = microBenchmark(
      [&] {
        Acts::SeedContainer seeds;
        createSeeds(cfg, ortho, cache, events[next].spacePoints, seeds);
        next = (next + 1) % numEvents;
        assumeRead(seeds);
      },
      1, numRuns);

  std::cout << "mean/event: seeds=" << perEvent(seedSummary.seeds)
            << " trueSeeds=" << perEvent(seedSummary.trueSeeds)
            << " efficiency="
            << static_cast<float>(seedSummary.matchedPrimaries) /
                   static_cast<float>(
                       std::max<std::size_t>(1, summary.seedablePrimaries))
            << "\n"
            << result << std::endl;

  return 0;
}
