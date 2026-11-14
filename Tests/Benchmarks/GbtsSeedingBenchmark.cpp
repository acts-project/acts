// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Throughput benchmark for the GBTS seeder on an ITk-like synthetic event.
///
/// The detector and the event come from `ActsFatras::Synthetic` through
/// SyntheticEventOptions.hpp, shared with the other seeding benchmarks. Only
/// the GBTS-specific geometry, namely the logical layer numbering and the layer
/// connection table, is built here.
///
/// Two timings are reported. The one-shot `createSeeds` is what an experiment
/// calls, and it includes filling the node storage. The second times that fill
/// alone, so that the graph stage is the difference between the two. A filled
/// storage cannot be reused instead: the graph builder writes its edge state
/// into the storage and only partly resets it.
///
/// The cuts are the ACTS defaults, which are what ATLAS runs on the ITk pixel
/// detector. The layer connection table below is not: it is written for these
/// layouts rather than trained on them, so the efficiency here is not a
/// statement about GBTS as ATLAS runs it. What it loses sits at the
/// barrel-endcap transition, where a hand-written table is weakest.
///
/// @note A connection table belongs to the layout it was written for. A layout
///       describing the rings of an endcap has many more discs per side than
///       one outlining it, and a table whose reach along z does not cover them
///       costs half the efficiency. See `kDiscReach`.
///
/// @note GBTS returns seeds of four to eleven space points, so it is scored on
///       sharing `evaluateSeeds`' minimum with one primary rather than on every
///       space point matching. Scoring it the strict way costs four points of
///       efficiency that say nothing about the seeder.

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsGraphBuilder.hpp"
#include "Acts/Seeding/GbtsLayerConnection.hpp"
#include "Acts/Seeding/GbtsNodeStorage.hpp"
#include "Acts/Seeding/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/EventGenerator.hpp"
#include "ActsFatras/Synthetic/SeedingTruth.hpp"
#include "ActsTests/CommonHelpers/BenchmarkTools.hpp"

#include <cstdint>
#include <iostream>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "SyntheticEventOptions.hpp"

namespace po = boost::program_options;
namespace Exp = Acts::Experimental;

using namespace Acts::UnitLiterals;
using namespace ActsTests;
using namespace ActsFatras::Synthetic;

namespace {

/// Rings a disc may carry before the endcap encoding below runs out of room.
constexpr int kRingsPerDisc = 8;

/// How many discs back along z the connection table lets a disc reach. One per
/// ring set of the ITk pixel endcap, see `makeConnectionTable`.
constexpr int kDiscReach = 9;

/// GBTS numbers its logical layers as `gbtsId * 1000 + subIndex` and derives
/// the barrel flag from the leading digit as `layerId / 10000 == 8`, so the
/// barrel takes the eighties, one group per cylinder, with the eta module as
/// the sub-index.
///
/// The endcaps cannot do the same: a layout describing the ring structure of a
/// real endcap has far more than ten discs per side, so a group per disc would
/// run out of the nineties and collide with the barrel. Each endcap is instead
/// one group, 90 and 70, with the disc and the ring packed into the sub-index
/// together.
Exp::GbtsExperimentLayerId gbtsLayerId(const DetectorLayer& layer) {
  const Exp::GbtsExperimentLayerId layerIdx = layer.layer;
  const Exp::GbtsExperimentLayerId moduleIdx = layer.moduleIndex;
  if (layer.shape == SurfaceShape::Cylinder) {
    return (80 + layerIdx) * 1000 + moduleIdx;
  }
  const Exp::GbtsExperimentLayerId group =
      layer.side == SurfaceSide::Positive ? 90 : 70;
  return group * 1000 + layerIdx * kRingsPerDisc + moduleIdx;
}

/// The `gbtsLayerId` encoding above only holds within limits: the barrel group
/// numbers have to stay in the eighties for `layerId / 10000 == 8` to mean
/// barrel, and both sub-indices have to fit below the thousands digit.
/// `GbtsLayerIndex` is also 16 bits, which bounds the number of layers.
/// @param layout the layout to check
void checkEncodingFits(const DetectorLayout& layout) {
  int numBarrel = 0;
  int maxBarrelModule = 0;
  int maxDiscSubIndex = 0;
  for (const DetectorLayer& layer : layout.layers) {
    const auto layerIdx = static_cast<int>(layer.layer);
    const auto moduleIdx = static_cast<int>(layer.moduleIndex);
    if (layer.shape == SurfaceShape::Cylinder) {
      numBarrel = std::max(numBarrel, layerIdx + 1);
      maxBarrelModule = std::max(maxBarrelModule, moduleIdx);
    } else {
      if (moduleIdx >= kRingsPerDisc) {
        throw std::invalid_argument(
            "the GBTS layer id encoding takes at most eight rings per disc");
      }
      maxDiscSubIndex =
          std::max(maxDiscSubIndex, layerIdx * kRingsPerDisc + moduleIdx);
    }
  }
  if (numBarrel > 10) {
    throw std::invalid_argument(
        "the GBTS layer id encoding takes at most ten barrel layers");
  }
  if (maxBarrelModule >= 1000) {
    throw std::invalid_argument(
        "the GBTS layer id encoding takes at most a thousand modules per "
        "cylinder");
  }
  if (maxDiscSubIndex >= 1000) {
    throw std::invalid_argument(
        "the GBTS layer id encoding takes at most a hundred and twenty-five "
        "discs per side");
  }
  if (layout.layers.size() > 65535) {
    throw std::invalid_argument("GbtsLayerIndex is 16 bits wide");
  }
}

std::vector<Exp::GbtsLayerDescription> makeLayerDescriptions(
    const DetectorLayout& layout) {
  std::vector<Exp::GbtsLayerDescription> descriptions;
  descriptions.reserve(layout.layers.size());
  for (const DetectorLayer& layer : layout.layers) {
    Exp::GbtsLayerDescription description;
    description.id = gbtsLayerId(layer);
    description.type = layer.shape == SurfaceShape::Cylinder
                           ? Exp::GbtsLayerType::Barrel
                           : Exp::GbtsLayerType::Endcap;
    description.refCoord = layer.refCoord;
    description.minBound = layer.minBound;
    description.maxBound = layer.maxBound;
    descriptions.push_back(description);
  }
  return descriptions;
}

/// The connection table lists which layer may feed which. `src` is the outer
/// layer of a doublet, `dst` the inner one. The per-eta-bin compatibility is
/// worked out by GbtsGeometry from the beam spot constraint, so listing a
/// layer pair here only says that the pair is allowed at all.
std::vector<Exp::GbtsLayerConnection> makeConnections(
    const DetectorLayout& layout) {
  // group layers by (side, layer) so that all modules of a surface are linked
  auto modulesOf = [&](SurfaceSide side, int layerIdx) {
    std::vector<Exp::GbtsExperimentLayerId> ids;
    for (const DetectorLayer& layer : layout.layers) {
      if (layer.side == side &&
          layer.layer == static_cast<std::uint32_t>(layerIdx)) {
        ids.push_back(gbtsLayerId(layer));
      }
    }
    return ids;
  };

  std::vector<Exp::GbtsLayerConnection> connections;
  auto connect = [&](SurfaceSide outerSide, int outerLayer,
                     SurfaceSide innerSide, int innerLayer) {
    for (const Exp::GbtsExperimentLayerId src :
         modulesOf(outerSide, outerLayer)) {
      for (const Exp::GbtsExperimentLayerId dst :
           modulesOf(innerSide, innerLayer)) {
        connections.push_back({src, dst});
      }
    }
  };

  int numBarrel = 0;
  int numDiscs = 0;
  for (const DetectorLayer& layer : layout.layers) {
    if (layer.shape == SurfaceShape::Cylinder) {
      numBarrel = std::max(numBarrel, static_cast<int>(layer.layer) + 1);
    } else if (layer.side == SurfaceSide::Positive) {
      numDiscs = std::max(numDiscs, static_cast<int>(layer.layer) + 1);
    }
  }

  // barrel to barrel, adjacent and one skipped layer
  for (int i = 0; i < numBarrel; ++i) {
    for (int step : {1, 2}) {
      if (i + step < numBarrel) {
        connect(SurfaceSide::Barrel, i + step, SurfaceSide::Barrel, i);
      }
    }
  }
  for (const SurfaceSide side :
       {SurfaceSide::Positive, SurfaceSide::Negative}) {
    // Disc to disc. Discs are numbered outwards along z but consecutive ones do
    // not sit at consecutive radii, a layout describing the rings of an endcap
    // interleaving the ring sets, so `kDiscReach` has to cover the number of
    // ring sets rather than just the neighbours.
    for (int j = 0; j < numDiscs; ++j) {
      for (int step = 1; step <= kDiscReach; ++step) {
        if (j + step < numDiscs) {
          connect(side, j + step, side, j);
        }
      }
    }
    // forward transition: the innermost barrel layers are fed by the first
    // discs of each side, and by as many of them as a disc reaches back over
    for (int j = 0; j < kDiscReach; ++j) {
      for (int i = 0; i < numBarrel; ++i) {
        connect(side, j, SurfaceSide::Barrel, i);
      }
    }
  }

  return connections;
}

/// GBTS reads the layer of a space point from a 16 bit column of its own name,
/// while the generator writes a 32 bit one. Copy it over once per event,
/// outside the timed region; `checkEncodingFits` makes sure that nothing is cut
/// off.
/// @param spacePoints the space points of the event
void addGbtsLayerIndex(Acts::SpacePointContainer& spacePoints) {
  auto layerIds = spacePoints.column<std::uint32_t>("layerId");
  auto gbtsLayers =
      spacePoints.createColumn<Exp::GbtsLayerIndex>("gbtsLayerIndex");
  for (auto sp : spacePoints) {
    sp.extra(gbtsLayers) = static_cast<Exp::GbtsLayerIndex>(sp.extra(layerIds));
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  SyntheticEventOptions::Values shared;
  std::size_t numRuns = 10;
  float minPt = 900.f;
  // Efficiency is counted over a harder threshold than the seeder is cut at, so
  // that the turn-on stays out of it.
  float truthPt = 1000.f;
  bool verbose = false;
  // The graph is capped at this many edges per space point. The synthetic event
  // is noise free, so more geometric doublet candidates survive the cuts than
  // in a real event; the default still leaves room for all of them on the ITk
  // at pile-up 200. Benchmarking against a saturated ceiling would hide any
  // work done past it, so lower it only to measure the truncated regime.
  float maxEdgesPerSP = Exp::GbtsGraphBuilder::Config{}.maxEdgesPerSP;

  try {
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message");
    SyntheticEventOptions::add(desc, shared);
    desc.add_options()("runs",
                       po::value<std::size_t>(&numRuns)->default_value(numRuns),
                       "number of benchmark runs")(
        "min-pt", po::value<float>(&minPt)->default_value(minPt),
        "seed momentum threshold in MeV")(
        "truth-pt", po::value<float>(&truthPt)->default_value(truthPt),
        "momentum in MeV a primary has to reach to be counted in the "
        "efficiency")(
        "max-edges-per-sp",
        po::value<float>(&maxEdgesPerSP)->default_value(maxEdgesPerSP),
        "ceiling on the number of graph edges per space point")(
        "verbose", po::bool_switch(&verbose),
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

  DetectorLayout layout;
  EventConfig eventConfig;
  try {
    layout = SyntheticEventOptions::makeLayout(shared);
    eventConfig = SyntheticEventOptions::makeConfig(shared);
    checkEncodingFits(layout);
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }
  Event event = generateEvent(layout, eventConfig);
  addGbtsLayerIndex(event.spacePoints);
  const Acts::SpacePointContainer& spacePoints = event.spacePoints;

  constexpr float etaBinWidth = 0.2f;
  const std::vector<Exp::GbtsLayerDescription> layers =
      makeLayerDescriptions(layout);
  const std::vector<Exp::GbtsLayerConnection> connections =
      makeConnections(layout);
  const auto geometry =
      std::make_shared<Exp::GbtsGeometry>(layers, connections, etaBinWidth);

  const auto makeLogger = [&](const std::string& name) {
    // Warnings are let through even when quiet: the graph builder stops once
    // `maxEdgesPerSP` is reached and says so at that level, and a truncated
    // graph costs efficiency without failing.
    return Acts::getDefaultLogger(name, verbose
                                            ? Acts::Logging::Level::DEBUG
                                            : Acts::Logging::Level::WARNING);
  };

  // The cuts are the ACTS defaults, which ATLAS runs the ITk pixel detector
  // with: `ActsTrk::GbtsSeedingTool` overrides only the connector file, the ML
  // lookup table and minPt. Two of its defaults differ and are set below.
  Exp::GraphBasedTrackSeeder::Config cfg;
  // ATLAS runs with the cluster width and its trained lookup table, which the
  // synthetic event has no cluster shapes to offer. It only adjusts the barrel
  // cot(theta) cuts, which is not where this event loses anything.
  cfg.useClusterWidthCuts = false;

  Exp::GbtsGraphBuilder::Config graphCfg;
  graphCfg.minPt = minPt * 1_MeV;
  graphCfg.minZ0 = -200.f;
  graphCfg.maxZ0 = 200.f;
  graphCfg.maxOuterRadius = 550.f;
  graphCfg.matchBeforeCreate = true;
  graphCfg.maxEdgesPerSP = maxEdgesPerSP;

  const Exp::GraphBasedTrackSeeder seeder(
      Exp::GraphBasedTrackSeeder::DerivedConfig(cfg), geometry,
      makeLogger("Gbts"));
  const Exp::GbtsGraphBuilder graphBuilder(graphCfg, geometry,
                                           makeLogger("GbtsGraphBuilder"));
  const Exp::GbtsTrackingFilter filter(Exp::GbtsTrackingFilter::Config{},
                                       geometry, makeLogger("GbtsFilter"));
  const Exp::GbtsRoiDescriptor roi(-5., 5., graphCfg.minZ0, graphCfg.maxZ0);
  const Exp::GraphBasedTrackSeeder::Options options{.bFieldInZ =
                                                        eventConfig.bFieldZ};

  const EventSummary summary = summarize(event, truthPt * 1_MeV / 1_GeV);
  std::cout << "layers=" << layout.layers.size()
            << " etaBins=" << geometry->numBins()
            << " binGroups=" << geometry->binGroups().size() << "\n"
            << "spacePoints=" << summary.spacePoints
            << " primaryHits=" << summary.primaryHits
            << " secondaryHits=" << summary.secondaryHits
            << " primaries=" << summary.primaries
            << " secondaries=" << summary.secondaries
            << " seedable=" << summary.seedablePrimaries << std::endl;

  // matched outside the timed region, so that the truth lookup does not show
  // up in the measurement
  Acts::SeedContainer reference;
  reference.assignSpacePointContainer(spacePoints);
  seeder.createSeeds(spacePoints, roi, graphBuilder, filter, options,
                     reference);
  const SeedingSummary seedSummary =
      evaluateSeeds(event, reference, truthPt * 1_MeV / 1_GeV);

  std::cout << "seeds=" << seedSummary.seeds
            << " trueSeeds=" << seedSummary.trueSeeds << " efficiency="
            << static_cast<float>(seedSummary.matchedPrimaries) /
                   static_cast<float>(
                       std::max<std::size_t>(1, summary.seedablePrimaries))
            << std::endl;

  // What an experiment calls: node building, graph and seed extraction in one.
  const auto full = microBenchmark(
      [&] {
        Acts::SeedContainer seeds;
        seeds.assignSpacePointContainer(spacePoints);
        seeder.createSeeds(spacePoints, roi, graphBuilder, filter, options,
                           seeds);
        assumeRead(seeds);
      },
      1, numRuns);
  std::cout << "full: " << full << std::endl;

  // The node storage fill alone, which copies every space point into the
  // storage and bins it. The graph stage is what the full call adds to it.
  const auto layerColumn =
      spacePoints.column<Exp::GbtsLayerIndex>("gbtsLayerIndex");
  const auto clusterWidthColumn = spacePoints.column<float>("clusterWidth");
  const auto localPositionColumn = spacePoints.column<float>("localPositionY");
  const auto nodes = microBenchmark(
      [&] {
        Exp::GbtsNodeStorage nodeStorage = seeder.makeNodeStorage();
        nodeStorage.extend(spacePoints, layerColumn, clusterWidthColumn,
                           localPositionColumn);
        nodeStorage.finalize();
        assumeRead(nodeStorage.numberOfNodes());
      },
      1, numRuns);
  std::cout << "nodes only: " << nodes << std::endl;

  return 0;
}
