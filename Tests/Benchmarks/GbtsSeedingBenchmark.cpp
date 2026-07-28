// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Throughput benchmark for the GBTS seeder on an ITk-like synthetic event.
///
/// The detector and the event come from SyntheticItkEvent.hpp, which is shared
/// with the other seeding benchmarks. Only the GBTS-specific geometry, namely
/// the logical layer numbering and the layer connection table, is built here.

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsLayerConnection.hpp"
#include "Acts/Seeding/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsTests/CommonHelpers/BenchmarkTools.hpp"

#include <cstdint>
#include <iostream>
#include <memory>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "SyntheticItkEvent.hpp"

namespace po = boost::program_options;
namespace Exp = Acts::Experimental;

using namespace Acts::UnitLiterals;
using namespace ActsTests;
using namespace ActsTests::SyntheticItk;

namespace {

/// GBTS numbers its logical layers as `gbtsId * 1000 + etaModule`, with the
/// barrel in the eighties, the positive endcap in the nineties and the negative
/// endcap in the seventies. GraphBasedTrackSeeder relies on that encoding: it
/// derives the barrel flag as `layerId / 10000 == 8`.
std::int32_t gbtsLayerId(const Layer& layer) {
  std::int32_t group = 0;
  if (layer.type == LayerType::Barrel) {
    group = 80 + layer.layer;
  } else if (layer.side > 0) {
    group = 90 + layer.layer;
  } else {
    group = 78 - layer.layer;
  }
  return group * 1000 + layer.module;
}

std::vector<Exp::GbtsLayerDescription> makeLayerDescriptions(
    const DetectorLayout& layout) {
  std::vector<Exp::GbtsLayerDescription> descriptions;
  descriptions.reserve(layout.layers.size());
  for (const Layer& layer : layout.layers) {
    Exp::GbtsLayerDescription description;
    description.id = gbtsLayerId(layer);
    description.type = layer.type == LayerType::Barrel
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
std::string makeConnectionTable(const DetectorLayout& layout,
                                float etaBinWidth) {
  // group layers by (side, layer) so that all modules of a surface are linked
  auto modulesOf = [&](int side, int layerIdx) {
    std::vector<std::int32_t> ids;
    for (const Layer& layer : layout.layers) {
      if (layer.side == side && layer.layer == layerIdx) {
        ids.push_back(gbtsLayerId(layer));
      }
    }
    return ids;
  };

  std::vector<std::pair<std::int32_t, std::int32_t>> connections;
  auto connect = [&](int outerSide, int outerLayer, int innerSide,
                     int innerLayer) {
    for (const std::int32_t src : modulesOf(outerSide, outerLayer)) {
      for (const std::int32_t dst : modulesOf(innerSide, innerLayer)) {
        connections.emplace_back(src, dst);
      }
    }
  };

  int numBarrel = 0;
  int numDisks = 0;
  for (const Layer& layer : layout.layers) {
    if (layer.type == LayerType::Barrel) {
      numBarrel = std::max(numBarrel, layer.layer + 1);
    } else if (layer.side > 0) {
      numDisks = std::max(numDisks, layer.layer + 1);
    }
  }

  // barrel to barrel, adjacent and one skipped layer
  for (int i = 0; i < numBarrel; ++i) {
    for (int step : {1, 2}) {
      if (i + step < numBarrel) {
        connect(0, i + step, 0, i);
      }
    }
  }
  for (int side : {+1, -1}) {
    // disk to disk, adjacent and one skipped disk
    for (int j = 0; j < numDisks; ++j) {
      for (int step : {1, 2}) {
        if (j + step < numDisks) {
          connect(side, j + step, side, j);
        }
      }
    }
    // forward transition: the innermost barrel layers are fed by the first
    // disks of each side
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        connect(side, j, 0, i);
      }
    }
  }

  std::ostringstream os;
  os << connections.size() << " " << etaBinWidth << "\n";
  for (std::size_t k = 0; k < connections.size(); ++k) {
    // index stage src dst height width entries
    os << k << " " << k << " " << connections[k].first << " "
       << connections[k].second << " 0 0 0\n";
  }
  return os.str();
}

}  // namespace

int main(int argc, char* argv[]) {
  std::size_t numTracks = EventConfig{}.numTracks;
  std::size_t numRuns = 10;
  bool verbose = false;
  // The synthetic event is noise free, so more of the geometric doublet
  // candidates survive the cuts than in a real event and the graph outgrows the
  // 2000000 edges that GraphBasedTrackSeeder allows by default. Benchmarking
  // against a saturated ceiling would hide any work done past it, so the
  // default here is raised until the graph fits. Pass the production value to
  // measure the truncated regime instead.
  std::uint32_t maxEdges = 6000000;

  try {
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")(
        "tracks", po::value<std::size_t>(&numTracks)->default_value(numTracks),
        "number of generated tracks")(
        "runs", po::value<std::size_t>(&numRuns)->default_value(numRuns),
        "number of benchmark runs")(
        "max-edges",
        po::value<std::uint32_t>(&maxEdges)->default_value(maxEdges),
        "ceiling on the number of graph edges")(
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

  const DetectorLayout layout = makePixelLayout();

  EventConfig eventConfig;
  eventConfig.numTracks = numTracks;
  const Acts::SpacePointContainer spacePoints =
      generateEvent(layout, eventConfig);

  constexpr float etaBinWidth = 0.2f;
  const std::string table = makeConnectionTable(layout, etaBinWidth);
  std::istringstream tableStream(table);
  const Exp::GbtsLayerConnectionMap connections =
      Exp::GbtsLayerConnectionMap::fromStream(tableStream, false);
  const auto geometry = std::make_shared<Exp::GbtsGeometry>(
      makeLayerDescriptions(layout), connections);

  Exp::GraphBasedTrackSeeder::Config cfg;
  cfg.minPt = 1_GeV;
  cfg.minZ0 = -200.f;
  cfg.maxZ0 = 200.f;
  cfg.maxOuterRadius = 550.f;
  cfg.useMl = false;
  cfg.nMaxEdges = maxEdges;

  const Exp::GraphBasedTrackSeeder::DerivedConfig derived(cfg);
  const Exp::GraphBasedTrackSeeder seeder(
      derived, geometry,
      Acts::getDefaultLogger("Gbts", verbose ? Acts::Logging::Level::DEBUG
                                             : Acts::Logging::Level::FATAL));
  const Exp::GbtsTrackingFilter filter(Exp::GbtsTrackingFilter::Config{},
                                       geometry);
  const Exp::GbtsRoiDescriptor roi(0, -5., 5., 0, -std::numbers::pi,
                                   std::numbers::pi, 0, cfg.minZ0, cfg.maxZ0);
  const Exp::GraphBasedTrackSeeder::Options options(eventConfig.bFieldZ * 1_T);
  const std::vector<bool> isPixelLayer(layout.layers.size(), true);

  std::size_t numSeeds = 0;

  std::cout << "layers=" << layout.layers.size()
            << " etaBins=" << geometry->numBins()
            << " binGroups=" << geometry->binGroups().size()
            << " tracks=" << numTracks << " spacePoints=" << spacePoints.size()
            << std::endl;

  const auto result = microBenchmark(
      [&] {
        Acts::SeedContainer seeds;
        seeds.assignSpacePointContainer(spacePoints);
        seeder.createSeeds(spacePoints, roi, isPixelLayer, filter, options,
                           seeds);
        numSeeds = seeds.size();
        assumeRead(seeds);
      },
      1, numRuns);

  std::cout << "seeds=" << numSeeds << " " << result << std::endl;

  return 0;
}
