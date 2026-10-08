// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Throughput benchmark for the grid triplet seeder on an ITk-like synthetic
/// event.
///
/// The detector and the event come from `ActsFatras::Synthetic` via
/// SyntheticEventOptions.hpp and the triplet cuts from
/// TripletSeedingConfig.hpp, both shared with the other seeding benchmarks.
/// What is added here is the grid and the driving loop, taken from how ATLAS
/// runs this seeder on the ITk pixel detector:
/// `ActsTrk::GridTripletSeedingTool` with the overrides
/// `ActsPixelSeedingToolCfg` applies, at the ITk main pass values of the
/// tracking flags it reads (900 MeV, 5 mm).
///
/// The same loop runs on either space point grid, so that they can be compared
/// on equal terms:
///
///  - `--grid cylindrical` bins in (phi, z, r), which is what ATLAS runs.
///  - `--grid spherical` bins in (phi, eta, r), whose middle axis measures a
///    direction rather than a position and so groups space points the way the
///    doublet cuts do. What it costs is that cot(theta) = z / r is the eta seen
///    from the origin and not the one of the track, so a vertex displaced by z0
///    moves a space point at radius r by z0 / r; `makeNeighbors` sizes the
///    neighbour windows from that offset.
///
/// Everything the experiment does per event is inside the benchmarked region:
/// filling the grid, sorting each bin in radius, copying into the packed
/// container the seeder consumes, and the seeding itself. `bottomPairs` and
/// `topPairs` report the space point pairs the binning lets through to the
/// doublet finder and `ranges` how often it is entered, so that a tuning run
/// can be read without a profiler.
///
/// @note Efficiency does not separate the two grids on this event: a forward
///       primary leaves fourteen space points, so finding one true seed is easy
///       for anything. Read the time and the candidate pair counts.

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding/CylindricalSpacePointGrid.hpp"
#include "Acts/Seeding/DoubletSeedFinder.hpp"
#include "Acts/Seeding/SphericalSpacePointGrid.hpp"
#include "Acts/Seeding/TripletSeedFinder.hpp"
#include "Acts/Seeding/TripletSeeder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/EventCsvWriter.hpp"
#include "ActsFatras/Synthetic/EventGenerator.hpp"
#include "ActsFatras/Synthetic/SeedingTruth.hpp"
#include "ActsTests/CommonHelpers/BenchmarkTools.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
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

/// The (phi, eta, r) grid, which has no ATLAS counterpart to copy from. Only
/// the middle axis differs from the cylindrical setup; everything else is
/// shared. Each value below was scanned on the synthetic event and sits where
/// efficiency stops paying for the time it costs.
struct SphericalGridConfig {
  /// Pseudorapidity acceptance, matching `cotThetaMax`
  float etaMin = -4.f;
  float etaMax = 4.f;
  /// Bin width on the middle axis, in cot(theta) rather than in eta. Set by the
  /// offset below and not by a resolution: a window has to reach that far
  /// regardless, so binning finer only enters the doublet finder more often
  /// without rejecting anything more.
  float cotBinWidth = 0.6f;
  /// Widening of the bins proportional to |cot(theta)|, giving the forward
  /// region a relative rather than an absolute resolution. Off: every value
  /// tried cost more in extra candidate pairs than it saved in bins.
  float cotBinGrowth = 0.f;
  /// Radial bin edges in mm. A single bin: splitting the range does cut
  /// candidate pairs, but a bin costs the doublet finder the same whether it
  /// holds one space point or a hundred, and the extra entries came out more
  /// expensive than the pairs they saved.
  std::vector<float> rBinEdges{0.f, 320.f};
  /// How far a bottom or top space point of the same track may sit from the
  /// middle one on the cot(theta) axis, which `makeNeighbors` turns into
  /// per-bin windows. A track from z0 crosses radius r at cot(theta) =
  /// cot(theta)_track + z0 / r, so two of its space points are at most |z0| *
  /// (1 / r_inner - 1 / r_outer) apart. The bottoms want about twice the reach
  /// of the tops, sitting at the small radii where the offset is largest.
  float cotSpreadBottom = 1.5f;
  float cotSpreadTop = 0.6f;
  /// Radial range middle candidates are taken from, as in `rRangeMiddleSP`
  float middleRMin = 40.f;
  float middleRMax = 260.f;
  /// Longitudinal reach of middle candidates in mm, i.e. how ATLAS avoids
  /// seeding from the last disc where nothing can sit beyond the middle space
  /// point. A bin here holds a range of cot(theta) rather than of z, so
  /// `makeSphericalMiddleAxis` turns this into a per-bin cap on the radius.
  float middleZMax = 2700.f;
  /// Longitudinal reach of a doublet in mm. The pixel configuration switches
  /// this off because a cylindrical grid's z bin windows already bound it; a
  /// spherical grid has no z bins, so it has to be made explicitly again.
  /// Without it the seeder returns 55k seeds instead of 33k for the same 7k
  /// true ones. The ATLAS default of 600 mm costs efficiency here.
  float deltaZMax = 900.f;
};

/// Turn a reach along an axis into the per-bin neighbour windows the grid bin
/// finder wants: the smallest window covering the bin itself grown by `below`
/// and `above`. The bins of the spherical grid are not uniform, so the reach is
/// stated once in the coordinate of the axis.
/// @param edges the bin edges of the axis
/// @param below how far below the lower edge of a bin to reach
/// @param above how far above the upper edge of a bin to reach
/// @return one window per bin
std::vector<std::pair<int, int>> makeNeighbors(const std::vector<float>& edges,
                                               float below, float above) {
  const int numBins = static_cast<int>(edges.size()) - 1;
  std::vector<std::pair<int, int>> neighbors(numBins);
  for (int bin = 0; bin < numBins; ++bin) {
    int down = 0;
    while (bin - down > 0 && edges[bin] - below < edges[bin - down]) {
      ++down;
    }
    int up = 0;
    while (bin + up + 1 < numBins &&
           edges[bin + 1] + above > edges[bin + up + 1]) {
      ++up;
    }
    neighbors[bin] = {-down, up};
  }
  return neighbors;
}

/// The binning of the middle axis, in whichever coordinate the grid bins it.
/// The driving loop uses it to look up the radial range middle candidates are
/// taken from without having to know which grid it is running on.
struct MiddleAxis {
  /// Bin edges, in mm for the cylindrical grid and in cot(theta) for the
  /// spherical one
  std::vector<float> edges;
  /// Radial range of middle candidates, one entry per bin
  std::vector<std::pair<float, float>> rRanges;
};

Acts::CylindricalSpacePointGrid::Config makeCylindricalGridConfig(
    const ItkPixelConfig& cfg) {
  Acts::CylindricalSpacePointGrid::Config gridCfg;
  gridCfg.minPt = cfg.minPt;
  gridCfg.rMin = 0;
  gridCfg.rMax = cfg.gridRMax;
  gridCfg.zMin = cfg.zMin;
  gridCfg.zMax = cfg.zMax;
  gridCfg.deltaRMax = cfg.deltaRMax;
  gridCfg.cotThetaMax = cfg.cotThetaMax;
  gridCfg.impactMax = cfg.impactMax;
  gridCfg.phiMin = cfg.phiMin;
  gridCfg.phiMax = cfg.phiMax;
  gridCfg.phiBinDeflectionCoverage = cfg.phiBinDeflectionCoverage;
  gridCfg.maxPhiBins = cfg.maxPhiBins;
  gridCfg.zBinEdges = cfg.zBinEdges;
  gridCfg.rBinEdges = cfg.rBinEdges;
  gridCfg.bFieldInZ = cfg.bFieldInZ;
  gridCfg.bottomBinFinder.emplace(cfg.numPhiNeighbors, cfg.zBinNeighborsBottom,
                                  cfg.rBinNeighborsBottom);
  gridCfg.topBinFinder.emplace(cfg.numPhiNeighbors, cfg.zBinNeighborsTop,
                               cfg.rBinNeighborsTop);
  gridCfg.navigation[0ul] = {};
  gridCfg.navigation[1ul] = cfg.zBinsCustomLooping;
  gridCfg.navigation[2ul] = cfg.rBinsCustomLooping;
  return gridCfg;
}

/// The middle axis of the cylindrical grid: its z bins, with the radial range
/// ATLAS takes middle candidates from in each of them.
/// @param cfg the seeding configuration
/// @return the middle axis binning
MiddleAxis makeCylindricalMiddleAxis(const ItkPixelConfig& cfg) {
  return {cfg.zBinEdges, cfg.rRangeMiddleSP};
}

/// The bin edges of the middle axis of the spherical grid, in cot(theta), which
/// is the coordinate space points are binned in. Bins may grow with
/// |cot(theta)|, so they are built outwards from the centre in both directions.
/// @param sph the spherical grid configuration
/// @return the bin edges in cot(theta)
std::vector<float> makeCotThetaBinEdges(const SphericalGridConfig& sph) {
  std::vector<float> positive{0.f};
  const float cotMax = std::sinh(std::max(-sph.etaMin, sph.etaMax));
  while (positive.back() < cotMax) {
    positive.push_back(positive.back() + sph.cotBinWidth +
                       sph.cotBinGrowth * positive.back());
  }

  std::vector<float> edges;
  edges.reserve(2 * positive.size() - 1);
  for (auto it = positive.rbegin(); it != positive.rend(); ++it) {
    edges.push_back(-*it);
  }
  edges.insert(edges.end(), std::next(positive.begin()), positive.end());
  return edges;
}

/// The same edges in pseudorapidity, which is how the grid takes them.
/// @param sph the spherical grid configuration
/// @return the bin edges in eta
std::vector<float> makeEtaBinEdges(const SphericalGridConfig& sph) {
  std::vector<float> edges = makeCotThetaBinEdges(sph);
  for (float& edge : edges) {
    edge = std::asinh(edge);
  }
  return edges;
}

Exp::SphericalSpacePointGrid::Config makeSphericalGridConfig(
    const ItkPixelConfig& cfg, const SphericalGridConfig& sph) {
  Exp::SphericalSpacePointGrid::Config gridCfg;
  gridCfg.minPt = cfg.minPt;
  gridCfg.rMin = 0;
  gridCfg.rMax = cfg.gridRMax;
  gridCfg.etaMin = sph.etaMin;
  gridCfg.etaMax = sph.etaMax;
  gridCfg.etaBinEdges = makeEtaBinEdges(sph);
  gridCfg.deltaRMax = cfg.deltaRMax;
  gridCfg.impactMax = cfg.impactMax;
  gridCfg.phiMin = cfg.phiMin;
  gridCfg.phiMax = cfg.phiMax;
  gridCfg.phiBinDeflectionCoverage = cfg.phiBinDeflectionCoverage;
  gridCfg.maxPhiBins = cfg.maxPhiBins;
  gridCfg.rBinEdges = sph.rBinEdges;
  gridCfg.bFieldInZ = cfg.bFieldInZ;

  // both axes take their windows from the reach of what looks for them: the
  // radial one from how far the doublet finders go, the middle one from how far
  // the beam spot can move a space point along it
  const std::vector<float> cotThetaEdges = makeCotThetaBinEdges(sph);
  gridCfg.bottomBinFinder.emplace(
      cfg.numPhiNeighbors,
      makeNeighbors(cotThetaEdges, sph.cotSpreadBottom, sph.cotSpreadBottom),
      makeNeighbors(sph.rBinEdges, cfg.deltaRMaxBottomSP, 0.f));
  gridCfg.topBinFinder.emplace(
      cfg.numPhiNeighbors,
      makeNeighbors(cotThetaEdges, sph.cotSpreadTop, sph.cotSpreadTop),
      makeNeighbors(sph.rBinEdges, 0.f, cfg.deltaRMaxTopSP));
  gridCfg.navigation[0ul] = {};
  gridCfg.navigation[1ul] = {};
  gridCfg.navigation[2ul] = {};
  return gridCfg;
}

/// The middle axis of the spherical grid. Every bin offers the same range,
/// except that the more forward ones have it cut back so that no middle
/// candidate sits beyond `middleZMax`.
/// @param sph the spherical grid configuration
/// @return the middle axis binning
MiddleAxis makeSphericalMiddleAxis(const SphericalGridConfig& sph) {
  const std::vector<float> edges = makeCotThetaBinEdges(sph);
  std::vector<std::pair<float, float>> ranges;
  ranges.reserve(edges.size() - 1);
  for (std::size_t bin = 0; bin + 1 < edges.size(); ++bin) {
    const float cotMax =
        std::max(std::abs(edges[bin]), std::abs(edges[bin + 1]));
    ranges.emplace_back(sph.middleRMin,
                        std::min(sph.middleRMax, sph.middleZMax / cotMax));
  }
  return {edges, ranges};
}

/// The radial range middle space point candidates are taken from, looked up in
/// the middle axis bin the candidates sit in.
/// @param cfg the seeding configuration
/// @param axis the binning of the middle axis
/// @param coordinate the middle axis coordinate of the middle space points
/// @param variableRange the range derived from the extent of the event
/// @return the radial range in mm
std::pair<float, float> radiusRangeForMiddle(
    const ItkPixelConfig& cfg, const MiddleAxis& axis, float coordinate,
    const Acts::Range1D<float>& variableRange) {
  if (cfg.useVariableMiddleSPRange) {
    return {variableRange.min(), variableRange.max()};
  }
  auto edge = std::ranges::lower_bound(axis.edges, coordinate);
  std::size_t bin = std::distance(axis.edges.begin(), edge);
  if (bin > 0) {
    --bin;
  }
  return axis.rRanges.at(std::min(bin, axis.rRanges.size() - 1));
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
  std::unique_ptr<const Acts::Logger> logger;
};

/// Run one full seeding pass over an event, the way the ATLAS tool does.
/// @tparam GridType the space point grid to bin the event in
/// @param cfg the seeding configuration
/// @param gridConfig the configuration of the grid
/// @param axis the binning of the middle axis of the grid
/// @param cache the reusable finders and their state
/// @param spacePoints the space points of the event
/// @param seeds receives the seeds
template <typename GridType>
void createSeeds(const ItkPixelConfig& cfg,
                 const typename GridType::Config& gridConfig,
                 const MiddleAxis& axis, SeedingCache& cache,
                 const Acts::SpacePointContainer& spacePoints,
                 Acts::SeedContainer& seeds) {
  constexpr bool spherical =
      std::is_same_v<GridType, Exp::SphericalSpacePointGrid>;

  GridType grid(gridConfig, cache.logger->clone());

  for (std::size_t i = 0; i < spacePoints.size(); ++i) {
    const auto& sp = spacePoints[i];
    if constexpr (spherical) {
      grid.insert(i, sp.phi(), sp.z() / sp.r(), sp.r());
    } else {
      grid.insert(i, sp.phi(), sp.z(), sp.r());
    }
  }
  // Not `grid.sortBinsByR`: that reads the packed ZR column, and the generated
  // container carries Z and R separately, as an experiment's own space point
  // EDM typically does before it is repacked below.
  for (std::size_t i = 0; i < grid.numberOfBins(); ++i) {
    std::ranges::sort(grid.at(i),
                      [&](Acts::SpacePointIndex a, Acts::SpacePointIndex b) {
                        return spacePoints[a].r() < spacePoints[b].r();
                      });
  }

  // the seeder works on a packed container laid out in grid bin order
  Acts::SpacePointContainer gridSpacePoints(
      Acts::SpacePointColumns::CopiedFromIndex |
      Acts::SpacePointColumns::PackedXY | Acts::SpacePointColumns::PackedZR |
      Acts::SpacePointColumns::VarianceZ | Acts::SpacePointColumns::VarianceR);
  gridSpacePoints.reserve(grid.numberOfSpacePoints());
  std::vector<Acts::SpacePointIndexRange> gridSpacePointRanges;
  gridSpacePointRanges.reserve(grid.numberOfBins());
  for (std::size_t i = 0; i < grid.numberOfBins(); ++i) {
    const std::uint32_t begin = gridSpacePoints.size();
    for (const Acts::SpacePointIndex spIndex : grid.at(i)) {
      const auto& sp = spacePoints[spIndex];
      auto newSp = gridSpacePoints.createSpacePoint();
      newSp.copiedFromIndex() = sp.index();
      newSp.xy() = std::array<float, 2>{sp.x(), sp.y()};
      newSp.zr() = std::array<float, 2>{sp.z(), sp.r()};
      newSp.varianceZ() = sp.varianceZ();
      newSp.varianceR() = sp.varianceR();
    }
    gridSpacePointRanges.emplace_back(begin, gridSpacePoints.size());
  }

  // Radial extent of the event, exploiting that every bin is sorted in radius.
  // Not `grid.computeRadiusRange`: as with the sort above, and because the grid
  // indexes the container it was filled from rather than this reordered copy.
  const Acts::Range1D<float> rRange = [&]() -> Acts::Range1D<float> {
    float minRange = std::numeric_limits<float>::max();
    float maxRange = std::numeric_limits<float>::lowest();
    for (const Acts::SpacePointIndexRange& range : gridSpacePointRanges) {
      if (range.first == range.second) {
        continue;
      }
      minRange = std::min(gridSpacePoints[range.first].zr()[1], minRange);
      maxRange = std::max(gridSpacePoints[range.second - 1].zr()[1], maxRange);
    }
    return {minRange, maxRange};
  }();
  const Acts::Range1D<float> variableMiddleRange(
      std::floor(rRange.min() / 2) * 2 + cfg.deltaRMiddleMinSPRange,
      std::floor(rRange.max() / 2) * 2 - cfg.deltaRMiddleMaxSPRange);

  Acts::BroadTripletSeedFilter::State filterState;
  Acts::BroadTripletSeedFilter::Cache filterCache;
  Acts::BroadTripletSeedFilter filter(cache.filterConfig, filterState,
                                      filterCache, *cache.logger);

  std::vector<Acts::SpacePointContainer::ConstRange> bottomSpRanges;
  std::vector<Acts::SpacePointContainer::ConstRange> topSpRanges;

  Acts::SeedContainer candidates;

  for (const auto [bottom, middle, top] : grid.binnedGroup()) {
    const Acts::SpacePointContainer::ConstRange middleSpRange =
        gridSpacePoints.range(gridSpacePointRanges[middle]).asConst();
    if (middleSpRange.empty()) {
      continue;
    }

    bottomSpRanges.clear();
    for (const std::size_t b : bottom) {
      bottomSpRanges.push_back(
          gridSpacePoints.range(gridSpacePointRanges[b]).asConst());
    }
    topSpRanges.clear();
    for (const std::size_t t : top) {
      topSpRanges.push_back(
          gridSpacePoints.range(gridSpacePointRanges[t]).asConst());
    }

    // all middle candidates of a group share a middle axis bin, so this is a
    // per-group lookup
    const std::array<float, 2> front = middleSpRange.front().zr();
    const std::pair<float, float> middleRange = radiusRangeForMiddle(
        cfg, axis, spherical ? front[0] / front[1] : front[0],
        variableMiddleRange);

    cache.seeder.createSeedsFromGroups(
        cache.seederCache, *cache.bottomDoubletFinder, *cache.topDoubletFinder,
        *cache.tripletFinder, filter, gridSpacePoints, bottomSpRanges,
        middleSpRange, topSpRanges, middleRange, candidates);
  }

  // Keep a seed only if it is the best one found for at least one of its three
  // space points. ATLAS applies this on top of the ACTS filter.
  seeds.assignSpacePointContainer(spacePoints);
  seeds.reserve(candidates.size());
  for (Acts::MutableSeedProxy candidate : candidates) {
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
        gridSpacePoints.at(indices[0]).copiedFromIndex(),
        gridSpacePoints.at(indices[1]).copiedFromIndex(),
        gridSpacePoints.at(indices[2]).copiedFromIndex()};
    auto seed = seeds.createSeed();
    seed.assignSpacePointIndices(original);
    seed.quality() = candidate.quality();
    seed.vertexZ() = candidate.vertexZ();
  }
}

/// Count the space point pairs the grid offers to the doublet finders, i.e. the
/// combinatorial load a grid configuration imposes. Everything downstream only
/// ever sees pairs that were offered here.
/// @tparam GridType the space point grid to bin the event in
/// @param cfg the seeding configuration
/// @param gridConfig the configuration of the grid
/// @param axis the binning of the middle axis of the grid
/// @param spacePoints the space points of the event
/// @param logger the logger to build the grid with
template <typename GridType>
void printGridStatistics(const ItkPixelConfig& cfg,
                         const typename GridType::Config& gridConfig,
                         const MiddleAxis& axis,
                         const Acts::SpacePointContainer& spacePoints,
                         const Acts::Logger& logger) {
  constexpr bool spherical =
      std::is_same_v<GridType, Exp::SphericalSpacePointGrid>;

  GridType grid(gridConfig, logger.clone());
  for (std::size_t i = 0; i < spacePoints.size(); ++i) {
    const auto& sp = spacePoints[i];
    if constexpr (spherical) {
      grid.insert(i, sp.phi(), sp.z() / sp.r(), sp.r());
    } else {
      grid.insert(i, sp.phi(), sp.z(), sp.r());
    }
  }

  std::size_t occupied = 0;
  for (std::size_t i = 0; i < grid.numberOfBins(); ++i) {
    occupied += grid.at(i).empty() ? 0 : 1;
  }

  std::size_t middles = 0;
  std::size_t bottomPairs = 0;
  std::size_t topPairs = 0;
  // the doublet finder is entered once per middle space point and neighbour
  // bin, so this counts the fixed cost the binning imposes independently of how
  // many space points the bins hold
  std::size_t ranges = 0;
  for (const auto [bottom, middle, top] : grid.binnedGroup()) {
    std::size_t numMiddle = 0;
    for (const Acts::SpacePointIndex index : grid.at(middle)) {
      const auto& sp = spacePoints[index];
      const std::pair<float, float> range = radiusRangeForMiddle(
          cfg, axis, spherical ? sp.z() / sp.r() : sp.z(), {0, 0});
      numMiddle += (sp.r() >= range.first && sp.r() <= range.second) ? 1 : 0;
    }
    middles += numMiddle;
    ranges += numMiddle * (bottom.size() + top.size());
    for (const std::size_t b : bottom) {
      bottomPairs += numMiddle * grid.at(b).size();
    }
    for (const std::size_t t : top) {
      topPairs += numMiddle * grid.at(t).size();
    }
  }

  std::cout << "bins=" << grid.numberOfBins() << " occupied=" << occupied
            << " middles=" << middles << " bottomPairs=" << bottomPairs
            << " topPairs=" << topPairs << " ranges=" << ranges << std::endl;
}

}  // namespace

int main(int argc, char* argv[]) {
  SyntheticEventOptions::Values shared;
  std::size_t numRuns = 10;
  // Number of distinct events the seeder is timed over. One keeps the cache hot
  // the way a repeated measurement does; more averages over the event-to-event
  // spread in multiplicity at the cost of a larger working set.
  std::size_t numEvents = 1;
  float minPt = ItkPixelConfig{}.minPt / 1_MeV;
  // Efficiency is counted over a harder threshold than the seeder is cut at, so
  // that the turn-on stays out of it.
  float truthPt = 1000.f;
  std::string gridName = "cylindrical";
  std::string dumpPrefix;
  std::vector<float> middleRange;
  bool allMiddleBins = false;
  bool forceDeltaZMax = false;
  bool verbose = false;

  SphericalGridConfig sph;

  try {
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message");
    SyntheticEventOptions::add(desc, shared);
    desc.add_options()(
        "grid", po::value<std::string>(&gridName)->default_value(gridName),
        "space point grid to bin the event in, cylindrical or spherical")(
        "runs", po::value<std::size_t>(&numRuns)->default_value(numRuns),
        "number of benchmark runs")(
        "events", po::value<std::size_t>(&numEvents)->default_value(numEvents),
        "number of distinct events to generate and cycle through")(
        "min-pt", po::value<float>(&minPt)->default_value(minPt),
        "seed momentum threshold in MeV")(
        "truth-pt", po::value<float>(&truthPt)->default_value(truthPt),
        "momentum in MeV a primary has to reach to be counted in the "
        "efficiency")(
        "cot-bin-width",
        po::value<float>(&sph.cotBinWidth)->default_value(sph.cotBinWidth),
        "spherical grid: width of the middle axis bins in cot(theta)")(
        "cot-bin-growth",
        po::value<float>(&sph.cotBinGrowth)->default_value(sph.cotBinGrowth),
        "spherical grid: widening of the bins per unit of |cot(theta)|")(
        "cot-spread-bottom",
        po::value<float>(&sph.cotSpreadBottom)
            ->default_value(sph.cotSpreadBottom),
        "spherical grid: cot(theta) reach towards bottom candidates")(
        "cot-spread-top",
        po::value<float>(&sph.cotSpreadTop)->default_value(sph.cotSpreadTop),
        "spherical grid: cot(theta) reach towards top candidates")(
        "delta-z-max",
        po::value<float>(&sph.deltaZMax)->default_value(sph.deltaZMax),
        "spherical grid: longitudinal reach of a doublet in mm")(
        "r-bins", po::value<std::vector<float>>(&sph.rBinEdges)->multitoken(),
        "spherical grid: radial bin edges in mm")(
        "middle-range",
        po::value<std::vector<float>>(&middleRange)->multitoken(),
        "override the radial range of middle candidates in every bin of either "
        "grid, so that both see the same middle space points")(
        "all-middle-bins", po::bool_switch(&allMiddleBins),
        "cylindrical grid: drop the custom looping order, which otherwise "
        "leaves the two outermost z bins out of the middle bins entirely")(
        "dump", po::value<std::string>(&dumpPrefix),
        "write the generated event to <arg>_spacepoints.csv and "
        "<arg>_particles.csv")("verbose", po::bool_switch(&verbose),
                               "log the seeder's own statistics");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help") > 0) {
      std::cout << desc << std::endl;
      return 0;
    }
    // the cut belongs to the spherical setup, but asking for it explicitly
    // applies it to either grid, so that the two can be told apart from it
    forceDeltaZMax = !vm["delta-z-max"].defaulted();
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }

  if (gridName != "cylindrical" && gridName != "spherical") {
    std::cerr << "error: unknown grid '" << gridName << "'" << std::endl;
    return 1;
  }
  const bool spherical = gridName == "spherical";

  // the middle axis is walked outwards one bin at a time, so a width that does
  // not advance would never terminate
  if (sph.cotBinWidth <= 0.f || sph.cotBinGrowth < 0.f) {
    std::cerr << "error: the middle axis bins have to grow" << std::endl;
    return 1;
  }
  if (sph.rBinEdges.size() < 2 || !std::ranges::is_sorted(sph.rBinEdges)) {
    std::cerr << "error: --r-bins takes at least two increasing edges"
              << std::endl;
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

  // Generated up front and then only read, so that no part of the generation
  // lands inside the timed region.
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
  // reaches the phi binning of the grid, both doublet finders and the triplet
  // finder, so it has to be set before any of them is built
  cfg.minPt = minPt * 1_MeV;
  if (spherical || forceDeltaZMax) {
    cfg.deltaZMax = sph.deltaZMax;
  }

  SeedingCache cache;
  cache.logger = Acts::getDefaultLogger(
      "GridTriplet",
      verbose ? Acts::Logging::Level::DEBUG : Acts::Logging::Level::FATAL);
  cache.filterConfig = makeFilterConfig(cfg);
  const Acts::DoubletSeedFinder::Config bottomConfig =
      makeBottomDoubletConfig(cfg);
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

  if (allMiddleBins) {
    cfg.zBinsCustomLooping.clear();
  }
  const Acts::CylindricalSpacePointGrid::Config cylindricalConfig =
      makeCylindricalGridConfig(cfg);
  const Exp::SphericalSpacePointGrid::Config sphericalConfig =
      makeSphericalGridConfig(cfg, sph);
  MiddleAxis axis =
      spherical ? makeSphericalMiddleAxis(sph) : makeCylindricalMiddleAxis(cfg);
  if (middleRange.size() == 2) {
    std::ranges::fill(axis.rRanges,
                      std::pair<float, float>{middleRange[0], middleRange[1]});
  } else if (!middleRange.empty()) {
    std::cerr << "error: --middle-range takes two values" << std::endl;
    return 1;
  }

  // dispatch once, so that the grid type is a compile time property of the
  // seeding pass and not a branch inside it
  const auto seed = [&](const Event& event, Acts::SeedContainer& seeds) {
    if (spherical) {
      createSeeds<Exp::SphericalSpacePointGrid>(
          cfg, sphericalConfig, axis, cache, event.spacePoints, seeds);
    } else {
      createSeeds<Acts::CylindricalSpacePointGrid>(
          cfg, cylindricalConfig, axis, cache, event.spacePoints, seeds);
    }
  };

  if (spherical) {
    printGridStatistics<Exp::SphericalSpacePointGrid>(
        cfg, sphericalConfig, axis, events.front().spacePoints, *cache.logger);
  } else {
    printGridStatistics<Acts::CylindricalSpacePointGrid>(
        cfg, cylindricalConfig, axis, events.front().spacePoints,
        *cache.logger);
  }

  // Truth and efficiency are accumulated over every event, and outside the
  // timed region so that the truth lookup does not show up in the measurement.
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
    seed(event, reference);
    const SeedingSummary oneSeeding =
        evaluateSeeds(event, reference, truthPt * 1_MeV / 1_GeV);
    seedSummary.seeds += oneSeeding.seeds;
    seedSummary.trueSeeds += oneSeeding.trueSeeds;
    seedSummary.matchedPrimaries += oneSeeding.matchedPrimaries;
  }

  const auto perEvent = [&](std::size_t value) { return value / numEvents; };
  std::cout << "grid=" << gridName << " events=" << numEvents
            << "\nmean/event: spacePoints=" << perEvent(summary.spacePoints)
            << " primaryHits=" << perEvent(summary.primaryHits)
            << " secondaryHits=" << perEvent(summary.secondaryHits)
            << " primaries=" << perEvent(summary.primaries)
            << " secondaries=" << perEvent(summary.secondaries)
            << " seedable=" << perEvent(summary.seedablePrimaries) << std::endl;

  // Each timed run seeds one event, cycling through them, so the reported
  // median is per event whatever `--events` is.
  std::size_t next = 0;
  const auto result = microBenchmark(
      [&] {
        Acts::SeedContainer seeds;
        seed(events[next], seeds);
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
