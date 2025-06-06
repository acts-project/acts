// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TripletSeedingAlgorithm2.hpp"

#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding2/TripletSeedFinder2.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <cmath>
#include <csignal>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace ActsExamples {

TripletSeedingAlgorithm2::TripletSeedingAlgorithm2(const Config& cfg,
                                                   Acts::Logging::Level lvl)
    : IAlgorithm("SeedingAlgorithm2", lvl), m_cfg(cfg) {
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputSeeds.initialize(m_cfg.outputSeeds);

  static_assert(std::numeric_limits<
                    decltype(m_cfg.finderConfig.deltaRMaxTopSP)>::has_quiet_NaN,
                "Value of deltaRMaxTopSP must support NaN values");

  static_assert(std::numeric_limits<
                    decltype(m_cfg.finderConfig.deltaRMinTopSP)>::has_quiet_NaN,
                "Value of deltaRMinTopSP must support NaN values");

  static_assert(
      std::numeric_limits<
          decltype(m_cfg.finderConfig.deltaRMaxBottomSP)>::has_quiet_NaN,
      "Value of deltaRMaxBottomSP must support NaN values");

  static_assert(
      std::numeric_limits<
          decltype(m_cfg.finderConfig.deltaRMinBottomSP)>::has_quiet_NaN,
      "Value of deltaRMinBottomSP must support NaN values");

  if (m_cfg.gridConfig.cotThetaMax != m_cfg.finderConfig.cotThetaMax) {
    throw std::invalid_argument("Inconsistent config cotThetaMax");
  }

  if (m_cfg.gridConfig.minPt != m_cfg.finderConfig.minPt) {
    throw std::invalid_argument("Inconsistent config minPt");
  }

  if (m_cfg.gridConfig.bFieldInZ != m_cfg.finderOptions.bFieldInZ) {
    throw std::invalid_argument("Inconsistent config bFieldInZ");
  }

  if (m_cfg.gridConfig.zBinEdges.size() - 1 != m_cfg.zBinNeighborsTop.size() &&
      m_cfg.zBinNeighborsTop.empty() == false) {
    throw std::invalid_argument("Inconsistent config zBinNeighborsTop");
  }

  if (m_cfg.gridConfig.zBinEdges.size() - 1 !=
          m_cfg.zBinNeighborsBottom.size() &&
      m_cfg.zBinNeighborsBottom.empty() == false) {
    throw std::invalid_argument("Inconsistent config zBinNeighborsBottom");
  }

  if (!m_cfg.zBinsCustomLooping.empty()) {
    // check that the bins required in the custom bin looping
    // are contained in the bins defined by the total number of edges

    for (std::size_t i : m_cfg.zBinsCustomLooping) {
      if (i >= m_cfg.gridConfig.zBinEdges.size()) {
        throw std::invalid_argument(
            "Inconsistent config zBinsCustomLooping does not contain a subset "
            "of bins defined by zBinEdges");
      }
    }
  }

  if (m_cfg.useExtraCuts) {
    // This function will be applied to select space points during grid filling
    m_spacePointSelector.connect<itkFastTrackingSPselect>();

    // This function will be applied to the doublet compatibility selection
    m_cfg.finderConfig.experimentCuts.connect<itkFastTrackingCuts>();
  }

  m_cfg.gridConfig.bottomBinFinder.emplace(m_cfg.numPhiNeighbors,
                                           m_cfg.zBinNeighborsBottom, 0);
  m_cfg.gridConfig.topBinFinder.emplace(m_cfg.numPhiNeighbors,
                                        m_cfg.zBinNeighborsTop, 0);
  m_cfg.gridConfig.navigation[1ul] = m_cfg.zBinsCustomLooping;

  m_cfg.finderConfig.filter = std::make_unique<Acts::TripletSeedFilter2>(
      m_cfg.filterConfig, logger().cloneWithSuffix("Filter"));

  m_seedFinder = Acts::TripletSeedFinder2(m_cfg.finderConfig.derive(),
                                          logger().cloneWithSuffix("Finder"));
}

ProcessCode TripletSeedingAlgorithm2::execute(
    const AlgorithmContext& ctx) const {
  const SimSpacePointContainer& spacePoints = m_inputSpacePoints(ctx);

  Acts::SpacePointContainer2 coreSpacePoints;
  auto& rColumn = coreSpacePoints.createDenseExtraColumn<float>("r");
  auto& phiColumn = coreSpacePoints.createDenseExtraColumn<float>("phi");
  auto& varianceRColumn =
      coreSpacePoints.createDenseExtraColumn<float>("varianceR");
  auto& varianceZColumn =
      coreSpacePoints.createDenseExtraColumn<float>("varianceZ");
  coreSpacePoints.reserve(spacePoints.size());
  for (const auto& sp : spacePoints) {
    // check if the space point passes the selection
    if (m_spacePointSelector(sp)) {
      auto newSp = coreSpacePoints.createSpacePoint(
          std::array<Acts::SourceLink, 1>{Acts::SourceLink(&sp)}, sp.x(),
          sp.y(), sp.z());
      newSp.extra(rColumn) = sp.r();
      newSp.extra(phiColumn) = std::atan2(sp.y(), sp.x());
      newSp.extra(varianceRColumn) = sp.varianceR();
      newSp.extra(varianceZColumn) = sp.varianceZ();
    }
  }

  Acts::CylindricalSpacePointGrid2 grid(m_cfg.gridConfig.derive(),
                                        logger().cloneWithSuffix("Grid"));

  grid.fill(coreSpacePoints, phiColumn, rColumn);

  // Compute radius Range
  // we rely on the fact the grid is storing the proxies
  // with a sorting in the radius
  const Acts::Range1D<float> rRange =
      grid.computeRadiusRange(coreSpacePoints, rColumn);

  Acts::TripletSeedFinder2::Options finderOptions = m_cfg.finderOptions;

  /// variable middle SP radial region of interest
  finderOptions.rMiddleSpRange = {
      std::floor(rRange.min() / 2) * 2 + m_cfg.deltaRMiddleMinSPRange,
      std::floor(rRange.max() / 2) * 2 - m_cfg.deltaRMiddleMaxSPRange};

  // run the seeding
  Acts::SeedContainer2 seeds;
  Acts::TripletSeedFinder2::State state;
  static thread_local Acts::TripletSeedFinder2::Cache cache;

  auto derivedOptions = finderOptions.derive(m_seedFinder->config());
  m_seedFinder->initialize(state, cache, derivedOptions);

  std::vector<Acts::SpacePointIndex2> bottomSp;
  std::vector<Acts::SpacePointIndex2> middleSp;
  std::vector<Acts::SpacePointIndex2> topSp;

  for (const auto [bottom, middle, top] : grid.binnedGround()) {
    ACTS_VERBOSE("Process middle " << middle);

    bottomSp.clear();
    for (const auto b : bottom) {
      bottomSp.insert(bottomSp.end(), grid.at(b).begin(), grid.at(b).end());
    }
    middleSp.clear();
    middleSp = grid.at(middle);
    topSp.clear();
    for (const auto t : top) {
      topSp.insert(topSp.end(), grid.at(t).begin(), grid.at(t).end());
    }

    // we compute this here since all middle space point candidates belong to
    // the same z-bin
    auto firstMiddleSp = coreSpacePoints.at(middleSp.front());
    auto [minRadiusRangeForMiddle, maxRadiusRangeForMiddle] =
        retrieveRadiusRangeForMiddle(firstMiddleSp,
                                     state.options.rMiddleSpRange);
    ACTS_VERBOSE("Validity range (radius) for the middle space point is ["
                 << minRadiusRangeForMiddle << ", " << maxRadiusRangeForMiddle
                 << "]");

    for (Acts::SpacePointIndex2 middleSpIndex : middleSp) {
      auto spM = coreSpacePoints.at(middleSpIndex);
      const float rM = spM.extra(rColumn);

      // check if spM is outside our radial region of interest
      if (rM < minRadiusRangeForMiddle) {
        continue;
      }
      if (rM > maxRadiusRangeForMiddle) {
        // break because SPs are sorted in r
        break;
      }

      m_seedFinder->createSeeds(
          state, cache,
          Acts::TripletSeedFinder2::ContainerPointers(
              coreSpacePoints, rColumn, varianceRColumn, varianceZColumn),
          bottomSp, middleSpIndex, topSp, seeds);
    }  // loop on middle space points
  }

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePoints.size() << " space points");

  // we have seeds of proxies
  // convert them to seed of external space points
  SimSeedContainer seedContainerForStorage;
  seedContainerForStorage.reserve(seeds.size());
  for (const auto& seed : seeds) {
    auto sps = seed.spacePointIndices();
    seedContainerForStorage.emplace_back(*coreSpacePoints.at(sps[0])
                                              .sourceLinks()[0]
                                              .get<const SimSpacePoint*>(),
                                         *coreSpacePoints.at(sps[1])
                                              .sourceLinks()[0]
                                              .get<const SimSpacePoint*>(),
                                         *coreSpacePoints.at(sps[2])
                                              .sourceLinks()[0]
                                              .get<const SimSpacePoint*>());
    seedContainerForStorage.back().setVertexZ(seed.vertexZ());
    seedContainerForStorage.back().setQuality(seed.quality());
  }

  m_outputSeeds(ctx, std::move(seedContainerForStorage));
  return ProcessCode::SUCCESS;
}

std::pair<float, float> TripletSeedingAlgorithm2::retrieveRadiusRangeForMiddle(
    const Acts::ConstSpacePointProxy2& spM,
    const Acts::Range1D<float>& rMiddleSpRange) const {
  if (m_cfg.useVariableMiddleSPRange) {
    return {rMiddleSpRange.min(), rMiddleSpRange.max()};
  }
  if (m_cfg.rRangeMiddleSP.empty()) {
    return {m_cfg.rMinMiddle, m_cfg.rMaxMiddle};
  }

  /// get zBin position of the middle SP
  auto pVal =
      std::lower_bound(m_cfg.zBinEdges.begin(), m_cfg.zBinEdges.end(), spM.z());
  int zBin = std::distance(m_cfg.zBinEdges.begin(), pVal);
  /// protects against zM at the limit of zBinEdges
  zBin == 0 ? zBin : --zBin;
  return {m_cfg.rRangeMiddleSP[zBin][0], m_cfg.rRangeMiddleSP[zBin][1]};
}

}  // namespace ActsExamples
