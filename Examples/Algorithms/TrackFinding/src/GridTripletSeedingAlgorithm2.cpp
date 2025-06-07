// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/GridTripletSeedingAlgorithm2.hpp"

#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding2/DoubletFinder2.hpp"
#include "Acts/Seeding2/GroupedTripletSeedFilter2.hpp"
#include "Acts/Seeding2/GroupedTripletSeedFinder2.hpp"
#include "Acts/Seeding2/SpacePointContainerPointers2.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <cmath>
#include <csignal>
#include <cstddef>
#include <stdexcept>

namespace ActsExamples {

GridTripletSeedingAlgorithm2::GridTripletSeedingAlgorithm2(
    const Config& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("GridTripletSeedingAlgorithm2", lvl), m_cfg(cfg) {
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputSeeds.initialize(m_cfg.outputSeeds);

  // check that the bins required in the custom bin looping
  // are contained in the bins defined by the total number of edges
  for (std::size_t i : m_cfg.zBinsCustomLooping) {
    if (i >= m_cfg.zBinEdges.size()) {
      throw std::invalid_argument(
          "Inconsistent config zBinsCustomLooping does not contain a subset "
          "of bins defined by zBinEdges");
    }
  }

  if (m_cfg.useExtraCuts) {
    // This function will be applied to select space points during grid filling
    m_spacePointSelector.connect<itkFastTrackingSPselect>();
  }

  m_gridConfig.minPt = m_cfg.minPt;
  m_gridConfig.rMin = m_cfg.rMin;
  m_gridConfig.rMax = m_cfg.rMax;
  m_gridConfig.zMin = m_cfg.zMin;
  m_gridConfig.zMax = m_cfg.zMax;
  m_gridConfig.deltaRMax = m_cfg.deltaRMax;
  m_gridConfig.cotThetaMax = m_cfg.cotThetaMax;
  m_gridConfig.impactMax = m_cfg.impactMax;
  m_gridConfig.phiMin = m_cfg.phiMin;
  m_gridConfig.phiMax = m_cfg.phiMax;
  m_gridConfig.phiBinDeflectionCoverage = m_cfg.phiBinDeflectionCoverage;
  m_gridConfig.maxPhiBins = m_cfg.maxPhiBins;
  m_gridConfig.bFieldInZ = m_cfg.bFieldInZ;
  m_gridConfig.bottomBinFinder.emplace(m_cfg.numPhiNeighbors,
                                       m_cfg.zBinNeighborsBottom, 0);
  m_gridConfig.topBinFinder.emplace(m_cfg.numPhiNeighbors,
                                    m_cfg.zBinNeighborsTop, 0);
  m_gridConfig.navigation[1ul] = m_cfg.zBinsCustomLooping;

  m_seedFinder =
      Acts::GroupedTripletSeedFinder2({}, logger().cloneWithSuffix("Finder"));

  Acts::GroupedTripletSeedFilter2::Config filterConfig;
  filterConfig.deltaInvHelixDiameter = m_cfg.deltaInvHelixDiameter;
  filterConfig.deltaRMin = m_cfg.deltaRMin;
  filterConfig.compatSeedWeight = m_cfg.compatSeedWeight;
  filterConfig.impactWeightFactor = m_cfg.impactWeightFactor;
  filterConfig.zOriginWeightFactor = m_cfg.zOriginWeightFactor;
  filterConfig.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  filterConfig.compatSeedLimit = m_cfg.compatSeedLimit;
  filterConfig.seedWeightIncrement = m_cfg.seedWeightIncrement;
  filterConfig.numSeedIncrement = m_cfg.numSeedIncrement;

  m_seedFilter = Acts::GroupedTripletSeedFilter2(
      filterConfig, logger().cloneWithSuffix("Filter"));
}

ProcessCode GridTripletSeedingAlgorithm2::execute(
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

  Acts::CylindricalSpacePointGrid2 grid(m_gridConfig,
                                        logger().cloneWithSuffix("Grid"));

  grid.fill(coreSpacePoints, phiColumn, rColumn);

  // Compute radius Range
  // we rely on the fact the grid is storing the proxies
  // with a sorting in the radius
  const Acts::Range1D<float> rRange =
      grid.computeRadiusRange(coreSpacePoints, rColumn);

  Acts::GroupedTripletSeedFinder2::Options finderOptions;
  finderOptions.bFieldInZ = m_cfg.bFieldInZ;

  Acts::GroupedTripletSeedFilter2::Options filterOptions;
  filterOptions.seedConfirmation = m_cfg.seedConfirmation;
  filterOptions.seedConfRange = m_cfg.centralSeedConfirmationRange;

  Acts::DoubletFinder2::Cuts doubletCuts;
  doubletCuts.deltaRMin = m_cfg.deltaRMin;
  doubletCuts.deltaRMax = m_cfg.deltaRMax;
  doubletCuts.deltaZMin = m_cfg.deltaZMin;
  doubletCuts.deltaZMax = m_cfg.deltaZMax;
  doubletCuts.impactMax = m_cfg.impactMax;
  doubletCuts.interactionPointCut = m_cfg.interactionPointCut;
  doubletCuts.collisionRegionMin = m_cfg.collisionRegionMin;
  doubletCuts.collisionRegionMax = m_cfg.collisionRegionMax;
  doubletCuts.cotThetaMax = m_cfg.cotThetaMax;
  doubletCuts.minPt = m_cfg.minPt;
  doubletCuts.helixCutTolerance = m_cfg.helixCutTolerance;
  if (m_cfg.useExtraCuts) {
    doubletCuts.experimentCuts.connect<itkFastTrackingCuts>();
  }
  auto derivedDoubletOptions = doubletCuts.derive(m_cfg.bFieldInZ);

  Acts::GroupedTripletSeedFinder2::TripletCuts tripletCuts;
  tripletCuts.minPt = m_cfg.minPt;
  tripletCuts.sigmaScattering = m_cfg.sigmaScattering;
  tripletCuts.radLengthPerSeed = m_cfg.radLengthPerSeed;
  tripletCuts.maxPtScattering = m_cfg.maxPtScattering;
  tripletCuts.impactMax = m_cfg.impactMax;
  tripletCuts.helixCutTolerance = m_cfg.helixCutTolerance;
  tripletCuts.toleranceParam = m_cfg.toleranceParam;
  auto derivedTripletOptions = tripletCuts.derive(finderOptions);

  /// variable middle SP radial region of interest
  Acts::Range1D<float> rMiddleSpRange = {
      std::floor(rRange.min() / 2) * 2 + m_cfg.deltaRMiddleMinSPRange,
      std::floor(rRange.max() / 2) * 2 - m_cfg.deltaRMiddleMaxSPRange};

  // run the seeding
  Acts::SeedContainer2 seeds;
  Acts::GroupedTripletSeedFinder2::State state;
  static thread_local Acts::GroupedTripletSeedFinder2::Cache cache;

  std::vector<Acts::SpacePointIndex2> bottomSps;
  std::vector<Acts::SpacePointIndex2> middleSps;
  std::vector<Acts::SpacePointIndex2> topSps;

  for (const auto [bottom, middle, top] : grid.binnedGround()) {
    ACTS_VERBOSE("Process middle " << middle);

    bottomSps.clear();
    for (const auto b : bottom) {
      bottomSps.insert(bottomSps.end(), grid.at(b).begin(), grid.at(b).end());
    }
    middleSps.clear();
    middleSps = grid.at(middle);
    topSps.clear();
    for (const auto t : top) {
      topSps.insert(topSps.end(), grid.at(t).begin(), grid.at(t).end());
    }

    std::ranges::sort(bottomSps, {}, [&](Acts::SpacePointIndex2 spIndex) {
      return coreSpacePoints.at(spIndex).extra(rColumn);
    });
    std::ranges::sort(middleSps, {}, [&](Acts::SpacePointIndex2 spIndex) {
      return coreSpacePoints.at(spIndex).extra(rColumn);
    });
    std::ranges::sort(topSps, {}, [&](Acts::SpacePointIndex2 spIndex) {
      return coreSpacePoints.at(spIndex).extra(rColumn);
    });

    // we compute this here since all middle space point candidates belong to
    // the same z-bin
    auto firstMiddleSp = coreSpacePoints.at(middleSps.front());
    auto [minRadiusRangeForMiddle, maxRadiusRangeForMiddle] =
        retrieveRadiusRangeForMiddle(firstMiddleSp, rMiddleSpRange);
    ACTS_VERBOSE("Validity range (radius) for the middle space point is ["
                 << minRadiusRangeForMiddle << ", " << maxRadiusRangeForMiddle
                 << "]");

    for (Acts::SpacePointIndex2 middleSp : middleSps) {
      auto spM = coreSpacePoints.at(middleSp);
      const float zM = spM.z();
      const float rM = spM.extra(rColumn);

      // check if spM is outside our radial region of interest
      if (rM < minRadiusRangeForMiddle) {
        continue;
      }
      if (rM > maxRadiusRangeForMiddle) {
        // break because SPs are sorted in r
        break;
      }

      // apply cut on the number of top SP if seedConfirmation is true
      if (m_cfg.seedConfirmation) {
        // check if middle SP is in the central or forward region
        filterOptions.seedConfRange =
            (zM > m_cfg.centralSeedConfirmationRange.zMaxSeedConf ||
             zM < m_cfg.centralSeedConfirmationRange.zMinSeedConf)
                ? m_cfg.forwardSeedConfirmationRange
                : m_cfg.centralSeedConfirmationRange;
        // set the minimum number of top SP depending on whether the middle SP
        // is in the central or forward region
        filterOptions.nTopSeedConf =
            rM > filterOptions.seedConfRange.rMaxSeedConf
                ? filterOptions.seedConfRange.nTopForLargeR
                : filterOptions.seedConfRange.nTopForSmallR;
      }

      m_seedFinder->createSeedsFromGroup(
          finderOptions, state, cache, derivedDoubletOptions,
          derivedDoubletOptions, derivedTripletOptions, *m_seedFilter,
          Acts::SpacePointContainerPointers2(coreSpacePoints, rColumn,
                                             varianceRColumn, varianceZColumn),
          bottomSps, middleSp, topSps, seeds);
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

std::pair<float, float>
GridTripletSeedingAlgorithm2::retrieveRadiusRangeForMiddle(
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
