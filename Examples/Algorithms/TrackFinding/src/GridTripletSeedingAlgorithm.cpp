// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/GridTripletSeedingAlgorithm.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding2/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding2/BroadTripletSeedFinder.hpp"
#include "Acts/Seeding2/DoubletSeedFinder.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <cmath>
#include <csignal>
#include <cstddef>
#include <stdexcept>

namespace ActsExamples {

GridTripletSeedingAlgorithm::GridTripletSeedingAlgorithm(
    const Config& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("GridTripletSeedingAlgorithm", lvl), m_cfg(cfg) {
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
  // TODO switch to `m_cfg.rMin`
  m_gridConfig.rMin = 0;
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
  m_gridConfig.rBinEdges = {};
  m_gridConfig.zBinEdges = m_cfg.zBinEdges;
  m_gridConfig.bFieldInZ = m_cfg.bFieldInZ;
  m_gridConfig.bottomBinFinder.emplace(m_cfg.numPhiNeighbors,
                                       m_cfg.zBinNeighborsBottom, 0);
  m_gridConfig.topBinFinder.emplace(m_cfg.numPhiNeighbors,
                                    m_cfg.zBinNeighborsTop, 0);
  m_gridConfig.navigation[0ul] = {};
  m_gridConfig.navigation[1ul] = m_cfg.zBinsCustomLooping;
  m_gridConfig.navigation[2ul] = {};

  m_seedFinder = Acts::Experimental::BroadTripletSeedFinder(
      logger().cloneWithSuffix("Finder"));

  Acts::Experimental::BroadTripletSeedFilter::Config filterConfig;
  filterConfig.deltaInvHelixDiameter = m_cfg.deltaInvHelixDiameter;
  filterConfig.deltaRMin = m_cfg.deltaRMin;
  filterConfig.compatSeedWeight = m_cfg.compatSeedWeight;
  filterConfig.impactWeightFactor = m_cfg.impactWeightFactor;
  filterConfig.zOriginWeightFactor = m_cfg.zOriginWeightFactor;
  filterConfig.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  filterConfig.compatSeedLimit = m_cfg.compatSeedLimit;
  filterConfig.seedWeightIncrement = m_cfg.seedWeightIncrement;
  filterConfig.numSeedIncrement = m_cfg.numSeedIncrement;
  filterConfig.seedConfirmation = m_cfg.seedConfirmation;
  filterConfig.centralSeedConfirmationRange =
      m_cfg.centralSeedConfirmationRange;
  filterConfig.forwardSeedConfirmationRange =
      m_cfg.forwardSeedConfirmationRange;
  filterConfig.maxSeedsPerSpMConf = std::numeric_limits<std::size_t>::max();
  filterConfig.maxQualitySeedsPerSpMConf =
      std::numeric_limits<std::size_t>::max();
  filterConfig.useDeltaRinsteadOfTopRadius = m_cfg.useDeltaRinsteadOfTopRadius;

  m_seedFilter = Acts::Experimental::BroadTripletSeedFilter(
      filterConfig, logger().cloneWithSuffix("Filter"));
}

ProcessCode GridTripletSeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const SimSpacePointContainer& spacePoints = m_inputSpacePoints(ctx);

  Acts::Experimental::SpacePointContainer2 coreSpacePoints;
  coreSpacePoints.createExtraColumns(
      Acts::Experimental::SpacePointKnownExtraColumn::R |
      Acts::Experimental::SpacePointKnownExtraColumn::Phi |
      Acts::Experimental::SpacePointKnownExtraColumn::VarianceR |
      Acts::Experimental::SpacePointKnownExtraColumn::VarianceZ);
  coreSpacePoints.reserve(spacePoints.size());
  for (const auto& sp : spacePoints) {
    // check if the space point passes the selection
    if (m_spacePointSelector(sp)) {
      auto newSp = coreSpacePoints.createSpacePoint(
          std::array<Acts::SourceLink, 1>{Acts::SourceLink(&sp)}, sp.x(),
          sp.y(), sp.z());
      newSp.r() = sp.r();
      newSp.phi() = std::atan2(sp.y(), sp.x());
      newSp.varianceR() = sp.varianceR();
      newSp.varianceZ() = sp.varianceZ();
    }
  }

  Acts::Experimental::CylindricalSpacePointGrid2 grid(
      m_gridConfig, logger().cloneWithSuffix("Grid"));

  grid.fill(coreSpacePoints);

  // Compute radius Range
  // we rely on the fact the grid is storing the proxies
  // with a sorting in the radius
  const Acts::Range1D<float> rRange = grid.computeRadiusRange(coreSpacePoints);

  Acts::Experimental::BroadTripletSeedFinder::Options finderOptions;
  finderOptions.bFieldInZ = m_cfg.bFieldInZ;

  Acts::Experimental::DoubletSeedFinder::Config bottomDoubletFinderConfig;
  bottomDoubletFinderConfig.candidateDirection = Acts::Direction::Backward();
  bottomDoubletFinderConfig.deltaRMin = std::isnan(m_cfg.deltaRMaxBottom)
                                            ? m_cfg.deltaRMin
                                            : m_cfg.deltaRMinBottom;
  bottomDoubletFinderConfig.deltaRMax = std::isnan(m_cfg.deltaRMaxBottom)
                                            ? m_cfg.deltaRMax
                                            : m_cfg.deltaRMaxBottom;
  bottomDoubletFinderConfig.deltaZMin = m_cfg.deltaZMin;
  bottomDoubletFinderConfig.deltaZMax = m_cfg.deltaZMax;
  bottomDoubletFinderConfig.impactMax = m_cfg.impactMax;
  bottomDoubletFinderConfig.interactionPointCut = m_cfg.interactionPointCut;
  bottomDoubletFinderConfig.collisionRegionMin = m_cfg.collisionRegionMin;
  bottomDoubletFinderConfig.collisionRegionMax = m_cfg.collisionRegionMax;
  bottomDoubletFinderConfig.cotThetaMax = m_cfg.cotThetaMax;
  bottomDoubletFinderConfig.minPt = m_cfg.minPt;
  bottomDoubletFinderConfig.helixCutTolerance = m_cfg.helixCutTolerance;
  if (m_cfg.useExtraCuts) {
    bottomDoubletFinderConfig.experimentCuts.connect<itkFastTrackingCuts>();
  }
  Acts::Experimental::DoubletSeedFinder bottomDoubletFinder(
      Acts::Experimental::DoubletSeedFinder::DerivedConfig(
          bottomDoubletFinderConfig, m_cfg.bFieldInZ));

  Acts::Experimental::DoubletSeedFinder::Config topDoubletFinderConfig =
      bottomDoubletFinderConfig;
  topDoubletFinderConfig.candidateDirection = Acts::Direction::Forward();
  topDoubletFinderConfig.deltaRMin =
      std::isnan(m_cfg.deltaRMaxTop) ? m_cfg.deltaRMin : m_cfg.deltaRMinTop;
  topDoubletFinderConfig.deltaRMax =
      std::isnan(m_cfg.deltaRMaxTop) ? m_cfg.deltaRMax : m_cfg.deltaRMaxTop;
  Acts::Experimental::DoubletSeedFinder topDoubletFinder(
      Acts::Experimental::DoubletSeedFinder::DerivedConfig(
          topDoubletFinderConfig, m_cfg.bFieldInZ));

  Acts::Experimental::BroadTripletSeedFinder::TripletCuts tripletCuts;
  tripletCuts.minPt = m_cfg.minPt;
  tripletCuts.sigmaScattering = m_cfg.sigmaScattering;
  tripletCuts.radLengthPerSeed = m_cfg.radLengthPerSeed;
  tripletCuts.maxPtScattering = m_cfg.maxPtScattering;
  tripletCuts.impactMax = m_cfg.impactMax;
  tripletCuts.helixCutTolerance = m_cfg.helixCutTolerance;
  tripletCuts.toleranceParam = m_cfg.toleranceParam;
  Acts::Experimental::BroadTripletSeedFinder::DerivedTripletCuts
      derivedTripletCuts(tripletCuts, m_cfg.bFieldInZ);

  // variable middle SP radial region of interest
  Acts::Range1D<float> rMiddleSpRange = {
      std::floor(rRange.min() / 2) * 2 + m_cfg.deltaRMiddleMinSPRange,
      std::floor(rRange.max() / 2) * 2 - m_cfg.deltaRMiddleMaxSPRange};

  // run the seeding
  Acts::Experimental::SeedContainer2 seeds;
  Acts::Experimental::BroadTripletSeedFinder::State state;
  static thread_local Acts::Experimental::BroadTripletSeedFinder::Cache cache;

  std::vector<std::span<const Acts::SpacePointIndex2>> bottomSpGroups;
  std::span<const Acts::SpacePointIndex2> middleSps;
  std::vector<std::span<const Acts::SpacePointIndex2>> topSpGroups;

  for (const auto [bottom, middle, top] : grid.binnedGroup()) {
    ACTS_VERBOSE("Process middle " << middle);

    bottomSpGroups.clear();
    for (const auto b : bottom) {
      bottomSpGroups.push_back(grid.at(b));
    }
    middleSps = grid.at(middle);
    topSpGroups.clear();
    for (const auto t : top) {
      topSpGroups.push_back(grid.at(t));
    }

    // we compute this here since all middle space point candidates belong to
    // the same z-bin
    auto firstMiddleSp = coreSpacePoints.at(middleSps.front());
    auto radiusRangeForMiddle = retrieveRadiusRangeForMiddle(
        Acts::Experimental::ConstSpacePointProxy2(firstMiddleSp),
        rMiddleSpRange);
    ACTS_VERBOSE("Validity range (radius) for the middle space point is ["
                 << radiusRangeForMiddle.first << ", "
                 << radiusRangeForMiddle.second << "]");

    m_seedFinder->createSeedsFromSortedGroups(
        finderOptions, state, cache, bottomDoubletFinder, topDoubletFinder,
        derivedTripletCuts, *m_seedFilter, coreSpacePoints, bottomSpGroups,
        middleSps, topSpGroups, radiusRangeForMiddle, seeds);
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
GridTripletSeedingAlgorithm::retrieveRadiusRangeForMiddle(
    const Acts::Experimental::ConstSpacePointProxy2& spM,
    const Acts::Range1D<float>& rMiddleSpRange) const {
  if (m_cfg.useVariableMiddleSPRange) {
    return {rMiddleSpRange.min(), rMiddleSpRange.max()};
  }
  if (m_cfg.rRangeMiddleSP.empty()) {
    return {m_cfg.rMinMiddle, m_cfg.rMaxMiddle};
  }

  // get zBin position of the middle SP
  auto pVal = std::ranges::lower_bound(m_cfg.zBinEdges, spM.z());
  std::size_t zBin = std::distance(m_cfg.zBinEdges.begin(), pVal);
  // protects against zM at the limit of zBinEdges
  zBin == 0 ? zBin : --zBin;
  return {m_cfg.rRangeMiddleSP[zBin][0], m_cfg.rRangeMiddleSP[zBin][1]};
}

}  // namespace ActsExamples
