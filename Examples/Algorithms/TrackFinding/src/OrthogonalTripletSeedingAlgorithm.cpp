// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/OrthogonalTripletSeedingAlgorithm.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding2/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding2/CylindricalSpacePointKDTree.hpp"
#include "Acts/Seeding2/DoubletSeedFinder.hpp"
#include "Acts/Seeding2/TripletSeedFinder.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <cmath>
#include <csignal>
#include <cstddef>

namespace ActsExamples {

namespace {

static inline bool itkFastTrackingCuts(
    const Acts::ConstSpacePointProxy2 & /*middle*/,
    const Acts::ConstSpacePointProxy2 &other, float cotTheta,
    bool isBottomCandidate) {
  static float rMin = 45;
  static float cotThetaMax = 1.5;

  if (isBottomCandidate && other.zr()[1] < rMin &&
      (cotTheta > cotThetaMax || cotTheta < -cotThetaMax)) {
    return false;
  }
  return true;
}

static inline bool itkFastTrackingSPselect(const SimSpacePoint &sp) {
  // At small r we remove points beyond |z| > 200.
  float r = sp.r();
  float zabs = std::abs(sp.z());
  if (zabs > 200. && r < 45.) {
    return false;
  }

  // Remove space points beyond eta=4 if their z is larger than the max seed
  // z0 (150.)
  float cotTheta = 27.2899;  // corresponds to eta=4
  if ((zabs - 150.) > cotTheta * r) {
    return false;
  }
  return true;
}

}  // namespace

OrthogonalTripletSeedingAlgorithm::OrthogonalTripletSeedingAlgorithm(
    const Config &cfg, Acts::Logging::Level lvl)
    : IAlgorithm(
          "OrthogonalTripletSeedingAlgorithm",
          Acts::getDefaultLogger("OrthogonalTripletSeedingAlgorithm", lvl)),
      m_cfg(cfg) {
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputSeeds.initialize(m_cfg.outputSeeds);

  if (m_cfg.useExtraCuts) {
    // This function will be applied to select space points during grid filling
    m_spacePointSelector.connect<itkFastTrackingSPselect>();
  }

  m_filterConfig.deltaInvHelixDiameter = m_cfg.deltaInvHelixDiameter;
  m_filterConfig.deltaRMin = m_cfg.deltaRMin;
  m_filterConfig.compatSeedWeight = m_cfg.compatSeedWeight;
  m_filterConfig.impactWeightFactor = m_cfg.impactWeightFactor;
  m_filterConfig.zOriginWeightFactor = m_cfg.zOriginWeightFactor;
  m_filterConfig.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  m_filterConfig.compatSeedLimit = m_cfg.compatSeedLimit;
  m_filterConfig.seedWeightIncrement = m_cfg.seedWeightIncrement;
  m_filterConfig.numSeedIncrement = m_cfg.numSeedIncrement;
  m_filterConfig.seedConfirmation = m_cfg.seedConfirmation;
  m_filterConfig.centralSeedConfirmationRange =
      m_cfg.centralSeedConfirmationRange;
  m_filterConfig.forwardSeedConfirmationRange =
      m_cfg.forwardSeedConfirmationRange;
  m_filterConfig.maxSeedsPerSpMConf = m_cfg.maxSeedsPerSpMConf;
  m_filterConfig.maxQualitySeedsPerSpMConf = m_cfg.maxQualitySeedsPerSpMConf;
  m_filterConfig.useDeltaRinsteadOfTopRadius =
      m_cfg.useDeltaRinsteadOfTopRadius;

  m_filterLogger = logger().cloneWithSuffix("Filter");

  m_seedFinder = Acts::TripletSeeder(logger().cloneWithSuffix("Finder"));
}

ProcessCode OrthogonalTripletSeedingAlgorithm::execute(
    const AlgorithmContext &ctx) const {
  const SimSpacePointContainer &spacePoints = m_inputSpacePoints(ctx);

  Acts::SpacePointContainer2 coreSpacePoints(
      Acts::SpacePointColumns::SourceLinks | Acts::SpacePointColumns::PackedXY |
      Acts::SpacePointColumns::PackedZR | Acts::SpacePointColumns::Phi |
      Acts::SpacePointColumns::VarianceZ | Acts::SpacePointColumns::VarianceR);
  coreSpacePoints.reserve(spacePoints.size());

  Acts::Experimental::CylindricalSpacePointKDTreeBuilder kdTreeBuilder;
  kdTreeBuilder.reserve(spacePoints.size());

  Acts::Extent rRangeSPExtent;

  for (std::size_t i = 0; i < spacePoints.size(); ++i) {
    const auto &sp = spacePoints[i];

    // check if the space point passes the selection
    if (m_spacePointSelector.connected() && !m_spacePointSelector(sp)) {
      continue;
    }

    Acts::SpacePointIndex2 newSpIndex = coreSpacePoints.size();
    auto newSp = coreSpacePoints.createSpacePoint();
    newSp.assignSourceLinks(
        std::array<Acts::SourceLink, 1>{Acts::SourceLink(&sp)});
    newSp.xy() = std::array<float, 2>{static_cast<float>(sp.x()),
                                      static_cast<float>(sp.y())};
    newSp.zr() = std::array<float, 2>{static_cast<float>(sp.z()),
                                      static_cast<float>(sp.r())};
    newSp.phi() = static_cast<float>(std::atan2(sp.y(), sp.x()));
    newSp.varianceZ() = static_cast<float>(sp.varianceZ());
    newSp.varianceR() = static_cast<float>(sp.varianceR());

    kdTreeBuilder.insert(newSpIndex, newSp.phi(), newSp.zr()[1], newSp.zr()[0]);

    rRangeSPExtent.extend({newSp.xy()[0], newSp.xy()[1], newSp.zr()[0]});
  }

  ACTS_VERBOSE("Created k-d tree populated with " << kdTreeBuilder.size()
                                                  << " space points");

  Acts::Experimental::CylindricalSpacePointKDTree kdTree =
      kdTreeBuilder.build();

  Acts::Experimental::CylindricalSpacePointKDTree::Options lhOptions;
  lhOptions.rMax = m_cfg.rMax;
  lhOptions.zMin = m_cfg.zMin;
  lhOptions.zMax = m_cfg.zMax;
  lhOptions.phiMin = m_cfg.phiMin;
  lhOptions.phiMax = m_cfg.phiMax;
  lhOptions.deltaRMin = std::isnan(m_cfg.deltaRMinBottom)
                            ? m_cfg.deltaRMin
                            : m_cfg.deltaRMinBottom;
  lhOptions.deltaRMax = std::isnan(m_cfg.deltaRMaxBottom)
                            ? m_cfg.deltaRMax
                            : m_cfg.deltaRMaxBottom;
  lhOptions.collisionRegionMin = m_cfg.collisionRegionMin;
  lhOptions.collisionRegionMax = m_cfg.collisionRegionMax;
  lhOptions.cotThetaMax = m_cfg.cotThetaMax;
  lhOptions.deltaPhiMax = m_cfg.deltaPhiMax;
  Acts::Experimental::CylindricalSpacePointKDTree::Options hlOptions =
      lhOptions;
  hlOptions.deltaRMin =
      std::isnan(m_cfg.deltaRMinTop) ? m_cfg.deltaRMin : m_cfg.deltaRMinTop;
  hlOptions.deltaRMax =
      std::isnan(m_cfg.deltaRMaxTop) ? m_cfg.deltaRMax : m_cfg.deltaRMaxTop;

  Acts::DoubletSeedFinder::Config bottomDoubletFinderConfig;
  bottomDoubletFinderConfig.spacePointsSortedByRadius = false;
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
  auto bottomDoubletFinder =
      Acts::DoubletSeedFinder::create(Acts::DoubletSeedFinder::DerivedConfig(
          bottomDoubletFinderConfig, m_cfg.bFieldInZ));

  Acts::DoubletSeedFinder::Config topDoubletFinderConfig =
      bottomDoubletFinderConfig;
  topDoubletFinderConfig.candidateDirection = Acts::Direction::Forward();
  topDoubletFinderConfig.deltaRMin =
      std::isnan(m_cfg.deltaRMaxTop) ? m_cfg.deltaRMin : m_cfg.deltaRMinTop;
  topDoubletFinderConfig.deltaRMax =
      std::isnan(m_cfg.deltaRMaxTop) ? m_cfg.deltaRMax : m_cfg.deltaRMaxTop;
  auto topDoubletFinder =
      Acts::DoubletSeedFinder::create(Acts::DoubletSeedFinder::DerivedConfig(
          topDoubletFinderConfig, m_cfg.bFieldInZ));

  Acts::TripletSeedFinder::Config tripletFinderConfig;
  tripletFinderConfig.useStripInfo = false;
  tripletFinderConfig.sortedByCotTheta = true;
  tripletFinderConfig.minPt = m_cfg.minPt;
  tripletFinderConfig.sigmaScattering = m_cfg.sigmaScattering;
  tripletFinderConfig.radLengthPerSeed = m_cfg.radLengthPerSeed;
  tripletFinderConfig.impactMax = m_cfg.impactMax;
  tripletFinderConfig.helixCutTolerance = m_cfg.helixCutTolerance;
  tripletFinderConfig.toleranceParam = m_cfg.toleranceParam;
  auto tripletFinder =
      Acts::TripletSeedFinder::create(Acts::TripletSeedFinder::DerivedConfig(
          tripletFinderConfig, m_cfg.bFieldInZ));

  // variable middle SP radial region of interest
  const Acts::Range1D<float> rMiddleSPRange(
      std::floor(rRangeSPExtent.min(Acts::AxisDirection::AxisR) / 2) * 2 +
          m_cfg.deltaRMiddleMinSPRange,
      std::floor(rRangeSPExtent.max(Acts::AxisDirection::AxisR) / 2) * 2 -
          m_cfg.deltaRMiddleMaxSPRange);

  // run the seeding
  Acts::SeedContainer2 seeds;
  Acts::BroadTripletSeedFilter::State filterState;
  Acts::BroadTripletSeedFilter::Cache filterCache;
  Acts::BroadTripletSeedFilter seedFilter(m_filterConfig, filterState,
                                          filterCache, *m_filterLogger);

  static thread_local Acts::TripletSeeder::Cache cache;
  static thread_local Acts::Experimental::CylindricalSpacePointKDTree::
      Candidates candidates;

  // Run the seeding algorithm by iterating over all the points in the tree
  // and seeing what happens if we take them to be our middle space point.
  for (const auto &middle : kdTree) {
    ACTS_VERBOSE("Process middle " << middle.second);

    const auto spM = coreSpacePoints.at(middle.second).asConst();

    // Cut: Ensure that the middle space point lies within a valid r-region for
    // middle points.
    const float rM = spM.zr()[1];
    if (m_cfg.useVariableMiddleSPRange) {
      if (rM < rMiddleSPRange.min() || rM > rMiddleSPRange.max()) {
        continue;
      }
    } else {
      if (rM > m_cfg.rMaxMiddle || rM < m_cfg.rMinMiddle) {
        continue;
      }
    }

    // remove all middle SPs outside phi and z region of interest
    const float zM = spM.zr()[0];
    if (zM < m_cfg.zOutermostLayers.first ||
        zM > m_cfg.zOutermostLayers.second) {
      continue;
    }
    if (const float phiM = spM.phi();
        phiM > m_cfg.phiMax || phiM < m_cfg.phiMin) {
      continue;
    }

    std::size_t nTopSeedConf = 0;
    if (m_cfg.seedConfirmation) {
      // check if middle SP is in the central or forward region
      Acts::SeedConfirmationRangeConfig seedConfRange =
          (zM > m_cfg.centralSeedConfirmationRange.zMaxSeedConf ||
           zM < m_cfg.centralSeedConfirmationRange.zMinSeedConf)
              ? m_cfg.forwardSeedConfirmationRange
              : m_cfg.centralSeedConfirmationRange;
      // set the minimum number of top SP depending on whether the middle SP is
      // in the central or forward region
      nTopSeedConf = rM > seedConfRange.rMaxSeedConf
                         ? seedConfRange.nTopForLargeR
                         : seedConfRange.nTopForSmallR;
    }

    candidates.clear();
    kdTree.validTuples(lhOptions, hlOptions, spM, nTopSeedConf, candidates);

    Acts::SpacePointContainer2::ConstSubset bottomSps =
        coreSpacePoints.subset(candidates.bottom_lh_v).asConst();
    Acts::SpacePointContainer2::ConstSubset topSps =
        coreSpacePoints.subset(candidates.top_lh_v).asConst();
    m_seedFinder->createSeedsFromGroup(
        cache, *bottomDoubletFinder, *topDoubletFinder, *tripletFinder,
        seedFilter, coreSpacePoints, bottomSps, spM, topSps, seeds);

    bottomSps = coreSpacePoints.subset(candidates.bottom_hl_v).asConst();
    topSps = coreSpacePoints.subset(candidates.top_hl_v).asConst();
    m_seedFinder->createSeedsFromGroup(
        cache, *bottomDoubletFinder, *topDoubletFinder, *tripletFinder,
        seedFilter, coreSpacePoints, bottomSps, spM, topSps, seeds);
  }

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePoints.size() << " space points");

  // we have seeds of proxies
  // convert them to seed of external space points
  SimSeedContainer seedContainerForStorage;
  seedContainerForStorage.reserve(seeds.size());
  for (const auto &seed : seeds) {
    auto sps = seed.spacePointIndices();
    seedContainerForStorage.emplace_back(*coreSpacePoints.at(sps[0])
                                              .sourceLinks()[0]
                                              .get<const SimSpacePoint *>(),
                                         *coreSpacePoints.at(sps[1])
                                              .sourceLinks()[0]
                                              .get<const SimSpacePoint *>(),
                                         *coreSpacePoints.at(sps[2])
                                              .sourceLinks()[0]
                                              .get<const SimSpacePoint *>());
    seedContainerForStorage.back().setVertexZ(seed.vertexZ());
    seedContainerForStorage.back().setQuality(seed.quality());
  }

  m_outputSeeds(ctx, std::move(seedContainerForStorage));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
