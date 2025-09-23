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
    const Acts::Experimental::ConstSpacePointProxy2 & /*middle*/,
    const Acts::Experimental::ConstSpacePointProxy2 &other, float cotTheta,
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
    : IAlgorithm("OrthogonalTripletSeedingAlgorithm", lvl), m_cfg(cfg) {
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

  m_seedFinder =
      Acts::Experimental::TripletSeeder(logger().cloneWithSuffix("Finder"));
}

ProcessCode OrthogonalTripletSeedingAlgorithm::execute(
    const AlgorithmContext &ctx) const {
  const SimSpacePointContainer &spacePoints = m_inputSpacePoints(ctx);

  Acts::Experimental::SpacePointContainer2 coreSpacePoints(
      Acts::Experimental::SpacePointColumns::SourceLinks |
      Acts::Experimental::SpacePointColumns::XY |
      Acts::Experimental::SpacePointColumns::ZR |
      Acts::Experimental::SpacePointColumns::Phi |
      Acts::Experimental::SpacePointColumns::VarianceZ |
      Acts::Experimental::SpacePointColumns::VarianceR);
  coreSpacePoints.reserve(spacePoints.size());

  std::vector<typename tree_t::pair_t> kdTreePoints;
  kdTreePoints.reserve(spacePoints.size());

  Acts::Extent rRangeSPExtent;

  for (std::size_t i = 0; i < spacePoints.size(); ++i) {
    const auto &sp = spacePoints[i];

    // check if the space point passes the selection
    if (m_spacePointSelector.connected() && !m_spacePointSelector(sp)) {
      continue;
    }

    Acts::Experimental::SpacePointIndex2 newSpIndex = coreSpacePoints.size();
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

    /*
     * For every input point, we create a coordinate-pointer pair, which we then
     * linearly pass to the k-d tree constructor. That constructor will take
     * care of sorting the pairs and splitting the space.
     */
    typename tree_t::coordinate_t kdTreePoint;

    kdTreePoint[DimPhi] = newSp.phi();
    kdTreePoint[DimR] = newSp.zr()[1];
    kdTreePoint[DimZ] = newSp.zr()[0];

    kdTreePoints.emplace_back(kdTreePoint, newSpIndex);

    rRangeSPExtent.extend({newSp.xy()[0], newSp.xy()[1], newSp.zr()[0]});
  }

  ACTS_VERBOSE("Created k-d tree populated with " << kdTreePoints.size()
                                                  << " space points");

  tree_t kdTree(std::move(kdTreePoints));

  Acts::Experimental::DoubletSeedFinder::Config bottomDoubletFinderConfig;
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
  auto bottomDoubletFinder = Acts::Experimental::DoubletSeedFinder::create(
      Acts::Experimental::DoubletSeedFinder::DerivedConfig(
          bottomDoubletFinderConfig, m_cfg.bFieldInZ));

  Acts::Experimental::DoubletSeedFinder::Config topDoubletFinderConfig =
      bottomDoubletFinderConfig;
  topDoubletFinderConfig.candidateDirection = Acts::Direction::Forward();
  topDoubletFinderConfig.deltaRMin =
      std::isnan(m_cfg.deltaRMaxTop) ? m_cfg.deltaRMin : m_cfg.deltaRMinTop;
  topDoubletFinderConfig.deltaRMax =
      std::isnan(m_cfg.deltaRMaxTop) ? m_cfg.deltaRMax : m_cfg.deltaRMaxTop;
  auto topDoubletFinder = Acts::Experimental::DoubletSeedFinder::create(
      Acts::Experimental::DoubletSeedFinder::DerivedConfig(
          topDoubletFinderConfig, m_cfg.bFieldInZ));

  Acts::Experimental::TripletSeedFinder::Config tripletFinderConfig;
  tripletFinderConfig.useStripInfo = false;
  tripletFinderConfig.sortedByCotTheta = true;
  tripletFinderConfig.minPt = m_cfg.minPt;
  tripletFinderConfig.sigmaScattering = m_cfg.sigmaScattering;
  tripletFinderConfig.radLengthPerSeed = m_cfg.radLengthPerSeed;
  tripletFinderConfig.impactMax = m_cfg.impactMax;
  tripletFinderConfig.helixCutTolerance = m_cfg.helixCutTolerance;
  tripletFinderConfig.toleranceParam = m_cfg.toleranceParam;
  auto tripletFinder = Acts::Experimental::TripletSeedFinder::create(
      Acts::Experimental::TripletSeedFinder::DerivedConfig(tripletFinderConfig,
                                                           m_cfg.bFieldInZ));

  // variable middle SP radial region of interest
  const Acts::Range1D<float> rMiddleSPRange(
      std::floor(rRangeSPExtent.min(Acts::AxisDirection::AxisR) / 2) * 2 +
          m_cfg.deltaRMiddleMinSPRange,
      std::floor(rRangeSPExtent.max(Acts::AxisDirection::AxisR) / 2) * 2 -
          m_cfg.deltaRMiddleMaxSPRange);

  // run the seeding
  Acts::Experimental::SeedContainer2 seeds;
  Acts::Experimental::BroadTripletSeedFilter::State filterState;
  Acts::Experimental::BroadTripletSeedFilter::Cache filterCache;
  Acts::Experimental::BroadTripletSeedFilter seedFilter(
      m_filterConfig, filterState, filterCache, *m_filterLogger);
  static thread_local Acts::Experimental::TripletSeeder::Cache cache;
  static thread_local SpacePointCandidates candidates;

  /*
   * Run the seeding algorithm by iterating over all the points in the tree
   * and seeing what happens if we take them to be our middle spacepoint.
   */
  for (const typename tree_t::pair_t &middle : kdTree) {
    ACTS_VERBOSE("Process middle " << middle.second);

    const auto spM = coreSpacePoints.at(middle.second).asConst();

    /*
     * Cut: Ensure that the middle spacepoint lies within a valid r-region for
     * middle points.
     */
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
    if (const float zM = spM.zr()[0]; zM < m_cfg.zOutermostLayers.first ||
                                      zM > m_cfg.zOutermostLayers.second) {
      continue;
    }
    if (const float phiM = spM.phi();
        phiM > m_cfg.phiMax || phiM < m_cfg.phiMin) {
      continue;
    }

    candidates.clear();
    findSpacePointCandidates(spM, kdTree, candidates);

    Acts::Experimental::SpacePointContainer2::ConstSubset bottomSps =
        coreSpacePoints.subset(candidates.bottom_lh_v).asConst();
    Acts::Experimental::SpacePointContainer2::ConstSubset topSps =
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

void OrthogonalTripletSeedingAlgorithm::findSpacePointCandidates(
    const Acts::Experimental::ConstSpacePointProxy2 &spM, const tree_t &tree,
    SpacePointCandidates &candidates) const {
  using range_t = typename tree_t::range_t;

  /*
   * Calculate the search ranges for bottom and top candidates for this middle
   * space point.
   */
  range_t bottom_r = validTupleOrthoRangeHL(spM);
  range_t top_r = validTupleOrthoRangeLH(spM);

  /*
   * Calculate the value of cot(θ) for this middle spacepoint.
   */
  float cotTheta =
      std::max(std::abs(spM.zr()[0] / spM.zr()[1]), m_cfg.cotThetaMax);

  /*
   * Calculate the maximum Δr, given that we have already constrained our
   * search space.
   */
  float deltaRMaxTop = top_r[DimR].max() - spM.zr()[1];
  float deltaRMaxBottom = spM.zr()[1] - bottom_r[DimR].min();

  /*
   * Create the search range for the bottom spacepoint assuming a
   * monotonically increasing z track, by calculating the minimum z value from
   * the cot(θ), and by setting the maximum to the z position of the middle
   * spacepoint - if the z position is higher than the middle point, then it
   * would be a decreasing z track!
   */
  range_t bottom_lh_r = bottom_r;
  bottom_lh_r[DimZ].shrink(spM.zr()[0] - cotTheta * deltaRMaxBottom,
                           spM.zr()[0]);

  /*
   * Calculate the search ranges for the other four sets of points in a
   * similar fashion.
   */
  range_t top_lh_r = top_r;
  top_lh_r[DimZ].shrink(spM.zr()[0], spM.zr()[0] + cotTheta * deltaRMaxTop);

  range_t bottom_hl_r = bottom_r;
  bottom_hl_r[DimZ].shrink(spM.zr()[0],
                           spM.zr()[0] + cotTheta * deltaRMaxBottom);
  range_t top_hl_r = top_r;
  top_hl_r[DimZ].shrink(spM.zr()[0] - cotTheta * deltaRMaxTop, spM.zr()[0]);

  /*
   * Now, we will actually search for the spaces. Remembering that we combine
   * bottom and top candidates for increasing and decreasing tracks
   * separately, we will first check whether both the search ranges for
   * increasing tracks are not degenerate - if they are, we will never find
   * any seeds and we do not need to bother doing the search.
   */
  if (!bottom_lh_r.degenerate() && !top_lh_r.degenerate()) {
    /*
     * Search the trees for points that lie in the given search range.
     */
    tree.rangeSearchMapDiscard(
        top_lh_r, [&candidates](const typename tree_t::coordinate_t &,
                                const typename tree_t::value_t &top) {
          candidates.top_lh_v.push_back(top);
        });
  }

  /*
   * Perform the same search for candidate bottom spacepoints, but for
   * monotonically decreasing z tracks.
   */
  if (!bottom_hl_r.degenerate() && !top_hl_r.degenerate()) {
    tree.rangeSearchMapDiscard(
        top_hl_r, [&candidates](const typename tree_t::coordinate_t &,
                                const typename tree_t::value_t &top) {
          candidates.top_hl_v.push_back(top);
        });
  }

  // apply cut on the number of top SP if seedConfirmation is true
  bool search_bot_hl = true;
  bool search_bot_lh = true;
  if (m_cfg.seedConfirmation) {
    // check if middle SP is in the central or forward region
    Acts::SeedConfirmationRangeConfig seedConfRange =
        (spM.zr()[0] > m_cfg.centralSeedConfirmationRange.zMaxSeedConf ||
         spM.zr()[0] < m_cfg.centralSeedConfirmationRange.zMinSeedConf)
            ? m_cfg.forwardSeedConfirmationRange
            : m_cfg.centralSeedConfirmationRange;
    // set the minimum number of top SP depending on whether the middle SP is
    // in the central or forward region
    std::size_t nTopSeedConf = spM.zr()[1] > seedConfRange.rMaxSeedConf
                                   ? seedConfRange.nTopForLargeR
                                   : seedConfRange.nTopForSmallR;
    // continue if number of top SPs is smaller than minimum
    if (candidates.top_lh_v.size() < nTopSeedConf) {
      search_bot_lh = false;
    }
    if (candidates.top_hl_v.size() < nTopSeedConf) {
      search_bot_hl = false;
    }
  }

  /*
   * Next, we perform a search for bottom candidates in increasing z tracks,
   * which only makes sense if we found any bottom candidates.
   */
  if (!candidates.top_lh_v.empty() && search_bot_lh) {
    tree.rangeSearchMapDiscard(
        bottom_lh_r, [&candidates](const typename tree_t::coordinate_t &,
                                   const typename tree_t::value_t &bottom) {
          candidates.bottom_lh_v.push_back(bottom);
        });
  }

  /*
   * And repeat for the top spacepoints for decreasing z tracks!
   */
  if (!candidates.top_hl_v.empty() && search_bot_hl) {
    tree.rangeSearchMapDiscard(
        bottom_hl_r, [&candidates](const typename tree_t::coordinate_t &,
                                   const typename tree_t::value_t &bottom) {
          candidates.bottom_hl_v.push_back(bottom);
        });
  }
}

auto OrthogonalTripletSeedingAlgorithm::validTupleOrthoRangeLH(
    const Acts::Experimental::ConstSpacePointProxy2 &low) const ->
    typename tree_t::range_t {
  float colMin = m_cfg.collisionRegionMin;
  float colMax = m_cfg.collisionRegionMax;
  float pL = low.phi();
  float rL = low.zr()[1];
  float zL = low.zr()[0];

  typename tree_t::range_t res;

  /*
   * Cut: Ensure that we search only in φ_min ≤ φ ≤ φ_max, as defined by the
   * seeding configuration.
   */
  res[DimPhi].shrinkMin(m_cfg.phiMin);
  res[DimPhi].shrinkMax(m_cfg.phiMax);

  /*
   * Cut: Ensure that we search only in r ≤ r_max, as defined by the seeding
   * configuration.
   */
  res[DimR].shrinkMax(m_cfg.rMax);

  /*
   * Cut: Ensure that we search only in z_min ≤ z ≤ z_max, as defined by the
   * seeding configuration.
   */
  res[DimZ].shrinkMin(m_cfg.zMin);
  res[DimZ].shrinkMax(m_cfg.zMax);

  /*
   * Cut: Ensure that we search only in Δr_min ≤ r - r_L ≤ Δr_max, as defined
   * by the seeding configuration and the given lower spacepoint.
   */
  res[DimR].shrinkMin(rL + m_cfg.deltaRMinTop);
  res[DimR].shrinkMax(rL + m_cfg.deltaRMaxTop);

  /*
   * Cut: Now that we have constrained r, we can use that new r range to
   * further constrain z.
   */
  float zMax = (res[DimR].max() / rL) * (zL - colMin) + colMin;
  float zMin = colMax - (res[DimR].max() / rL) * (colMax - zL);

  /*
   * This cut only works if z_low is outside the collision region for z.
   */
  if (zL > colMin) {
    res[DimZ].shrinkMax(zMax);
  } else if (zL < colMax) {
    res[DimZ].shrinkMin(zMin);
  }

  /*
   * Cut: Shrink the z-range using the maximum cot(θ).
   */
  res[DimZ].shrinkMin(zL - m_cfg.cotThetaMax * (res[DimR].max() - rL));
  res[DimZ].shrinkMax(zL + m_cfg.cotThetaMax * (res[DimR].max() - rL));

  /*
   * Cut: Shrink the φ range, such that Δφ_min ≤ φ - φ_L ≤ Δφ_max
   */
  res[DimPhi].shrinkMin(pL - m_cfg.deltaPhiMax);
  res[DimPhi].shrinkMax(pL + m_cfg.deltaPhiMax);

  // Cut: Ensure that z-distance between SPs is within max and min values.
  res[DimZ].shrinkMin(zL - m_cfg.deltaZMax);
  res[DimZ].shrinkMax(zL + m_cfg.deltaZMax);

  return res;
}

auto OrthogonalTripletSeedingAlgorithm::validTupleOrthoRangeHL(
    const Acts::Experimental::ConstSpacePointProxy2 &high) const ->
    typename tree_t::range_t {
  float pM = high.phi();
  float rM = high.zr()[1];
  float zM = high.zr()[0];

  typename tree_t::range_t res;

  /*
   * Cut: Ensure that we search only in φ_min ≤ φ ≤ φ_max, as defined by the
   * seeding configuration.
   */
  res[DimPhi].shrinkMin(m_cfg.phiMin);
  res[DimPhi].shrinkMax(m_cfg.phiMax);

  /*
   * Cut: Ensure that we search only in r ≤ r_max, as defined by the seeding
   * configuration.
   */
  res[DimR].shrinkMax(m_cfg.rMax);

  /*
   * Cut: Ensure that we search only in z_min ≤ z ≤ z_max, as defined by the
   * seeding configuration.
   */
  res[DimZ].shrinkMin(m_cfg.zMin);
  res[DimZ].shrinkMax(m_cfg.zMax);

  /*
   * Cut: Ensure that we search only in Δr_min ≤ r_H - r ≤ Δr_max, as defined
   * by the seeding configuration and the given higher spacepoint.
   */
  res[DimR].shrinkMin(rM - m_cfg.deltaRMaxBottom);
  res[DimR].shrinkMax(rM - m_cfg.deltaRMinBottom);

  /*
   * Cut: Now that we have constrained r, we can use that new r range to
   * further constrain z.
   */
  float fracR = res[DimR].min() / rM;

  float zMin =
      (zM - m_cfg.collisionRegionMin) * fracR + m_cfg.collisionRegionMin;
  float zMax =
      (zM - m_cfg.collisionRegionMax) * fracR + m_cfg.collisionRegionMax;

  res[DimZ].shrinkMin(std::min(zMin, zM));
  res[DimZ].shrinkMax(std::max(zMax, zM));

  /*
   * Cut: Shrink the φ range, such that Δφ_min ≤ φ - φ_H ≤ Δφ_max
   */
  res[DimPhi].shrinkMin(pM - m_cfg.deltaPhiMax);
  res[DimPhi].shrinkMax(pM + m_cfg.deltaPhiMax);

  // Cut: Ensure that z-distance between SPs is within max and min values.
  res[DimZ].shrinkMin(zM - m_cfg.deltaZMax);
  res[DimZ].shrinkMax(zM + m_cfg.deltaZMax);

  return res;
}

}  // namespace ActsExamples
