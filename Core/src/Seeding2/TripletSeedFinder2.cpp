// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/TripletSeedFinder2.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding2/TripletSeedFilter2.hpp"

#include <numeric>

#include <Eigen/src/Core/Matrix.h>

namespace Acts {

using namespace UnitLiterals;

TripletSeedFinder2::DerivedConfig TripletSeedFinder2::Config::derive() const {
  DerivedConfig result;

  static_cast<Config&>(result) = *this;

  // similar to `theta0Highland` in `Core/src/Material/Interactions.cpp`
  {
    const float xOverX0 = result.radLengthPerSeed;
    const float q2OverBeta2 = 1;  // q^2=1, beta^2~1
    // RPP2018 eq. 33.15 (treats beta and q² consistently)
    const float t = std::sqrt(xOverX0 * q2OverBeta2);
    // log((x/X0) * (q²/beta²)) = log((sqrt(x/X0) * (q/beta))²)
    //                          = 2 * log(sqrt(x/X0) * (q/beta))
    result.highland = 13.6_MeV * t * (1.0f + 0.038f * 2 * std::log(t));
  }

  const float maxScatteringAngle = result.highland / result.minPt;
  result.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;

  return result;
}

TripletSeedFinder2::DerivedOptions TripletSeedFinder2::Options::derive(
    const DerivedConfig& config) const {
  DerivedOptions result;

  static_cast<Options&>(result) = *this;

  // bFieldInZ is in (pT/radius) natively, no need for conversion
  result.pTPerHelixRadius = result.bFieldInZ;
  result.minHelixDiameter2 =
      std::pow(config.minPt * 2 / result.pTPerHelixRadius, 2) *
      config.helixCutTolerance;
  result.pT2perRadius = std::pow(config.highland / result.pTPerHelixRadius, 2);
  result.sigmapT2perRadius =
      result.pT2perRadius * std::pow(2 * config.sigmaScattering, 2);
  result.multipleScattering2 =
      config.maxScatteringAngle2 * std::pow(config.sigmaScattering, 2);

  return result;
}

TripletSeedFinder2::TripletSeedFinder2(const DerivedConfig& config,
                                       std::unique_ptr<const Logger> logger_)
    : m_cfg(config), m_logger(std::move(logger_)) {}

TripletSeedFinder2::DoubletCuts TripletSeedFinder2::deriveDoubletCuts(
    const DerivedOptions& options,
    SpacePointCandidateType spacePointCandidateType) const {
  DoubletCuts cuts;

  cuts.impactMax = m_cfg.impactMax;
  cuts.interactionPointCut = m_cfg.interactionPointCut;
  cuts.deltaRMin = (spacePointCandidateType == SpacePointCandidateType::eBottom)
                       ? m_cfg.deltaRMinBottomSP
                       : m_cfg.deltaRMinTopSP;
  cuts.deltaRMax = (spacePointCandidateType == SpacePointCandidateType::eBottom)
                       ? m_cfg.deltaRMaxBottomSP
                       : m_cfg.deltaRMaxTopSP;
  cuts.collisionRegionMin = m_cfg.collisionRegionMin;
  cuts.collisionRegionMax = m_cfg.collisionRegionMax;
  cuts.cotThetaMax = m_cfg.cotThetaMax;
  cuts.deltaZMax = m_cfg.deltaZMax;
  cuts.minHelixDiameter2 = options.minHelixDiameter2;
  cuts.experimentCuts = m_cfg.experimentCuts;

  return cuts;
}

TripletSeedFinder2::TripletCuts TripletSeedFinder2::deriveTripletCuts(
    const DerivedOptions& options) const {
  TripletCuts cuts;

  cuts.minPt = m_cfg.minPt;
  cuts.sigmaScattering = m_cfg.sigmaScattering;
  cuts.radLengthPerSeed = m_cfg.radLengthPerSeed;
  cuts.maxPtScattering = m_cfg.maxPtScattering;
  cuts.impactMax = m_cfg.impactMax;
  cuts.helixCutTolerance = m_cfg.helixCutTolerance;
  cuts.highland = m_cfg.highland;
  cuts.pTPerHelixRadius = options.pTPerHelixRadius;
  cuts.sigmapT2perRadius = options.sigmapT2perRadius;
  cuts.multipleScattering2 = options.multipleScattering2;
  cuts.minHelixDiameter2 = options.minHelixDiameter2;
  cuts.seedConfirmation = m_cfg.seedConfirmation;
  cuts.centralSeedConfirmationRange = m_cfg.centralSeedConfirmationRange;
  cuts.forwardSeedConfirmationRange = m_cfg.forwardSeedConfirmationRange;
  cuts.toleranceParam = m_cfg.toleranceParam;

  return cuts;
}

TripletSeedFinder2::MiddleSpacePointInfo
TripletSeedFinder2::computeMiddleSpacePointInfo(
    const ConstSpacePointProxy2& spM,
    const SpacePointContainer2::DenseColumn<float>& rColumn) {
  const float rM = spM.extra(rColumn);
  const float uIP = -1. / rM;
  const float cosPhiM = -spM.x() * uIP;
  const float sinPhiM = -spM.y() * uIP;
  const float uIP2 = uIP * uIP;

  return {uIP, uIP2, cosPhiM, sinPhiM};
}

std::pair<float, float> TripletSeedFinder2::retrieveRadiusRangeForMiddle(
    const ConstSpacePointProxy2& spM,
    const Range1D<float>& rMiddleSpRange) const {
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

void TripletSeedFinder2::initialize(State& state,
                                    const DerivedOptions& options) const {
  state.options = options;

  state.bottomDoubletCuts =
      deriveDoubletCuts(options, SpacePointCandidateType::eBottom);
  state.topDoubletCuts =
      deriveDoubletCuts(options, SpacePointCandidateType::eTop);

  state.tripletCuts = deriveTripletCuts(options);
}

void TripletSeedFinder2::createSeeds(
    State& state, Cache& cache, const ContainerPointers& containerPointers,
    std::span<const SpacePointIndex2> bottomSps,
    std::span<const SpacePointIndex2> middleSps,
    std::span<const SpacePointIndex2> topSps,
    SeedContainer2& outputSeeds) const {
  TripletSeedFilter2::Options filterOptions;
  filterOptions.seedConfirmation = m_cfg.seedConfirmation;

  cache.candidatesCollector.reserve(m_cfg.maxSeedsPerSpMConf,
                                    m_cfg.maxQualitySeedsPerSpMConf);

  if (middleSps.empty()) {
    ACTS_VERBOSE("No middle space points, skipping");
    return;
  }

  auto firstMiddleSp = containerPointers.spacePoints().at(middleSps.front());

  // we compute this here since all middle space point candidates belong to
  // the same z-bin
  auto [minRadiusRangeForMiddle, maxRadiusRangeForMiddle] =
      retrieveRadiusRangeForMiddle(firstMiddleSp, state.options.rMiddleSpRange);
  ACTS_VERBOSE("Validity range (radius) for the middle space point is ["
               << minRadiusRangeForMiddle << ", " << maxRadiusRangeForMiddle
               << "]");

  for (SpacePointIndex2 middleSpIndex : middleSps) {
    auto spM = containerPointers.spacePoints().at(middleSpIndex);

    const float rM = spM.extra(containerPointers.rColumn());
    const float zM = spM.z();

    // check if spM is outside our radial region of interest
    if (rM < minRadiusRangeForMiddle) {
      continue;
    }
    if (rM > maxRadiusRangeForMiddle) {
      // break because SPs are sorted in r
      break;
    }

    MiddleSpacePointInfo middleSpacePointInfo =
        computeMiddleSpacePointInfo(spM, containerPointers.rColumn());

    // Iterate over middle-top doublets
    cache.topDoublets.clear();
    createDoublets<SpacePointCandidateType::eTop>(
        state.topDoubletCuts, containerPointers, spM, middleSpacePointInfo,
        topSps, cache.topDoublets);

    // no top SP found -> try next spM
    if (cache.topDoublets.empty()) {
      ACTS_VERBOSE("No compatible Tops, moving to next middle candidate");
      continue;
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
      // continue if number of top SPs is smaller than minimum
      if (cache.topDoublets.size() < filterOptions.nTopSeedConf) {
        ACTS_VERBOSE("Number of top SPs is "
                     << cache.topDoublets.size()
                     << " and is smaller than minimum, moving to next middle "
                        "candidate");
        continue;
      }
    }

    // Iterate over middle-bottom doublets
    cache.bottomDoublets.clear();
    createDoublets<SpacePointCandidateType::eBottom>(
        state.bottomDoubletCuts, containerPointers, spM, middleSpacePointInfo,
        bottomSps, cache.bottomDoublets);

    // no bottom SP found -> try next spM
    if (cache.bottomDoublets.empty()) {
      ACTS_VERBOSE("No compatible Bottoms, moving to next middle candidate");
      continue;
    }

    ACTS_VERBOSE("Candidates: " << cache.bottomDoublets.size()
                                << " bottoms and " << cache.topDoublets.size()
                                << " tops for middle candidate indexed "
                                << spM.index());

    // filter candidates
    cache.candidatesCollector.clear();
    if (m_cfg.useDetailedDoubleMeasurementInfo) {
      createTriplets<MeasurementInfo::eDetailed>(
          cache.tripletCache, state.tripletCuts, *m_cfg.filter, filterOptions,
          state.filter, cache.filter, containerPointers, spM,
          cache.bottomDoublets, cache.topDoublets, cache.tripletTopCandidates,
          cache.candidatesCollector);
    } else {
      createTriplets<MeasurementInfo::eDefault>(
          cache.tripletCache, state.tripletCuts, *m_cfg.filter, filterOptions,
          state.filter, cache.filter, containerPointers, spM,
          cache.bottomDoublets, cache.topDoublets, cache.tripletTopCandidates,
          cache.candidatesCollector);
    }

    // retrieve all candidates
    // this collection is already sorted, higher weights first
    const std::size_t numQualitySeeds =
        cache.candidatesCollector.nHighQualityCandidates();
    cache.candidatesCollector.toSortedCandidates(
        containerPointers.spacePoints(), cache.sortedCandidates);
    m_cfg.filter->filter1SpFixed(
        filterOptions, state.filter, containerPointers.spacePoints(),
        containerPointers.copiedFromIndexColumn, cache.sortedCandidates,
        numQualitySeeds, outputSeeds);
  }  // loop on middle space points
}

template <TripletSeedFinder2::SpacePointCandidateType candidate_type>
void TripletSeedFinder2::createDoublets(
    const DoubletCuts& cuts, const ContainerPointers& containerPointers,
    const ConstSpacePointProxy2& middleSp,
    const MiddleSpacePointInfo& middleSpInfo,
    std::span<const SpacePointIndex2> candidateSps,
    Doublets& compatibleDoublets) {
  constexpr bool isBottomCandidate =
      candidate_type == SpacePointCandidateType::eBottom;

  const float impactMax = isBottomCandidate ? -cuts.impactMax : cuts.impactMax;

  const float xM = middleSp.x();
  const float yM = middleSp.y();
  const float zM = middleSp.z();
  const float rM = middleSp.extra(containerPointers.rColumn());

  float vIPAbs = 0;
  if (cuts.interactionPointCut) {
    // equivalent to m_cfg.impactMax / (rM * rM);
    vIPAbs = impactMax * middleSpInfo.uIP2;
  }

  float deltaR = 0.;
  float deltaZ = 0.;

  const auto outsideRangeCheck = [](const float value, const float min,
                                    const float max) -> bool {
    // intentionally using `|` after profiling. faster due to better branch
    // prediction
    return static_cast<bool>(static_cast<int>(value < min) |
                             static_cast<int>(value > max));
  };

  const auto calculateError = [&](const ConstSpacePointProxy2& otherSp,
                                  float iDeltaR2, float cotTheta) {
    // TOD use some reasonable defaults

    float varianceZM = containerPointers.hasVarianceColumns()
                           ? middleSp.extra(containerPointers.varianceZColumn())
                           : 0;
    float varianceZO = containerPointers.hasVarianceColumns()
                           ? otherSp.extra(containerPointers.varianceZColumn())
                           : 0;
    float varianceRM = containerPointers.hasVarianceColumns()
                           ? middleSp.extra(containerPointers.varianceRColumn())
                           : 0;
    float varianceRO = containerPointers.hasVarianceColumns()
                           ? otherSp.extra(containerPointers.varianceRColumn())
                           : 0;

    return iDeltaR2 * ((varianceZM + varianceZO) +
                       (cotTheta * cotTheta) * (varianceRM + varianceRO));
  };

  for (SpacePointIndex2 otherSpIndex : candidateSps) {
    ConstSpacePointProxy2 otherSp =
        containerPointers.spacePoints().at(otherSpIndex);

    if constexpr (isBottomCandidate) {
      deltaR = rM - otherSp.extra(containerPointers.rColumn());
    } else {
      deltaR = otherSp.extra(containerPointers.rColumn()) - rM;
    }

    if (outsideRangeCheck(deltaR, cuts.deltaRMin, cuts.deltaRMax)) {
      continue;
    }

    if constexpr (isBottomCandidate) {
      deltaZ = zM - otherSp.z();
    } else {
      deltaZ = otherSp.z() - zM;
    }

    // the longitudinal impact parameter zOrigin is defined as (zM - rM *
    // cotTheta) where cotTheta is the ratio Z/R (forward angle) of space
    // point duplet but instead we calculate (zOrigin * deltaR) and multiply
    // collisionRegion by deltaR to avoid divisions
    const float zOriginTimesDeltaR = zM * deltaR - rM * deltaZ;
    // check if duplet origin on z axis within collision region
    if (outsideRangeCheck(zOriginTimesDeltaR, cuts.collisionRegionMin * deltaR,
                          cuts.collisionRegionMax * deltaR)) {
      continue;
    }

    // if interactionPointCut is false we apply z cuts before coordinate
    // transformation to avoid unnecessary calculations. If
    // interactionPointCut is true we apply the curvature cut first because it
    // is more frequent but requires the coordinate transformation
    if (!cuts.interactionPointCut) {
      // check if duplet cotTheta is within the region of interest
      // cotTheta is defined as (deltaZ / deltaR) but instead we multiply
      // cotThetaMax by deltaR to avoid division
      if (outsideRangeCheck(deltaZ, -cuts.cotThetaMax * deltaR,
                            cuts.cotThetaMax * deltaR)) {
        continue;
      }
      // if z-distance between SPs is within max and min values
      if (outsideRangeCheck(deltaZ, -cuts.deltaZMax, cuts.deltaZMax)) {
        continue;
      }

      // transform SP coordinates to the u-v reference frame
      const float deltaX = otherSp.x() - xM;
      const float deltaY = otherSp.y() - yM;

      const float xNewFrame =
          deltaX * middleSpInfo.cosPhiM + deltaY * middleSpInfo.sinPhiM;
      const float yNewFrame =
          deltaY * middleSpInfo.cosPhiM - deltaX * middleSpInfo.sinPhiM;

      const float deltaR2 = deltaX * deltaX + deltaY * deltaY;
      const float iDeltaR2 = 1. / deltaR2;

      const float uT = xNewFrame * iDeltaR2;
      const float vT = yNewFrame * iDeltaR2;

      const float iDeltaR = std::sqrt(iDeltaR2);
      const float cotTheta = deltaZ * iDeltaR;

      const float er = calculateError(otherSp, iDeltaR2, cotTheta);

      // fill output vectors
      compatibleDoublets.emplace_back(
          otherSp.index(),
          {cotTheta, iDeltaR, er, uT, vT, xNewFrame, yNewFrame});
      continue;
    }

    // transform SP coordinates to the u-v reference frame
    const float deltaX = otherSp.x() - xM;
    const float deltaY = otherSp.y() - yM;

    const float xNewFrame =
        deltaX * middleSpInfo.cosPhiM + deltaY * middleSpInfo.sinPhiM;
    const float yNewFrame =
        deltaY * middleSpInfo.cosPhiM - deltaX * middleSpInfo.sinPhiM;

    const float deltaR2 = deltaX * deltaX + deltaY * deltaY;
    const float iDeltaR2 = 1. / deltaR2;

    const float uT = xNewFrame * iDeltaR2;
    const float vT = yNewFrame * iDeltaR2;

    // We check the interaction point by evaluating the minimal distance
    // between the origin and the straight line connecting the two points in
    // the doublets. Using a geometric similarity, the Im is given by
    // yNewFrame * rM / deltaR <= config.impactMax
    // However, we make here an approximation of the impact parameter
    // which is valid under the assumption yNewFrame / xNewFrame is small
    // The correct computation would be:
    // yNewFrame * yNewFrame * rM * rM <= config.impactMax *
    // config.impactMax * deltaR2
    if (std::abs(rM * yNewFrame) <= impactMax * xNewFrame) {
      // check if duplet cotTheta is within the region of interest
      // cotTheta is defined as (deltaZ / deltaR) but instead we multiply
      // cotThetaMax by deltaR to avoid division
      if (outsideRangeCheck(deltaZ, -cuts.cotThetaMax * deltaR,
                            cuts.cotThetaMax * deltaR)) {
        continue;
      }

      const float iDeltaR = std::sqrt(iDeltaR2);
      const float cotTheta = deltaZ * iDeltaR;

      // discard bottom-middle doublets in a certain (r, eta) region according
      // to detector specific cuts
      if constexpr (isBottomCandidate) {
        if (cuts.experimentCuts.connected() &&
            !cuts.experimentCuts(otherSp.extra(containerPointers.rColumn()),
                                 cotTheta)) {
          continue;
        }
      }

      const float er = calculateError(otherSp, iDeltaR2, cotTheta);

      // fill output vectors
      compatibleDoublets.emplace_back(
          otherSp.index(),
          {cotTheta, iDeltaR, er, uT, vT, xNewFrame, yNewFrame});
      continue;
    }

    // in the rotated frame the interaction point is positioned at x = -rM
    // and y ~= impactParam
    const float vIP = (yNewFrame > 0) ? -vIPAbs : vIPAbs;

    // we can obtain aCoef as the slope dv/du of the linear function,
    // estimated using du and dv between the two SP bCoef is obtained by
    // inserting aCoef into the linear equation
    const float aCoef = (vT - vIP) / (uT - middleSpInfo.uIP);
    const float bCoef = vIP - aCoef * middleSpInfo.uIP;
    // the distance of the straight line from the origin (radius of the
    // circle) is related to aCoef and bCoef by d^2 = bCoef^2 / (1 +
    // aCoef^2) = 1 / (radius^2) and we can apply the cut on the curvature
    if ((bCoef * bCoef) * cuts.minHelixDiameter2 > 1 + aCoef * aCoef) {
      continue;
    }

    // check if duplet cotTheta is within the region of interest
    // cotTheta is defined as (deltaZ / deltaR) but instead we multiply
    // cotThetaMax by deltaR to avoid division
    if (outsideRangeCheck(deltaZ, -cuts.cotThetaMax * deltaR,
                          cuts.cotThetaMax * deltaR)) {
      continue;
    }

    const float iDeltaR = std::sqrt(iDeltaR2);
    const float cotTheta = deltaZ * iDeltaR;

    // discard bottom-middle doublets in a certain (r, eta) region according
    // to detector specific cuts
    if constexpr (isBottomCandidate) {
      if (cuts.experimentCuts.connected() &&
          !cuts.experimentCuts(otherSp.extra(containerPointers.rColumn()),
                               cotTheta)) {
        continue;
      }
    }

    const float er = calculateError(otherSp, iDeltaR2, cotTheta);

    // fill output vectors
    compatibleDoublets.emplace_back(
        otherSp.index(), {cotTheta, iDeltaR, er, uT, vT, xNewFrame, yNewFrame});
  }
}

template <TripletSeedFinder2::MeasurementInfo measurement_info>
void TripletSeedFinder2::createTriplets(
    TripletCache& cache, const TripletCuts& cuts,
    const TripletSeedFilter2& filter,
    const TripletSeedFilter2::Options& filterOptions,
    TripletSeedFilter2::State& filterState,
    TripletSeedFilter2::Cache& filterCache,
    const ContainerPointers& containerPointers,
    const ConstSpacePointProxy2& spM, const Doublets& bottomDoublets,
    const Doublets& topDoublets, TripletTopCandidates& tripletTopCandidates,
    CandidatesForMiddleSp2& candidatesCollector) {
  constexpr bool isDetailedMeasurement =
      measurement_info == MeasurementInfo::eDetailed;

  const float rM = spM.extra(containerPointers.rColumn());
  const float cosPhiM = spM.x() / rM;
  const float sinPhiM = spM.y() / rM;
  const float varianceRM = containerPointers.hasVarianceColumns()
                               ? spM.extra(containerPointers.varianceRColumn())
                               : 0;
  const float varianceZM = containerPointers.hasVarianceColumns()
                               ? spM.extra(containerPointers.varianceZColumn())
                               : 0;

  // make index vectors for sorting
  cache.sortedBottoms.resize(bottomDoublets.size());
  std::iota(cache.sortedBottoms.begin(), cache.sortedBottoms.end(), 0);

  cache.sortedTops.resize(topDoublets.size());
  std::iota(cache.sortedTops.begin(), cache.sortedTops.end(), 0);

  if constexpr (!isDetailedMeasurement) {
    // sorting becomes less expensive when we copy the cotTheta values into
    // their own arrays due to more optimal usage of the cache

    std::ranges::sort(cache.sortedBottoms, {},
                      [&bottomDoublets](const std::size_t s) {
                        return bottomDoublets.cotTheta[s];
                      });
    std::ranges::sort(cache.sortedTops, {},
                      [&topDoublets](const std::size_t s) {
                        return topDoublets.cotTheta[s];
                      });
  }

  // Reserve enough space, in case current capacity is too little
  tripletTopCandidates.resize(topDoublets.size());

  std::size_t t0 = 0;

  for (const std::size_t b : cache.sortedBottoms) {
    // break if we reached the last top SP
    if (t0 >= topDoublets.size()) {
      break;
    }

    auto spB =
        containerPointers.spacePoints().at(bottomDoublets.spacePoints[b]);

    const auto& lb = bottomDoublets.linCircles[b];
    float cotThetaB = lb.cotTheta;
    float Vb = lb.V;
    float Ub = lb.U;
    float ErB = lb.Er;
    float iDeltaRB = lb.iDeltaR;

    // 1+(cot^2(theta)) = 1/sin^2(theta)
    float iSinTheta2 = 1. + cotThetaB * cotThetaB;
    float sigmaSquaredPtDependent = iSinTheta2 * cuts.sigmapT2perRadius;
    // calculate max scattering for min momentum at the seed's theta angle
    // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
    // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
    // scattering
    // but to avoid trig functions we approximate cot by scaling by
    // 1/sin^4(theta)
    // resolving with pT to p scaling --> only divide by sin^2(theta)
    // max approximation error for allowed scattering angles of 0.04 rad at
    // eta=infinity: ~8.5%
    float scatteringInRegion2 = cuts.multipleScattering2 * iSinTheta2;

    float sinTheta = 1 / std::sqrt(iSinTheta2);
    float cosTheta = cotThetaB * sinTheta;

    // clear all vectors used in each inner for loop
    tripletTopCandidates.clear();

    // coordinate transformation and checks for middle spacepoint
    // x and y terms for the rotation from UV to XY plane
    Eigen::Vector2f rotationTermsUVtoXY = {0, 0};
    if constexpr (isDetailedMeasurement) {
      rotationTermsUVtoXY[0] = cosPhiM * sinTheta;
      rotationTermsUVtoXY[1] = sinPhiM * sinTheta;
    }

    // minimum number of compatible top SPs to trigger the filter for a certain
    // middle bottom pair if seedConfirmation is false we always ask for at
    // least one compatible top to trigger the filter
    std::size_t minCompatibleTopSPs = 2;
    if (!cuts.seedConfirmation ||
        spB.extra(containerPointers.rColumn()) >
            filterOptions.seedConfRange.rMaxSeedConf) {
      minCompatibleTopSPs = 1;
    }
    if (cuts.seedConfirmation &&
        candidatesCollector.nHighQualityCandidates() > 0) {
      minCompatibleTopSPs++;
    }

    for (std::size_t indexSortedTop = t0; indexSortedTop < topDoublets.size();
         ++indexSortedTop) {
      const std::size_t t = cache.sortedTops[indexSortedTop];

      auto spT = containerPointers.spacePoints().at(topDoublets.spacePoints[t]);

      const auto& lt = topDoublets.linCircles[t];

      float cotThetaT = lt.cotTheta;
      float rMxy = 0.;
      float ub = 0.;
      float vb = 0.;
      float ut = 0.;
      float vt = 0.;
      Vector3 rMTransf;
      float xB = 0.;
      float yB = 0.;
      float xT = 0.;
      float yT = 0.;
      float iDeltaRB2 = 0.;
      float iDeltaRT2 = 0.;

      if constexpr (isDetailedMeasurement) {
        // protects against division by 0
        float dU = lt.U - Ub;
        if (dU == 0.) {
          continue;
        }
        // A and B are evaluated as a function of the circumference parameters
        // x_0 and y_0
        float A0 = (lt.V - Vb) / dU;

        float zPositionMiddle = cosTheta * std::sqrt(1 + A0 * A0);

        // position of Middle SP converted from UV to XY assuming cotTheta
        // evaluated from the Bottom and Middle SPs double
        Vector3 positionMiddle = {
            rotationTermsUVtoXY[0] - rotationTermsUVtoXY[1] * A0,
            rotationTermsUVtoXY[0] * A0 + rotationTermsUVtoXY[1],
            zPositionMiddle};

        if (!stripCoordinateCheck(cuts.toleranceParam, spM, containerPointers,
                                  positionMiddle, rMTransf)) {
          continue;
        }

        // coordinate transformation and checks for bottom spacepoint
        float B0 = 2. * (Vb - A0 * Ub);
        float Cb = 1. - B0 * lb.y;
        float Sb = A0 + B0 * lb.x;
        Vector3 positionBottom = {
            rotationTermsUVtoXY[0] * Cb - rotationTermsUVtoXY[1] * Sb,
            rotationTermsUVtoXY[0] * Sb + rotationTermsUVtoXY[1] * Cb,
            zPositionMiddle};

        Vector3 rBTransf;
        if (!stripCoordinateCheck(cuts.toleranceParam, spB, containerPointers,
                                  positionBottom, rBTransf)) {
          continue;
        }

        // coordinate transformation and checks for top spacepoint
        float Ct = 1. - B0 * lt.y;
        float St = A0 + B0 * lt.x;
        Vector3 positionTop = {
            rotationTermsUVtoXY[0] * Ct - rotationTermsUVtoXY[1] * St,
            rotationTermsUVtoXY[0] * St + rotationTermsUVtoXY[1] * Ct,
            zPositionMiddle};

        Vector3 rTTransf;
        if (!stripCoordinateCheck(cuts.toleranceParam, spT, containerPointers,
                                  positionTop, rTTransf)) {
          continue;
        }

        // bottom and top coordinates in the spM reference frame
        xB = rBTransf[0] - rMTransf[0];
        yB = rBTransf[1] - rMTransf[1];
        float zB = rBTransf[2] - rMTransf[2];
        xT = rTTransf[0] - rMTransf[0];
        yT = rTTransf[1] - rMTransf[1];
        float zT = rTTransf[2] - rMTransf[2];

        iDeltaRB2 = 1. / (xB * xB + yB * yB);
        iDeltaRT2 = 1. / (xT * xT + yT * yT);

        cotThetaB = -zB * std::sqrt(iDeltaRB2);
        cotThetaT = zT * std::sqrt(iDeltaRT2);
      }

      // use geometric average
      float cotThetaAvg2 = cotThetaB * cotThetaT;
      if constexpr (isDetailedMeasurement) {
        // use arithmetic average
        float averageCotTheta = 0.5 * (cotThetaB + cotThetaT);
        cotThetaAvg2 = averageCotTheta * averageCotTheta;
      } else if (cotThetaAvg2 <= 0) {
        continue;
      }

      // add errors of spB-spM and spM-spT pairs and add the correlation term
      // for errors on spM
      float error2 =
          lt.Er + ErB +
          2 * (cotThetaAvg2 * varianceRM + varianceZM) * iDeltaRB * lt.iDeltaR;

      float deltaCotTheta = cotThetaB - cotThetaT;
      float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;

      // Apply a cut on the compatibility between the r-z slope of the two
      // seed segments. This is done by comparing the squared difference
      // between slopes, and comparing to the squared uncertainty in this
      // difference - we keep a seed if the difference is compatible within
      // the assumed uncertainties. The uncertainties get contribution from
      // the  space-point-related squared error (error2) and a scattering term
      // calculated assuming the minimum pt we expect to reconstruct
      // (scatteringInRegion2). This assumes gaussian error propagation which
      // allows just adding the two errors if they are uncorrelated (which is
      // fair for scattering and measurement uncertainties)
      if (deltaCotTheta2 > error2 + scatteringInRegion2) {
        // skip top SPs based on cotTheta sorting when producing triplets
        if constexpr (isDetailedMeasurement) {
          continue;
        }
        // break if cotTheta from bottom SP < cotTheta from top SP because
        // the SP are sorted by cotTheta
        if (cotThetaB < cotThetaT) {
          break;
        }
        t0 = indexSortedTop + 1;
        continue;
      }

      if constexpr (isDetailedMeasurement) {
        rMxy = std::sqrt(rMTransf[0] * rMTransf[0] + rMTransf[1] * rMTransf[1]);
        double irMxy = 1 / rMxy;
        float Ax = rMTransf[0] * irMxy;
        float Ay = rMTransf[1] * irMxy;

        ub = (xB * Ax + yB * Ay) * iDeltaRB2;
        vb = (yB * Ax - xB * Ay) * iDeltaRB2;
        ut = (xT * Ax + yT * Ay) * iDeltaRT2;
        vt = (yT * Ax - xT * Ay) * iDeltaRT2;
      }

      float dU = 0;
      float A = 0;
      float S2 = 0;
      float B = 0;
      float B2 = 0;

      if constexpr (isDetailedMeasurement) {
        dU = ut - ub;
        // protects against division by 0
        if (dU == 0.) {
          continue;
        }
        A = (vt - vb) / dU;
        S2 = 1. + A * A;
        B = vb - A * ub;
        B2 = B * B;
      } else {
        dU = lt.U - Ub;
        // protects against division by 0
        if (dU == 0.) {
          continue;
        }
        // A and B are evaluated as a function of the circumference parameters
        // x_0 and y_0
        A = (lt.V - Vb) / dU;
        S2 = 1. + A * A;
        B = Vb - A * Ub;
        B2 = B * B;
      }

      // sqrt(S2)/B = 2 * helixradius
      // calculated radius must not be smaller than minimum radius
      if (S2 < B2 * cuts.minHelixDiameter2) {
        continue;
      }

      // refinement of the cut on the compatibility between the r-z slope of
      // the two seed segments using a scattering term scaled by the actual
      // measured pT (p2scatterSigma)
      float iHelixDiameter2 = B2 / S2;
      // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
      // from rad to deltaCotTheta
      float p2scatterSigma = iHelixDiameter2 * sigmaSquaredPtDependent;
      if (!std::isinf(cuts.maxPtScattering)) {
        // if pT > maxPtScattering, calculate allowed scattering angle using
        // maxPtScattering instead of pt.
        // To avoid 0-divison the pT check is skipped in case of B2==0, and
        // p2scatterSigma is calculated directly from maxPtScattering
        if (B2 == 0 || cuts.pTPerHelixRadius * std::sqrt(S2 / B2) >
                           2. * cuts.maxPtScattering) {
          float pTscatterSigma =
              (cuts.highland / cuts.maxPtScattering) * cuts.sigmaScattering;
          p2scatterSigma = pTscatterSigma * pTscatterSigma * iSinTheta2;
        }
      }

      // if deltaTheta larger than allowed scattering for calculated pT, skip
      if (deltaCotTheta2 > error2 + p2scatterSigma) {
        if constexpr (isDetailedMeasurement) {
          continue;
        }
        if (cotThetaB < cotThetaT) {
          break;
        }
        t0 = indexSortedTop;
        continue;
      }
      // A and B allow calculation of impact params in U/V plane with linear
      // function
      // (in contrast to having to solve a quadratic function in x/y plane)
      float Im = 0;
      if constexpr (isDetailedMeasurement) {
        Im = std::abs((A - B * rMxy) * rMxy);
      } else {
        Im = std::abs((A - B * rM) * rM);
      }

      if (Im > cuts.impactMax) {
        continue;
      }

      // inverse diameter is signed depending on if the curvature is
      // positive/negative in phi
      tripletTopCandidates.emplace_back(spT.index(), B / std::sqrt(S2), Im);
    }  // loop on tops

    // continue if number of top SPs is smaller than minimum required for filter
    if (tripletTopCandidates.topSpacePoints.size() < minCompatibleTopSPs) {
      continue;
    }

    float zOrigin = spM.z() - rM * lb.cotTheta;
    filter.filter2SpFixed(
        filterOptions, filterState, filterCache,
        containerPointers.spacePoints(), containerPointers.rColumn(),
        containerPointers.copiedFromIndexColumn, spB.index(), spM.index(),
        tripletTopCandidates.topSpacePoints, tripletTopCandidates.curvatures,
        tripletTopCandidates.impactParameters, zOrigin, candidatesCollector);
  }  // loop on bottoms
}

bool TripletSeedFinder2::stripCoordinateCheck(
    double tolerance, const ConstSpacePointProxy2& sp,
    const ContainerPointers& containerPointers,
    const Vector3& spacePointPosition, Vector3& outputCoordinates) {
  const Vector3& topStripVector =
      sp.extra(containerPointers.topStripVectorColumn());
  const Vector3& bottomStripVector =
      sp.extra(containerPointers.bottomStripVectorColumn());
  const Vector3& stripCenterDistance =
      sp.extra(containerPointers.stripCenterDistanceColumn());

  // cross product between top strip vector and spacepointPosition
  Vector3 d1 = topStripVector.cross(spacePointPosition);

  // scalar product between bottom strip vector and d1
  double bd1 = bottomStripVector.dot(d1);

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the bottom detector element
  double s1 = stripCenterDistance.dot(d1);
  if (std::abs(s1) > std::abs(bd1) * tolerance) {
    return false;
  }

  // cross product between bottom strip vector and spacepointPosition
  Vector3 d0 = bottomStripVector.cross(spacePointPosition);

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the top detector element
  double s0 = stripCenterDistance.dot(d0);
  if (std::abs(s0) > std::abs(bd1) * tolerance) {
    return false;
  }

  // if arrive here spacepointPosition is compatible with strip directions and
  // detector elements

  const Vector3& topStripCenterPosition =
      sp.extra(containerPointers.topStripCenterPositionColumn());

  // spacepointPosition corrected with respect to the top strip position and
  // direction and the distance between the strips
  s0 = s0 / bd1;
  outputCoordinates = topStripCenterPosition + topStripVector * s0;
  return true;
}

}  // namespace Acts
