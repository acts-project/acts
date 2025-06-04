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

namespace Acts {

using namespace UnitLiterals;

TripletSeedFinder2::DerivedConfig TripletSeedFinder2::Config::derive() const {
  DerivedConfig result;

  static_cast<Config&>(result) = *this;

  // TODO get rid of unit conversions
  {
    result.minPt /= 1_MeV;
    result.deltaRMinTopSP /= 1_mm;
    result.deltaRMaxTopSP /= 1_mm;
    result.deltaRMinBottomSP /= 1_mm;
    result.deltaRMaxBottomSP /= 1_mm;
    result.impactMax /= 1_mm;
    result.maxPtScattering /= 1_MeV;
    result.collisionRegionMin /= 1_mm;
    result.collisionRegionMax /= 1_mm;
    result.deltaZMax /= 1_mm;
    result.toleranceParam /= 1_mm;
  }

  // TODO use Interactions.hpp for highland

  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  result.highland = 13.6 * std::sqrt(result.radLengthPerSeed) *
                    (1 + 0.038 * std::log(result.radLengthPerSeed));
  const float maxScatteringAngle = result.highland / result.minPt;
  result.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;

  return result;
}

TripletSeedFinder2::DerivedOptions TripletSeedFinder2::Options::derive(
    const DerivedConfig& config) const {
  DerivedOptions result;

  static_cast<Options&>(result) = *this;

  // TODO get rid of unit conversions
  {
    result.beamPos[0] /= 1_mm;
    result.beamPos[1] /= 1_mm;
    result.bFieldInZ /= 1000. * 1_T;
  }

  result.pTPerHelixRadius = 1_T * 1e6 * result.bFieldInZ;
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

void TripletSeedFinder2::deriveDoubletCuts(
    DoubletCuts& cuts, const ConstSpacePointProxy2& spM,
    const SpacePointContainer2::DenseColumn<float>& rColumn) const {
  const float rM = spM.extra(rColumn);
  const float uIP = -1. / rM;
  const float cosPhiM = -spM.x() * uIP;
  const float sinPhiM = -spM.y() * uIP;
  const float uIP2 = uIP * uIP;

  cuts.uIP = uIP;
  cuts.uIP2 = uIP2;
  cuts.cosPhiM = cosPhiM;
  cuts.sinPhiM = sinPhiM;
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

void TripletSeedFinder2::createSeeds(
    const DerivedOptions& options, State& state,
    const SpacePointContainer2& spacePoints,
    const SpacePointContainer2::DenseColumn<float>& rColumn,
    const SpacePointContainer2::DenseColumn<float>* varianceRColumn,
    const SpacePointContainer2::DenseColumn<float>* varianceZColumn,
    std::vector<SpacePointIndex2>& bottomSps,
    std::vector<SpacePointIndex2>& middleSps,
    std::vector<SpacePointIndex2>& topSps, SeedContainer2& outputSeeds) const {
  std::ranges::sort(bottomSps, {}, [&](SpacePointIndex2 index) {
    return -spacePoints.at(index).z() / spacePoints.at(index).extra(rColumn);
  });
  std::ranges::sort(middleSps, {}, [&](SpacePointIndex2 index) {
    return spacePoints.at(index).extra(rColumn);
  });
  std::ranges::sort(topSps, {}, [&](SpacePointIndex2 index) {
    return spacePoints.at(index).z() / spacePoints.at(index).extra(rColumn);
  });

  TripletSeedFilter2::Options filterOptions;
  filterOptions.seedConfirmation = m_cfg.seedConfirmation;

  state.candidatesCollector.reserve(m_cfg.maxSeedsPerSpMConf,
                                    m_cfg.maxQualitySeedsPerSpMConf);

  if (middleSps.empty()) {
    ACTS_VERBOSE("No middle space points, skipping");
    return;
  }

  auto firstMiddleSp = spacePoints.at(middleSps.front());

  // we compute this here since all middle space point candidates belong to
  // the same z-bin
  auto [minRadiusRangeForMiddle, maxRadiusRangeForMiddle] =
      retrieveRadiusRangeForMiddle(firstMiddleSp, options.rMiddleSpRange);
  ACTS_VERBOSE("Validity range (radius) for the middle space point is ["
               << minRadiusRangeForMiddle << ", " << maxRadiusRangeForMiddle
               << "]");

  for (SpacePointIndex2 middleSpIndex : middleSps) {
    auto spM = spacePoints.at(middleSpIndex);

    const float rM = spM.extra(rColumn);
    const float zM = spM.z();

    // check if spM is outside our radial region of interest
    if (rM < minRadiusRangeForMiddle) {
      continue;
    }
    if (rM > maxRadiusRangeForMiddle) {
      // break because SPs are sorted in r
      break;
    }

    DoubletCuts doubletCuts;
    doubletCuts.deltaRMin = m_cfg.deltaRMinTopSP;
    doubletCuts.deltaRMax = m_cfg.deltaRMaxTopSP;
    deriveDoubletCuts(doubletCuts, spM, rColumn);

    // Iterate over middle-top doublets
    state.compatibleTopSp.clear();
    state.linCirclesTop.clear();
    state.linCircleCotThetaTop.clear();
    createDoublets<SpacePointCandidateType::eTop>(
        m_cfg, options, doubletCuts, spacePoints, rColumn, varianceRColumn,
        varianceZColumn, spM, topSps, state.compatibleTopSp,
        state.linCirclesTop, state.linCircleCotThetaTop);

    // no top SP found -> try next spM
    if (state.compatibleTopSp.empty()) {
      ACTS_VERBOSE("No compatible Tops, moving to next middle candidate");
      continue;
    }

    // apply cut on the number of top SP if seedConfirmation is true
    if (m_cfg.seedConfirmation) {
      // check if middle SP is in the central or forward region
      //
      // intentionally using `|` after profiling. faster due to better branch
      // prediction
      filterOptions.seedConfRange =
          (zM > m_cfg.centralSeedConfirmationRange.zMaxSeedConf |
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
      if (state.compatibleTopSp.size() < filterOptions.nTopSeedConf) {
        ACTS_VERBOSE("Number of top SPs is "
                     << state.compatibleTopSp.size()
                     << " and is smaller than minimum, moving to next middle "
                        "candidate");
        continue;
      }
    }

    doubletCuts.deltaRMin = m_cfg.deltaRMinBottomSP;
    doubletCuts.deltaRMax = m_cfg.deltaRMaxBottomSP;

    // Iterate over middle-bottom doublets
    state.compatibleBottomSp.clear();
    state.linCirclesBottom.clear();
    state.linCircleCotThetaBottom.clear();
    createDoublets<SpacePointCandidateType::eBottom>(
        m_cfg, options, doubletCuts, spacePoints, rColumn, varianceRColumn,
        varianceZColumn, spM, bottomSps, state.compatibleBottomSp,
        state.linCirclesBottom, state.linCircleCotThetaBottom);

    // no bottom SP found -> try next spM
    if (state.compatibleBottomSp.empty()) {
      ACTS_VERBOSE("No compatible Bottoms, moving to next middle candidate");
      continue;
    }

    ACTS_VERBOSE("Candidates: "
                 << state.compatibleBottomSp.size() << " bottoms and "
                 << state.compatibleTopSp.size()
                 << " tops for middle candidate indexed " << spM.index());

    // filter candidates
    state.candidatesCollector.clear();
    if (m_cfg.useDetailedDoubleMeasurementInfo) {
      createTriplets<MeasurementInfo::eDetailed>(
          m_cfg, options, filterOptions, state.filterState, spacePoints,
          rColumn, varianceRColumn, varianceZColumn, nullptr, nullptr, nullptr,
          nullptr, spM, state.compatibleBottomSp, state.linCirclesBottom,
          state.linCircleCotThetaBottom, state.compatibleTopSp,
          state.linCirclesTop, state.linCircleCotThetaTop, state.sortedBottoms,
          state.sortedTops, state.topSpVec, state.curvatures,
          state.impactParameters, state.candidatesCollector);
    } else {
      createTriplets<MeasurementInfo::eDefault>(
          m_cfg, options, filterOptions, state.filterState, spacePoints,
          rColumn, varianceRColumn, varianceZColumn, nullptr, nullptr, nullptr,
          nullptr, spM, state.compatibleBottomSp, state.linCirclesBottom,
          state.linCircleCotThetaBottom, state.compatibleTopSp,
          state.linCirclesTop, state.linCircleCotThetaTop, state.sortedBottoms,
          state.sortedTops, state.topSpVec, state.curvatures,
          state.impactParameters, state.candidatesCollector);
    }

    // retrieve all candidates
    // this collection is already sorted, higher weights first
    const std::size_t numQualitySeeds =
        state.candidatesCollector.nHighQualityCandidates();
    state.candidatesCollector.toSortedCandidates(spacePoints,
                                                 state.sortedCandidates);
    m_cfg.seedFilter->filter1SpFixed(filterOptions, state.filterState,
                                     state.sortedCandidates, numQualitySeeds,
                                     outputSeeds);
  }  // loop on middle space points
}

template <TripletSeedFinder2::SpacePointCandidateType candidate_type>
void TripletSeedFinder2::createDoublets(
    const DerivedConfig& config, const DerivedOptions& options,
    const DoubletCuts& cuts, const SpacePointContainer2& spacePoints,
    const SpacePointContainer2::DenseColumn<float>& rColumn,
    const SpacePointContainer2::DenseColumn<float>* varianceRColumn,
    const SpacePointContainer2::DenseColumn<float>* varianceZColumn,
    const ConstSpacePointProxy2& middleSp,
    const std::vector<SpacePointIndex2>& candidateSps,
    std::vector<SpacePointIndex2>& compatibleSp,
    std::vector<LinCircle>& linCircles, std::vector<float>& cotThetas) {
  constexpr bool isBottomCandidate =
      candidate_type == SpacePointCandidateType::eBottom;

  const float impactMax =
      isBottomCandidate ? -config.impactMax : config.impactMax;

  const float xM = middleSp.x();
  const float yM = middleSp.y();
  const float zM = middleSp.z();
  const float rM = middleSp.extra(rColumn);

  float vIPAbs = 0;
  if (config.interactionPointCut) {
    // equivalent to m_cfg.impactMax / (rM * rM);
    vIPAbs = impactMax * cuts.uIP2;
  }

  float deltaR = 0.;
  float deltaZ = 0.;

  auto calculateError = [&](const ConstSpacePointProxy2& otherSp,
                            float iDeltaR2, float cotTheta) {
    // TOD use some reasonable defaults

    float varianceZM =
        varianceZColumn != nullptr ? middleSp.extra(*varianceZColumn) : 0;
    float varianceZO =
        varianceZColumn != nullptr ? otherSp.extra(*varianceZColumn) : 0;
    float varianceRM =
        varianceRColumn != nullptr ? middleSp.extra(*varianceRColumn) : 0;
    float varianceRO =
        varianceRColumn != nullptr ? otherSp.extra(*varianceRColumn) : 0;

    return iDeltaR2 * ((varianceZM + varianceZO) +
                       (cotTheta * cotTheta) * (varianceRM + varianceRO));
  };

  for (SpacePointIndex2 otherSpIndex : candidateSps) {
    ConstSpacePointProxy2 otherSp = spacePoints.at(otherSpIndex);

    if constexpr (isBottomCandidate) {
      deltaR = rM - otherSp.extra(rColumn);
    } else {
      deltaR = otherSp.extra(rColumn) - rM;
    }

    if (deltaR<cuts.deltaRMin | deltaR> cuts.deltaRMax) {
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
    //
    // intentionally using `|` after profiling. faster due to better branch
    // prediction
    if (zOriginTimesDeltaR<config.collisionRegionMin * deltaR |
                           zOriginTimesDeltaR>
            config.collisionRegionMax *
        deltaR) {
      continue;
    }

    // if interactionPointCut is false we apply z cuts before coordinate
    // transformation to avoid unnecessary calculations. If
    // interactionPointCut is true we apply the curvature cut first because it
    // is more frequent but requires the coordinate transformation
    if (!config.interactionPointCut) {
      // check if duplet cotTheta is within the region of interest
      // cotTheta is defined as (deltaZ / deltaR) but instead we multiply
      // cotThetaMax by deltaR to avoid division
      //
      // intentionally using `|` after profiling. faster due to better branch
      // prediction
      if (deltaZ > config.cotThetaMax * deltaR |
          deltaZ < -config.cotThetaMax * deltaR) {
        continue;
      }
      // if z-distance between SPs is within max and min values
      //
      // intentionally using `|` after profiling. faster due to better branch
      // prediction
      if (deltaZ > config.deltaZMax | deltaZ < -config.deltaZMax) {
        continue;
      }

      // transform SP coordinates to the u-v reference frame
      const float deltaX = otherSp.x() - xM;
      const float deltaY = otherSp.y() - yM;

      const float xNewFrame = deltaX * cuts.cosPhiM + deltaY * cuts.sinPhiM;
      const float yNewFrame = deltaY * cuts.cosPhiM - deltaX * cuts.sinPhiM;

      const float deltaR2 = deltaX * deltaX + deltaY * deltaY;
      const float iDeltaR2 = 1. / deltaR2;

      const float uT = xNewFrame * iDeltaR2;
      const float vT = yNewFrame * iDeltaR2;

      const float iDeltaR = std::sqrt(iDeltaR2);
      const float cotTheta = deltaZ * iDeltaR;

      const float er = calculateError(otherSp, iDeltaR2, cotTheta);

      // fill output vectors
      compatibleSp.emplace_back(otherSp.index());
      linCircles.emplace_back(cotTheta, iDeltaR, er, uT, vT, xNewFrame,
                              yNewFrame);
      cotThetas.emplace_back(cotTheta);
      continue;
    }

    // transform SP coordinates to the u-v reference frame
    const float deltaX = otherSp.x() - xM;
    const float deltaY = otherSp.y() - yM;

    const float xNewFrame = deltaX * cuts.cosPhiM + deltaY * cuts.sinPhiM;
    const float yNewFrame = deltaY * cuts.cosPhiM - deltaX * cuts.sinPhiM;

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
      //
      // intentionally using `|` after profiling. faster due to better branch
      // prediction
      if (deltaZ > config.cotThetaMax * deltaR |
          deltaZ < -config.cotThetaMax * deltaR) {
        continue;
      }

      const float iDeltaR = std::sqrt(iDeltaR2);
      const float cotTheta = deltaZ * iDeltaR;

      // discard bottom-middle doublets in a certain (r, eta) region according
      // to detector specific cuts
      if constexpr (isBottomCandidate) {
        if (!config.experimentCuts(otherSp.extra(rColumn), cotTheta)) {
          continue;
        }
      }

      const float er = calculateError(otherSp, iDeltaR2, cotTheta);

      // fill output vectors
      compatibleSp.emplace_back(otherSp.index());
      linCircles.emplace_back(cotTheta, iDeltaR, er, uT, vT, xNewFrame,
                              yNewFrame);
      continue;
    }

    // in the rotated frame the interaction point is positioned at x = -rM
    // and y ~= impactParam
    const float vIP = (yNewFrame > 0) ? -vIPAbs : vIPAbs;

    // we can obtain aCoef as the slope dv/du of the linear function,
    // estimated using du and dv between the two SP bCoef is obtained by
    // inserting aCoef into the linear equation
    const float aCoef = (vT - vIP) / (uT - cuts.uIP);
    const float bCoef = vIP - aCoef * cuts.uIP;
    // the distance of the straight line from the origin (radius of the
    // circle) is related to aCoef and bCoef by d^2 = bCoef^2 / (1 +
    // aCoef^2) = 1 / (radius^2) and we can apply the cut on the curvature
    if ((bCoef * bCoef) * options.minHelixDiameter2 > 1 + aCoef * aCoef) {
      continue;
    }

    // check if duplet cotTheta is within the region of interest
    // cotTheta is defined as (deltaZ / deltaR) but instead we multiply
    // cotThetaMax by deltaR to avoid division
    //
    // intentionally using `|` after profiling. faster due to better branch
    // prediction
    if (deltaZ > config.cotThetaMax * deltaR |
        deltaZ < -config.cotThetaMax * deltaR) {
      continue;
    }

    const float iDeltaR = std::sqrt(iDeltaR2);
    const float cotTheta = deltaZ * iDeltaR;

    // discard bottom-middle doublets in a certain (r, eta) region according
    // to detector specific cuts
    if constexpr (isBottomCandidate) {
      if (!config.experimentCuts(otherSp.extra(rColumn), cotTheta)) {
        continue;
      }
    }

    const float er = calculateError(otherSp, iDeltaR2, cotTheta);

    // fill output vectors
    compatibleSp.emplace_back(otherSp.index());
    linCircles.emplace_back(cotTheta, iDeltaR, er, uT, vT, xNewFrame,
                            yNewFrame);
  }
}

template <TripletSeedFinder2::MeasurementInfo measurement_info>
void TripletSeedFinder2::createTriplets(
    const DerivedConfig& config, const DerivedOptions& options,
    const TripletSeedFilter2::Options& filterOptions,
    TripletSeedFilter2::State& filterState,
    const SpacePointContainer2& spacePoints,
    const SpacePointContainer2::DenseColumn<float>& rColumn,
    const SpacePointContainer2::DenseColumn<float>* varianceRColumn,
    const SpacePointContainer2::DenseColumn<float>* varianceZColumn,
    const SpacePointContainer2::DenseColumn<Vector3>* topStripVectorColumn,
    const SpacePointContainer2::DenseColumn<Vector3>* bottomStripVectorColumn,
    const SpacePointContainer2::DenseColumn<Vector3>* stripCenterDistanceColumn,
    const SpacePointContainer2::DenseColumn<Vector3>*
        topStripCenterPositionColumn,
    const ConstSpacePointProxy2& spM,
    const std::vector<SpacePointIndex2>& bottomSps,
    const std::vector<LinCircle>& bottomLinCircles,
    const std::vector<float>& bottomCotThetas,
    const std::vector<SpacePointIndex2>& topSps,
    const std::vector<LinCircle>& topLinCircles,
    const std::vector<float>& topCotThetas,
    std::vector<std::uint32_t>& sortedBottoms,
    std::vector<std::uint32_t>& sortedTops,
    std::vector<SpacePointIndex2>& topSpVec, std::vector<float>& curvatures,
    std::vector<float>& impactParameters,
    CandidatesForMiddleSp2& candidatesCollector) {
  constexpr bool isDetailedMeasurement =
      measurement_info == MeasurementInfo::eDetailed;

  const float rM = spM.extra(rColumn);
  const float cosPhiM = spM.x() / rM;
  const float sinPhiM = spM.y() / rM;
  const float varianceRM =
      varianceRColumn != nullptr ? spM.extra(*varianceRColumn) : 0;
  const float varianceZM =
      varianceZColumn != nullptr ? spM.extra(*varianceZColumn) : 0;

  // make index vectors for sorting
  sortedBottoms.resize(bottomSps.size());
  std::iota(sortedBottoms.begin(), sortedBottoms.end(), 0);

  sortedTops.resize(topSps.size());
  std::iota(sortedTops.begin(), sortedTops.end(), 0);

  if constexpr (!isDetailedMeasurement) {
    // sorting becomes less expensive when we copy the cotTheta values into
    // their own arrays due to more optimal usage of the cache

    std::ranges::sort(
        sortedBottoms, {},
        [&bottomCotThetas](const std::size_t s) { return bottomCotThetas[s]; });
    std::ranges::sort(sortedTops, {}, [&topCotThetas](const std::size_t s) {
      return topCotThetas[s];
    });
  }

  // Reserve enough space, in case current capacity is too little
  topSpVec.reserve(topSps.size());
  curvatures.reserve(topSps.size());
  impactParameters.reserve(topSps.size());

  std::size_t t0 = 0;

  for (const std::size_t b : sortedBottoms) {
    // break if we reached the last top SP
    if (t0 >= topSps.size()) {
      break;
    }

    auto spB = spacePoints.at(bottomSps[b]);

    const auto& lb = bottomLinCircles[b];
    float cotThetaB = lb.cotTheta;
    float Vb = lb.V;
    float Ub = lb.U;
    float ErB = lb.Er;
    float iDeltaRB = lb.iDeltaR;

    // 1+(cot^2(theta)) = 1/sin^2(theta)
    float iSinTheta2 = 1. + cotThetaB * cotThetaB;
    float sigmaSquaredPtDependent = iSinTheta2 * options.sigmapT2perRadius;
    // calculate max scattering for min momentum at the seed's theta angle
    // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
    // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
    // scattering
    // but to avoid trig functions we approximate cot by scaling by
    // 1/sin^4(theta)
    // resolving with pT to p scaling --> only divide by sin^2(theta)
    // max approximation error for allowed scattering angles of 0.04 rad at
    // eta=infinity: ~8.5%
    float scatteringInRegion2 = options.multipleScattering2 * iSinTheta2;

    float sinTheta = 1 / std::sqrt(iSinTheta2);
    float cosTheta = cotThetaB * sinTheta;

    // clear all vectors used in each inner for loop
    topSpVec.clear();
    curvatures.clear();
    impactParameters.clear();

    // coordinate transformation and checks for middle spacepoint
    // x and y terms for the rotation from UV to XY plane
    float rotationTermsUVtoXY[2] = {0, 0};
    if constexpr (isDetailedMeasurement) {
      rotationTermsUVtoXY[0] = cosPhiM * sinTheta;
      rotationTermsUVtoXY[1] = sinPhiM * sinTheta;
    }

    // minimum number of compatible top SPs to trigger the filter for a certain
    // middle bottom pair if seedConfirmation is false we always ask for at
    // least one compatible top to trigger the filter
    std::size_t minCompatibleTopSPs = 2;
    if (!config.seedConfirmation ||
        spB.extra(rColumn) > filterOptions.seedConfRange.rMaxSeedConf) {
      minCompatibleTopSPs = 1;
    }
    if (config.seedConfirmation &&
        candidatesCollector.nHighQualityCandidates() > 0) {
      minCompatibleTopSPs++;
    }

    for (std::size_t indexSortedTop = t0; indexSortedTop < topSps.size();
         ++indexSortedTop) {
      const std::size_t t = sortedTops[indexSortedTop];

      auto spT = spacePoints.at(topSps[t]);

      const auto& lt = topLinCircles[t];

      float cotThetaT = lt.cotTheta;
      float rMxy = 0.;
      float ub = 0.;
      float vb = 0.;
      float ut = 0.;
      float vt = 0.;
      double rMTransf[3];
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
        double positionMiddle[3] = {
            rotationTermsUVtoXY[0] - rotationTermsUVtoXY[1] * A0,
            rotationTermsUVtoXY[0] * A0 + rotationTermsUVtoXY[1],
            zPositionMiddle};

        if (!xyzCoordinateCheck(
                config.toleranceParam, spM, *topStripVectorColumn,
                *bottomStripVectorColumn, *stripCenterDistanceColumn,
                *topStripCenterPositionColumn, positionMiddle, rMTransf)) {
          continue;
        }

        // coordinate transformation and checks for bottom spacepoint
        float B0 = 2. * (Vb - A0 * Ub);
        float Cb = 1. - B0 * lb.y;
        float Sb = A0 + B0 * lb.x;
        double positionBottom[3] = {
            rotationTermsUVtoXY[0] * Cb - rotationTermsUVtoXY[1] * Sb,
            rotationTermsUVtoXY[0] * Sb + rotationTermsUVtoXY[1] * Cb,
            zPositionMiddle};

        double rBTransf[3];
        if (!xyzCoordinateCheck(
                config.toleranceParam, spB, *topStripVectorColumn,
                *bottomStripVectorColumn, *stripCenterDistanceColumn,
                *topStripCenterPositionColumn, positionBottom, rBTransf)) {
          continue;
        }

        // coordinate transformation and checks for top spacepoint
        float Ct = 1. - B0 * lt.y;
        float St = A0 + B0 * lt.x;
        double positionTop[3] = {
            rotationTermsUVtoXY[0] * Ct - rotationTermsUVtoXY[1] * St,
            rotationTermsUVtoXY[0] * St + rotationTermsUVtoXY[1] * Ct,
            zPositionMiddle};

        double rTTransf[3];
        if (!xyzCoordinateCheck(
                config.toleranceParam, spT, *topStripVectorColumn,
                *bottomStripVectorColumn, *stripCenterDistanceColumn,
                *topStripCenterPositionColumn, positionTop, rTTransf)) {
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
      if (S2 < B2 * options.minHelixDiameter2) {
        continue;
      }

      // refinement of the cut on the compatibility between the r-z slope of
      // the two seed segments using a scattering term scaled by the actual
      // measured pT (p2scatterSigma)
      float iHelixDiameter2 = B2 / S2;
      // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
      // from rad to deltaCotTheta
      float p2scatterSigma = iHelixDiameter2 * sigmaSquaredPtDependent;
      if (!std::isinf(config.maxPtScattering)) {
        // if pT > maxPtScattering, calculate allowed scattering angle using
        // maxPtScattering instead of pt.
        // To avoid 0-divison the pT check is skipped in case of B2==0, and
        // p2scatterSigma is calculated directly from maxPtScattering
        if (B2 == 0 || options.pTPerHelixRadius * std::sqrt(S2 / B2) >
                           2. * config.maxPtScattering) {
          float pTscatterSigma = (config.highland / config.maxPtScattering) *
                                 config.sigmaScattering;
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

      if (Im > config.impactMax) {
        continue;
      }

      topSpVec.push_back(spT.index());
      // inverse diameter is signed depending on if the curvature is
      // positive/negative in phi
      curvatures.push_back(B / std::sqrt(S2));
      impactParameters.push_back(Im);
    }  // loop on tops

    // continue if number of top SPs is smaller than minimum required for filter
    if (topSpVec.size() < minCompatibleTopSPs) {
      continue;
    }

    float zOrigin = spM.z() - rM * lb.cotTheta;
    config.seedFilter->filter2SpFixed(filterOptions, filterState, spacePoints,
                                      rColumn, spB.index(), spM.index(),
                                      topSpVec, curvatures, impactParameters,
                                      zOrigin, candidatesCollector);
  }  // loop on bottoms
}

bool TripletSeedFinder2::xyzCoordinateCheck(
    double toleranceParam, const ConstSpacePointProxy2& sp,
    const SpacePointContainer2::DenseColumn<Vector3>& topStripVectorColumn,
    const SpacePointContainer2::DenseColumn<Vector3>& bottomStripVectorColumn,
    const SpacePointContainer2::DenseColumn<Vector3>& stripCenterDistanceColumn,
    const SpacePointContainer2::DenseColumn<Vector3>&
        topStripCenterPositionColumn,
    const double* spacepointPosition, double* outputCoordinates) {
  const Vector3& topStripVector = sp.extra(topStripVectorColumn);
  const Vector3& bottomStripVector = sp.extra(bottomStripVectorColumn);
  const Vector3& stripCenterDistance = sp.extra(stripCenterDistanceColumn);

  const double xTopStripVector = topStripVector[0];
  const double yTopStripVector = topStripVector[1];
  const double zTopStripVector = topStripVector[2];
  const double xBottomStripVector = bottomStripVector[0];
  const double yBottomStripVector = bottomStripVector[1];
  const double zBottomStripVector = bottomStripVector[2];

  // cross product between top strip vector and spacepointPosition
  double d1[3] = {yTopStripVector * spacepointPosition[2] -
                      zTopStripVector * spacepointPosition[1],
                  zTopStripVector * spacepointPosition[0] -
                      xTopStripVector * spacepointPosition[2],
                  xTopStripVector * spacepointPosition[1] -
                      yTopStripVector * spacepointPosition[0]};

  // scalar product between bottom strip vector and d1
  double bd1 = xBottomStripVector * d1[0] + yBottomStripVector * d1[1] +
               zBottomStripVector * d1[2];

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the bottom detector element
  double s1 = stripCenterDistance[0] * d1[0] + stripCenterDistance[1] * d1[1] +
              stripCenterDistance[2] * d1[2];
  if (std::abs(s1) > std::abs(bd1) * toleranceParam) {
    return false;
  }

  // cross product between bottom strip vector and spacepointPosition
  double d0[3] = {yBottomStripVector * spacepointPosition[2] -
                      zBottomStripVector * spacepointPosition[1],
                  zBottomStripVector * spacepointPosition[0] -
                      xBottomStripVector * spacepointPosition[2],
                  xBottomStripVector * spacepointPosition[1] -
                      yBottomStripVector * spacepointPosition[0]};

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the top detector element
  double s0 = stripCenterDistance[0] * d0[0] + stripCenterDistance[1] * d0[1] +
              stripCenterDistance[2] * d0[2];
  if (std::abs(s0) > std::abs(bd1) * toleranceParam) {
    return false;
  }

  // if arrive here spacepointPosition is compatible with strip directions and
  // detector elements

  const Vector3& topStripCenterPosition =
      sp.extra(topStripCenterPositionColumn);

  // spacepointPosition corrected with respect to the top strip position and
  // direction and the distance between the strips
  s0 = s0 / bd1;
  outputCoordinates[0] = topStripCenterPosition[0] + xTopStripVector * s0;
  outputCoordinates[1] = topStripCenterPosition[1] + yTopStripVector * s0;
  outputCoordinates[2] = topStripCenterPosition[2] + zTopStripVector * s0;
  return true;
}

}  // namespace Acts
