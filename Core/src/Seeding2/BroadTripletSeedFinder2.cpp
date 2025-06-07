// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/BroadTripletSeedFinder2.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding2/BroadTripletSeedFilter2.hpp"
#include "Acts/Seeding2/DoubletSeedFinder2.hpp"

#include <numeric>

#include <Eigen/src/Core/Matrix.h>

namespace Acts {

using namespace UnitLiterals;

BroadTripletSeedFinder2::DerivedTripletCuts
BroadTripletSeedFinder2::TripletCuts::derive(float bFieldInZ) const {
  DerivedTripletCuts result;

  static_cast<TripletCuts&>(result) = *this;

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
  const float maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;

  // bFieldInZ is in (pT/radius) natively, no need for conversion
  result.pTPerHelixRadius = bFieldInZ;
  result.minHelixDiameter2 =
      std::pow(result.minPt * 2 / result.pTPerHelixRadius, 2) *
      result.helixCutTolerance;
  const float pT2perRadius =
      std::pow(result.highland / result.pTPerHelixRadius, 2);
  result.sigmapT2perRadius =
      pT2perRadius * std::pow(2 * result.sigmaScattering, 2);
  result.multipleScattering2 =
      maxScatteringAngle2 * std::pow(result.sigmaScattering, 2);

  return result;
}

BroadTripletSeedFinder2::BroadTripletSeedFinder2(
    const Config& config, std::unique_ptr<const Logger> logger_)
    : m_cfg(config), m_logger(std::move(logger_)) {}

void BroadTripletSeedFinder2::createSeedsFromGroup(
    const Options& options, State& state, Cache& cache,
    const DoubletSeedFinder2::DerivedCuts& bottomCuts,
    const DoubletSeedFinder2::DerivedCuts& topCuts,
    const DerivedTripletCuts& tripletCuts,
    const BroadTripletSeedFilter2& filter,
    const SpacePointContainerPointers2& containerPointers,
    std::span<const SpacePointIndex2> bottomSps, SpacePointIndex2 middleSp,
    std::span<const SpacePointIndex2> topSps,
    SeedContainer2& outputSeeds) const {
  cache.candidatesCollector.setMaxElements(
      options.filter.maxSeedsPerSpMConf,
      options.filter.maxQualitySeedsPerSpMConf);

  auto spM = containerPointers.spacePoints().at(middleSp);

  DoubletSeedFinder2::MiddleSpInfo middleSpInfo =
      DoubletSeedFinder2::computeMiddleSpInfo(spM, containerPointers.rColumn());

  // create middle-top doublets
  cache.topDoublets.clear();
  DoubletSeedFinder2::createDoublets<DoubletSeedFinder2::eTop>(
      topCuts, containerPointers, spM, middleSpInfo, topSps, cache.topDoublets);

  // no top SP found -> try next spM
  if (cache.topDoublets.empty()) {
    ACTS_VERBOSE("No compatible Tops, returning");
    return;
  }

  // create middle-bottom doublets
  cache.bottomDoublets.clear();
  DoubletSeedFinder2::createDoublets<DoubletSeedFinder2::eBottom>(
      bottomCuts, containerPointers, spM, middleSpInfo, bottomSps,
      cache.bottomDoublets);

  // no bottom SP found -> try next spM
  if (cache.bottomDoublets.empty()) {
    ACTS_VERBOSE("No compatible Bottoms, returning");
    return;
  }

  ACTS_VERBOSE("Candidates: " << cache.bottomDoublets.size() << " bottoms and "
                              << cache.topDoublets.size()
                              << " tops for middle candidate indexed "
                              << spM.index());

  // combine doublets to triplets
  cache.candidatesCollector.clear();
  if (options.useDetailedDoubleMeasurementInfo) {
    createTripletsDetailed(
        tripletCuts, filter, options.filter, state.filter, cache.filter,
        containerPointers, spM, cache.bottomDoublets, cache.topDoublets,
        cache.tripletTopCandidates, cache.candidatesCollector);
  } else {
    createTriplets(cache.tripletCache, tripletCuts, filter, options.filter,
                   state.filter, cache.filter, containerPointers, spM,
                   cache.bottomDoublets, cache.topDoublets,
                   cache.tripletTopCandidates, cache.candidatesCollector);
  }

  // retrieve all candidates
  // this collection is already sorted, higher weights first
  const std::size_t numQualitySeeds =
      cache.candidatesCollector.nHighQualityCandidates();
  cache.candidatesCollector.toSortedCandidates(containerPointers.spacePoints(),
                                               cache.sortedCandidates);
  filter.filter1SpFixed(options.filter, state.filter,
                        containerPointers.spacePoints(),
                        containerPointers.copiedFromIndexColumn,
                        cache.sortedCandidates, numQualitySeeds, outputSeeds);
}

void BroadTripletSeedFinder2::createTriplets(
    TripletCache& cache, const DerivedTripletCuts& cuts,
    const BroadTripletSeedFilter2& filter,
    const BroadTripletSeedFilter2::Options& filterOptions,
    BroadTripletSeedFilter2::State& filterState,
    BroadTripletSeedFilter2::Cache& filterCache,
    const SpacePointContainerPointers2& containerPointers,
    const ConstSpacePointProxy2& spM,
    const DoubletSeedFinder2::DoubletsForMiddleSp& bottomDoublets,
    const DoubletSeedFinder2::DoubletsForMiddleSp& topDoublets,
    TripletTopCandidates& tripletTopCandidates,
    CandidatesForMiddleSp2& candidatesCollector) {
  const float rM = spM.extra(containerPointers.rColumn());
  // TOD use some reasonable defaults
  const float varianceRM = containerPointers.hasVarianceColumns()
                               ? spM.extra(containerPointers.varianceRColumn())
                               : 0;
  const float varianceZM = containerPointers.hasVarianceColumns()
                               ? spM.extra(containerPointers.varianceZColumn())
                               : 0;

  // make index vectors for sorting
  cache.sortedBottoms.resize(bottomDoublets.size());
  std::iota(cache.sortedBottoms.begin(), cache.sortedBottoms.end(), 0);
  std::ranges::sort(cache.sortedBottoms, {},
                    [&bottomDoublets](const std::size_t s) {
                      return bottomDoublets.cotTheta[s];
                    });

  cache.sortedTops.resize(topDoublets.size());
  std::iota(cache.sortedTops.begin(), cache.sortedTops.end(), 0);
  std::ranges::sort(cache.sortedTops, {}, [&topDoublets](const std::size_t s) {
    return topDoublets.cotTheta[s];
  });

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

    // clear all vectors used in each inner for loop
    tripletTopCandidates.clear();

    // minimum number of compatible top SPs to trigger the filter for a certain
    // middle bottom pair if seedConfirmation is false we always ask for at
    // least one compatible top to trigger the filter
    std::size_t minCompatibleTopSPs = 2;
    if (!filterOptions.seedConfirmation ||
        spB.extra(containerPointers.rColumn()) >
            filterOptions.seedConfRange.rMaxSeedConf) {
      minCompatibleTopSPs = 1;
    }
    if (filterOptions.seedConfirmation &&
        candidatesCollector.nHighQualityCandidates() > 0) {
      minCompatibleTopSPs++;
    }

    for (std::size_t indexSortedTop = t0; indexSortedTop < topDoublets.size();
         ++indexSortedTop) {
      const std::size_t t = cache.sortedTops[indexSortedTop];
      auto spT = containerPointers.spacePoints().at(topDoublets.spacePoints[t]);
      const auto& lt = topDoublets.linCircles[t];
      float cotThetaT = lt.cotTheta;

      // use geometric average
      float cotThetaAvg2 = cotThetaB * cotThetaT;
      if (cotThetaAvg2 <= 0) {
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
        // break if cotTheta from bottom SP < cotTheta from top SP because
        // the SP are sorted by cotTheta
        if (cotThetaB < cotThetaT) {
          break;
        }
        t0 = indexSortedTop + 1;
        continue;
      }

      float dU = lt.U - Ub;
      // protects against division by 0
      if (dU == 0.) {
        continue;
      }
      // A and B are evaluated as a function of the circumference parameters
      // x_0 and y_0
      float A = (lt.V - Vb) / dU;
      float S2 = 1. + A * A;
      float B = Vb - A * Ub;
      float B2 = B * B;

      // sqrt(S2)/B = 2 * helixradius
      // calculated radius must not be smaller than minimum radius
      if (S2 < B2 * cuts.minHelixDiameter2) {
        continue;
      }

      // refinement of the cut on the compatibility between the r-z slope of
      // the two seed segments using a scattering term scaled by the actual
      // measured pT (p2scatterSigma)
      float iHelixDiameter2 = B2 / S2;
      float sigmaSquaredPtDependent = iSinTheta2 * cuts.sigmapT2perRadius;
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
        if (cotThetaB < cotThetaT) {
          break;
        }
        t0 = indexSortedTop + 1;
        continue;
      }

      // A and B allow calculation of impact params in U/V plane with linear
      // function
      // (in contrast to having to solve a quadratic function in x/y plane)
      float Im = std::abs((A - B * rM) * rM);
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

void BroadTripletSeedFinder2::createTripletsDetailed(
    const DerivedTripletCuts& cuts, const BroadTripletSeedFilter2& filter,
    const BroadTripletSeedFilter2::Options& filterOptions,
    BroadTripletSeedFilter2::State& filterState,
    BroadTripletSeedFilter2::Cache& filterCache,
    const SpacePointContainerPointers2& containerPointers,
    const ConstSpacePointProxy2& spM,
    const DoubletSeedFinder2::DoubletsForMiddleSp& bottomDoublets,
    const DoubletSeedFinder2::DoubletsForMiddleSp& topDoublets,
    TripletTopCandidates& tripletTopCandidates,
    CandidatesForMiddleSp2& candidatesCollector) {
  const float rM = spM.extra(containerPointers.rColumn());
  const float cosPhiM = spM.x() / rM;
  const float sinPhiM = spM.y() / rM;
  // TOD use some reasonable defaults
  const float varianceRM = containerPointers.hasVarianceColumns()
                               ? spM.extra(containerPointers.varianceRColumn())
                               : 0;
  const float varianceZM = containerPointers.hasVarianceColumns()
                               ? spM.extra(containerPointers.varianceZColumn())
                               : 0;

  // Reserve enough space, in case current capacity is too little
  tripletTopCandidates.resize(topDoublets.size());

  for (std::size_t b = 0; b < bottomDoublets.size(); ++b) {
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
    Eigen::Vector2f rotationTermsUVtoXY = {cosPhiM * sinTheta,
                                           sinPhiM * sinTheta};

    // minimum number of compatible top SPs to trigger the filter for a certain
    // middle bottom pair if seedConfirmation is false we always ask for at
    // least one compatible top to trigger the filter
    std::size_t minCompatibleTopSPs = 2;
    if (!filterOptions.seedConfirmation ||
        spB.extra(containerPointers.rColumn()) >
            filterOptions.seedConfRange.rMaxSeedConf) {
      minCompatibleTopSPs = 1;
    }
    if (filterOptions.seedConfirmation &&
        candidatesCollector.nHighQualityCandidates() > 0) {
      minCompatibleTopSPs++;
    }

    for (std::size_t t = 0; t < topDoublets.size(); ++t) {
      auto spT = containerPointers.spacePoints().at(topDoublets.spacePoints[t]);
      const auto& lt = topDoublets.linCircles[t];

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

      Vector3 rMTransf;
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
      float xB = rBTransf[0] - rMTransf[0];
      float yB = rBTransf[1] - rMTransf[1];
      float zB = rBTransf[2] - rMTransf[2];
      float xT = rTTransf[0] - rMTransf[0];
      float yT = rTTransf[1] - rMTransf[1];
      float zT = rTTransf[2] - rMTransf[2];

      float iDeltaRB2 = 1. / (xB * xB + yB * yB);
      float iDeltaRT2 = 1. / (xT * xT + yT * yT);

      cotThetaB = -zB * std::sqrt(iDeltaRB2);
      float cotThetaT = zT * std::sqrt(iDeltaRT2);

      // use arithmetic average
      float averageCotTheta = 0.5 * (cotThetaB + cotThetaT);
      float cotThetaAvg2 = averageCotTheta * averageCotTheta;

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
        continue;
      }

      float rMxy =
          std::sqrt(rMTransf[0] * rMTransf[0] + rMTransf[1] * rMTransf[1]);
      float irMxy = 1 / rMxy;
      float Ax = rMTransf[0] * irMxy;
      float Ay = rMTransf[1] * irMxy;

      float ub = (xB * Ax + yB * Ay) * iDeltaRB2;
      float vb = (yB * Ax - xB * Ay) * iDeltaRB2;
      float ut = (xT * Ax + yT * Ay) * iDeltaRT2;
      float vt = (yT * Ax - xT * Ay) * iDeltaRT2;

      dU = ut - ub;
      // protects against division by 0
      if (dU == 0.) {
        continue;
      }
      float A = (vt - vb) / dU;
      float S2 = 1. + A * A;
      float B = vb - A * ub;
      float B2 = B * B;

      // sqrt(S2)/B = 2 * helixradius
      // calculated radius must not be smaller than minimum radius
      if (S2 < B2 * cuts.minHelixDiameter2) {
        continue;
      }

      // refinement of the cut on the compatibility between the r-z slope of
      // the two seed segments using a scattering term scaled by the actual
      // measured pT (p2scatterSigma)
      float iHelixDiameter2 = B2 / S2;
      float sigmaSquaredPtDependent = iSinTheta2 * cuts.sigmapT2perRadius;
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
        continue;
      }

      // A and B allow calculation of impact params in U/V plane with linear
      // function
      // (in contrast to having to solve a quadratic function in x/y plane)
      float Im = std::abs((A - B * rMxy) * rMxy);
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

bool BroadTripletSeedFinder2::stripCoordinateCheck(
    float tolerance, const ConstSpacePointProxy2& sp,
    const SpacePointContainerPointers2& containerPointers,
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
  float bd1 = bottomStripVector.dot(d1);

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the bottom detector element
  float s1 = stripCenterDistance.dot(d1);
  if (std::abs(s1) > std::abs(bd1) * tolerance) {
    return false;
  }

  // cross product between bottom strip vector and spacepointPosition
  Vector3 d0 = bottomStripVector.cross(spacePointPosition);

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the top detector element
  float s0 = stripCenterDistance.dot(d0);
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
