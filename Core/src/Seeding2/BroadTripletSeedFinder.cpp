// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/BroadTripletSeedFinder.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding2/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding2/DoubletSeedFinder.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <numeric>

#include <Eigen/src/Core/Matrix.h>

using namespace Acts::UnitLiterals;

namespace Acts::Experimental {

namespace {

/// Check the compatibility of strip space point coordinates in xyz assuming
/// the Bottom-Middle direction with the strip measurement details
bool stripCoordinateCheck(float tolerance, const ConstSpacePointProxy2& sp,
                          const Eigen::Vector3f& spacePointPosition,
                          Eigen::Vector3f& outputCoordinates) {
  const Eigen::Vector3f& topStripVector = sp.topStripVector();
  const Eigen::Vector3f& bottomStripVector = sp.bottomStripVector();
  const Eigen::Vector3f& stripCenterDistance = sp.stripCenterDistance();

  // cross product between top strip vector and spacepointPosition
  Eigen::Vector3f d1 = topStripVector.cross(spacePointPosition);

  // scalar product between bottom strip vector and d1
  float bd1 = bottomStripVector.dot(d1);

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the bottom detector element
  float s1 = stripCenterDistance.dot(d1);
  if (std::abs(s1) > std::abs(bd1) * tolerance) {
    return false;
  }

  // cross product between bottom strip vector and spacepointPosition
  Eigen::Vector3f d0 = bottomStripVector.cross(spacePointPosition);

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the top detector element
  float s0 = stripCenterDistance.dot(d0);
  if (std::abs(s0) > std::abs(bd1) * tolerance) {
    return false;
  }

  // if arrive here spacepointPosition is compatible with strip directions and
  // detector elements

  const Eigen::Vector3f& topStripCenter = sp.topStripCenter();

  // spacepointPosition corrected with respect to the top strip position and
  // direction and the distance between the strips
  s0 = s0 / bd1;
  outputCoordinates = topStripCenter + topStripVector * s0;
  return true;
}

}  // namespace

BroadTripletSeedFinder::DerivedTripletCuts::DerivedTripletCuts(
    const TripletCuts& cuts, float bFieldInZ_)
    : TripletCuts(cuts), bFieldInZ(bFieldInZ_) {
  // similar to `theta0Highland` in `Core/src/Material/Interactions.cpp`
  {
    const double xOverX0 = radLengthPerSeed;
    const double q2OverBeta2 = 1;  // q^2=1, beta^2~1
    // RPP2018 eq. 33.15 (treats beta and q² consistently)
    const double t = std::sqrt(xOverX0 * q2OverBeta2);
    // log((x/X0) * (q²/beta²)) = log((sqrt(x/X0) * (q/beta))²)
    //                          = 2 * log(sqrt(x/X0) * (q/beta))
    highland =
        static_cast<float>(13.6_MeV * t * (1.0 + 0.038 * 2 * std::log(t)));
  }

  const float maxScatteringAngle = highland / minPt;
  const float maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;

  // bFieldInZ is in (pT/radius) natively, no need for conversion
  pTPerHelixRadius = bFieldInZ;
  minHelixDiameter2 = square(minPt * 2 / pTPerHelixRadius) * helixCutTolerance;
  const float pT2perRadius = square(highland / pTPerHelixRadius);
  sigmapT2perRadius = pT2perRadius * square(2 * sigmaScattering);
  multipleScattering2 = maxScatteringAngle2 * square(sigmaScattering);
}

BroadTripletSeedFinder::BroadTripletSeedFinder(
    std::unique_ptr<const Logger> logger_)
    : m_logger(std::move(logger_)) {}

void BroadTripletSeedFinder::createSeedsFromGroup(
    const Options& options, State& state, Cache& cache,
    const DoubletSeedFinder& bottomFinder, const DoubletSeedFinder& topFinder,
    const DerivedTripletCuts& tripletCuts, const BroadTripletSeedFilter& filter,
    const SpacePointContainer2& spacePoints,
    std::span<const SpacePointIndex2> bottomSps, SpacePointIndex2 middleSp,
    std::span<const SpacePointIndex2> topSps,
    SeedContainer2& outputSeeds) const {
  cache.candidatesCollector.setMaxElements(
      filter.config().maxSeedsPerSpMConf,
      filter.config().maxQualitySeedsPerSpMConf);

  auto spM = spacePoints[middleSp];

  DoubletSeedFinder::MiddleSpInfo middleSpInfo =
      DoubletSeedFinder::computeMiddleSpInfo(spM);

  // create middle-top doublets
  cache.topDoublets.clear();
  topFinder.createDoublets(spacePoints, spM, middleSpInfo, topSps,
                           cache.topDoublets);

  // no top SP found -> cannot form any triplet
  if (cache.topDoublets.empty()) {
    ACTS_VERBOSE("No compatible Tops, returning");
    return;
  }

  // apply cut on the number of top SP if seedConfirmation is true
  float rMaxSeedConf = 0;
  if (filter.config().seedConfirmation) {
    // check if middle SP is in the central or forward region
    const bool isForwardRegion =
        spM.z() > filter.config().centralSeedConfirmationRange.zMaxSeedConf ||
        spM.z() < filter.config().centralSeedConfirmationRange.zMinSeedConf;
    SeedConfirmationRangeConfig seedConfRange =
        isForwardRegion ? filter.config().forwardSeedConfirmationRange
                        : filter.config().centralSeedConfirmationRange;
    // set the minimum number of top SP depending on whether the middle SP is
    // in the central or forward region
    std::size_t nTopSeedConf = spM.r() > seedConfRange.rMaxSeedConf
                                   ? seedConfRange.nTopForLargeR
                                   : seedConfRange.nTopForSmallR;
    // set max bottom radius for seed confirmation
    rMaxSeedConf = seedConfRange.rMaxSeedConf;
    // continue if number of top SPs is smaller than minimum
    if (cache.topDoublets.size() < nTopSeedConf) {
      ACTS_VERBOSE("Number of top SPs is "
                   << cache.topDoublets.size()
                   << " and is smaller than minimum, returning");
      return;
    }
  }

  // create middle-bottom doublets
  cache.bottomDoublets.clear();
  bottomFinder.createDoublets(spacePoints, spM, middleSpInfo, bottomSps,
                              cache.bottomDoublets);

  // no bottom SP found -> cannot form any triplet
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
  if (options.useStripMeasurementInfo) {
    createStripTriplets(tripletCuts, rMaxSeedConf, filter, state.filter,
                        cache.filter, spacePoints, spM, cache.bottomDoublets,
                        cache.topDoublets, cache.tripletTopCandidates,
                        cache.candidatesCollector);
  } else {
    createTriplets(cache.tripletCache, tripletCuts, rMaxSeedConf, filter,
                   state.filter, cache.filter, spacePoints, spM,
                   cache.bottomDoublets, cache.topDoublets,
                   cache.tripletTopCandidates, cache.candidatesCollector);
  }

  // retrieve all candidates
  // this collection is already sorted, higher weights first
  const std::size_t numQualitySeeds =
      cache.candidatesCollector.nHighQualityCandidates();
  cache.candidatesCollector.toSortedCandidates(spacePoints,
                                               cache.sortedCandidates);
  filter.filter1SpFixed(state.filter, spacePoints, cache.sortedCandidates,
                        numQualitySeeds, outputSeeds);
}

void BroadTripletSeedFinder::createSeedsFromSortedGroups(
    const Options& options, State& state, Cache& cache,
    const DoubletSeedFinder& bottomFinder, const DoubletSeedFinder& topFinder,
    const DerivedTripletCuts& tripletCuts, const BroadTripletSeedFilter& filter,
    const SpacePointContainer2& spacePoints,
    const std::vector<std::span<const SpacePointIndex2>>& bottomSpGroups,
    std::span<const SpacePointIndex2> middleSps,
    const std::vector<std::span<const SpacePointIndex2>>& topSpGroups,
    const std::pair<float, float>& radiusRangeForMiddle,
    SeedContainer2& outputSeeds) const {
  if (middleSps.empty()) {
    return;
  }

  // initialize cache
  cache.candidatesCollector.setMaxElements(
      filter.config().maxSeedsPerSpMConf,
      filter.config().maxQualitySeedsPerSpMConf);

  cache.bottomSpOffsets.clear();
  cache.topSpOffsets.clear();

  // Initialize initial offsets for bottom and top space points with binary
  // search. This requires at least one middle space point to be present which
  // is already checked above.
  auto firstMiddleSp = spacePoints[middleSps.front()];
  float firstMiddleSpR = firstMiddleSp.r();

  std::ranges::transform(
      bottomSpGroups, std::back_inserter(cache.bottomSpOffsets),
      [&](const std::span<const SpacePointIndex2>& bottomSps) {
        auto low = std::ranges::lower_bound(
            bottomSps, firstMiddleSpR - bottomFinder.config().deltaRMax, {},
            [&](const SpacePointIndex2& spIndex) {
              auto sp = spacePoints[spIndex];
              return sp.r();
            });
        return low - bottomSps.begin();
      });
  std::ranges::transform(topSpGroups, std::back_inserter(cache.topSpOffsets),
                         [&](const std::span<const SpacePointIndex2>& topSps) {
                           auto low = std::ranges::lower_bound(
                               topSps,
                               firstMiddleSpR + topFinder.config().deltaRMin,
                               {}, [&](const SpacePointIndex2& spIndex) {
                                 auto sp = spacePoints[spIndex];
                                 return sp.r();
                               });
                           return low - topSps.begin();
                         });

  for (SpacePointIndex2 middleSp : middleSps) {
    auto spM = spacePoints[middleSp];
    const float rM = spM.r();

    // check if spM is outside our radial region of interest
    if (rM < radiusRangeForMiddle.first) {
      continue;
    }
    if (rM > radiusRangeForMiddle.second) {
      // break because SPs are sorted in r
      break;
    }

    DoubletSeedFinder::MiddleSpInfo middleSpInfo =
        DoubletSeedFinder::computeMiddleSpInfo(spM);

    // create middle-top doublets
    cache.topDoublets.clear();
    for (std::size_t i = 0; i < topSpGroups.size(); ++i) {
      topFinder.createSortedDoublets(spacePoints, spM, middleSpInfo,
                                     topSpGroups[i], cache.topSpOffsets[i],
                                     cache.topDoublets);
    }

    // no top SP found -> try next spM
    if (cache.topDoublets.empty()) {
      ACTS_VERBOSE("No compatible Tops, moving to next middle candidate");
      continue;
    }

    // apply cut on the number of top SP if seedConfirmation is true
    float rMaxSeedConf = 0;
    if (filter.config().seedConfirmation) {
      // check if middle SP is in the central or forward region
      const bool isForwardRegion =
          spM.z() > filter.config().centralSeedConfirmationRange.zMaxSeedConf ||
          spM.z() < filter.config().centralSeedConfirmationRange.zMinSeedConf;
      SeedConfirmationRangeConfig seedConfRange =
          isForwardRegion ? filter.config().forwardSeedConfirmationRange
                          : filter.config().centralSeedConfirmationRange;
      // set the minimum number of top SP depending on whether the middle SP is
      // in the central or forward region
      std::size_t nTopSeedConf = spM.r() > seedConfRange.rMaxSeedConf
                                     ? seedConfRange.nTopForLargeR
                                     : seedConfRange.nTopForSmallR;
      // set max bottom radius for seed confirmation
      rMaxSeedConf = seedConfRange.rMaxSeedConf;
      // continue if number of top SPs is smaller than minimum
      if (cache.topDoublets.size() < nTopSeedConf) {
        ACTS_VERBOSE(
            "Number of top SPs is "
            << cache.topDoublets.size()
            << " and is smaller than minimum, moving to next middle candidate");
        continue;
      }
    }

    // create middle-bottom doublets
    cache.bottomDoublets.clear();
    for (std::size_t i = 0; i < bottomSpGroups.size(); ++i) {
      bottomFinder.createSortedDoublets(
          spacePoints, spM, middleSpInfo, bottomSpGroups[i],
          cache.bottomSpOffsets[i], cache.bottomDoublets);
    }

    // no bottom SP found -> try next spM
    if (cache.bottomDoublets.empty()) {
      ACTS_VERBOSE("No compatible Bottoms, moving to next middle candidate");
      continue;
    }

    ACTS_VERBOSE("Candidates: " << cache.bottomDoublets.size()
                                << " bottoms and " << cache.topDoublets.size()
                                << " tops for middle candidate indexed "
                                << spM.index());

    // combine doublets to triplets
    cache.candidatesCollector.clear();
    if (options.useStripMeasurementInfo) {
      createStripTriplets(tripletCuts, rMaxSeedConf, filter, state.filter,
                          cache.filter, spacePoints, spM, cache.bottomDoublets,
                          cache.topDoublets, cache.tripletTopCandidates,
                          cache.candidatesCollector);
    } else {
      createTriplets(cache.tripletCache, tripletCuts, rMaxSeedConf, filter,
                     state.filter, cache.filter, spacePoints, spM,
                     cache.bottomDoublets, cache.topDoublets,
                     cache.tripletTopCandidates, cache.candidatesCollector);
    }

    // retrieve all candidates
    // this collection is already sorted, higher weights first
    const std::size_t numQualitySeeds =
        cache.candidatesCollector.nHighQualityCandidates();
    cache.candidatesCollector.toSortedCandidates(spacePoints,
                                                 cache.sortedCandidates);
    filter.filter1SpFixed(state.filter, spacePoints, cache.sortedCandidates,
                          numQualitySeeds, outputSeeds);
  }
}

void BroadTripletSeedFinder::createTriplets(
    TripletCache& cache, const DerivedTripletCuts& cuts, float rMaxSeedConf,
    const BroadTripletSeedFilter& filter,
    BroadTripletSeedFilter::State& filterState,
    BroadTripletSeedFilter::Cache& filterCache,
    const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
    const DoubletSeedFinder::DoubletsForMiddleSp& bottomDoublets,
    const DoubletSeedFinder::DoubletsForMiddleSp& topDoublets,
    TripletTopCandidates& tripletTopCandidates,
    CandidatesForMiddleSp2& candidatesCollector) {
  const float rM = spM.r();
  const float varianceRM = spM.varianceR();
  const float varianceZM = spM.varianceZ();

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

    auto spB = spacePoints[bottomDoublets.spacePoints[b]];
    const auto& lb = bottomDoublets.linCircles[b];

    float cotThetaB = lb.cotTheta;
    float Vb = lb.V;
    float Ub = lb.U;
    float ErB = lb.Er;
    float iDeltaRB = lb.iDeltaR;

    // 1+(cot^2(theta)) = 1/sin^2(theta)
    float iSinTheta2 = 1 + cotThetaB * cotThetaB;
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
    if (!filter.config().seedConfirmation || spB.r() > rMaxSeedConf) {
      minCompatibleTopSPs = 1;
    }
    if (filter.config().seedConfirmation &&
        candidatesCollector.nHighQualityCandidates() > 0) {
      minCompatibleTopSPs++;
    }

    for (std::size_t indexSortedTop = t0; indexSortedTop < topDoublets.size();
         ++indexSortedTop) {
      const std::size_t t = cache.sortedTops[indexSortedTop];
      auto spT = spacePoints[topDoublets.spacePoints[t]];
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
      float S2 = 1 + A * A;
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
        t0 = indexSortedTop;
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
    if (tripletTopCandidates.size() < minCompatibleTopSPs) {
      continue;
    }

    float zOrigin = spM.z() - rM * lb.cotTheta;
    filter.filter2SpFixed(
        filterState, filterCache, spacePoints, spB.index(), spM.index(),
        tripletTopCandidates.topSpacePoints, tripletTopCandidates.curvatures,
        tripletTopCandidates.impactParameters, zOrigin, candidatesCollector);
  }  // loop on bottoms
}

void BroadTripletSeedFinder::createStripTriplets(
    const DerivedTripletCuts& cuts, float rMaxSeedConf,
    const BroadTripletSeedFilter& filter,
    BroadTripletSeedFilter::State& filterState,
    BroadTripletSeedFilter::Cache& filterCache,
    const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
    const DoubletSeedFinder::DoubletsForMiddleSp& bottomDoublets,
    const DoubletSeedFinder::DoubletsForMiddleSp& topDoublets,
    TripletTopCandidates& tripletTopCandidates,
    CandidatesForMiddleSp2& candidatesCollector) {
  const float rM = spM.r();
  const float cosPhiM = spM.x() / rM;
  const float sinPhiM = spM.y() / rM;
  const float varianceRM = spM.varianceR();
  const float varianceZM = spM.varianceZ();

  // Reserve enough space, in case current capacity is too little
  tripletTopCandidates.resize(topDoublets.size());

  for (std::size_t b = 0; b < bottomDoublets.size(); ++b) {
    auto spB = spacePoints[bottomDoublets.spacePoints[b]];
    const auto& lb = bottomDoublets.linCircles[b];

    float cotThetaB = lb.cotTheta;
    float Vb = lb.V;
    float Ub = lb.U;
    float ErB = lb.Er;
    float iDeltaRB = lb.iDeltaR;

    // 1+(cot^2(theta)) = 1/sin^2(theta)
    float iSinTheta2 = 1 + cotThetaB * cotThetaB;
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
    if (!filter.config().seedConfirmation || spB.r() > rMaxSeedConf) {
      minCompatibleTopSPs = 1;
    }
    if (filter.config().seedConfirmation &&
        candidatesCollector.nHighQualityCandidates() > 0) {
      minCompatibleTopSPs++;
    }

    for (std::size_t t = 0; t < topDoublets.size(); ++t) {
      auto spT = spacePoints[topDoublets.spacePoints[t]];
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
      Eigen::Vector3f positionMiddle = {
          rotationTermsUVtoXY[0] - rotationTermsUVtoXY[1] * A0,
          rotationTermsUVtoXY[0] * A0 + rotationTermsUVtoXY[1],
          zPositionMiddle};

      Eigen::Vector3f rMTransf;
      if (!stripCoordinateCheck(cuts.toleranceParam, spM, positionMiddle,
                                rMTransf)) {
        continue;
      }

      // coordinate transformation and checks for bottom spacepoint
      float B0 = 2 * (Vb - A0 * Ub);
      float Cb = 1 - B0 * lb.y;
      float Sb = A0 + B0 * lb.x;
      Eigen::Vector3f positionBottom = {
          rotationTermsUVtoXY[0] * Cb - rotationTermsUVtoXY[1] * Sb,
          rotationTermsUVtoXY[0] * Sb + rotationTermsUVtoXY[1] * Cb,
          zPositionMiddle};

      Eigen::Vector3f rBTransf;
      if (!stripCoordinateCheck(cuts.toleranceParam, spB, positionBottom,
                                rBTransf)) {
        continue;
      }

      // coordinate transformation and checks for top spacepoint
      float Ct = 1 - B0 * lt.y;
      float St = A0 + B0 * lt.x;
      Eigen::Vector3f positionTop = {
          rotationTermsUVtoXY[0] * Ct - rotationTermsUVtoXY[1] * St,
          rotationTermsUVtoXY[0] * St + rotationTermsUVtoXY[1] * Ct,
          zPositionMiddle};

      Eigen::Vector3f rTTransf;
      if (!stripCoordinateCheck(cuts.toleranceParam, spT, positionTop,
                                rTTransf)) {
        continue;
      }

      // bottom and top coordinates in the spM reference frame
      float xB = rBTransf[0] - rMTransf[0];
      float yB = rBTransf[1] - rMTransf[1];
      float zB = rBTransf[2] - rMTransf[2];
      float xT = rTTransf[0] - rMTransf[0];
      float yT = rTTransf[1] - rMTransf[1];
      float zT = rTTransf[2] - rMTransf[2];

      float iDeltaRB2 = 1 / (xB * xB + yB * yB);
      float iDeltaRT2 = 1 / (xT * xT + yT * yT);

      cotThetaB = -zB * std::sqrt(iDeltaRB2);
      float cotThetaT = zT * std::sqrt(iDeltaRT2);

      // use arithmetic average
      float averageCotTheta = 0.5f * (cotThetaB + cotThetaT);
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
      float S2 = 1 + A * A;
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
    if (tripletTopCandidates.size() < minCompatibleTopSPs) {
      continue;
    }

    float zOrigin = spM.z() - rM * lb.cotTheta;
    filter.filter2SpFixed(
        filterState, filterCache, spacePoints, spB.index(), spM.index(),
        tripletTopCandidates.topSpacePoints, tripletTopCandidates.curvatures,
        tripletTopCandidates.impactParameters, zOrigin, candidatesCollector);
  }  // loop on bottoms
}

}  // namespace Acts::Experimental
