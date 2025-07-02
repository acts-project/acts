// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/DoubletSeedFinder.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"

namespace Acts::Experimental {

namespace {

enum class SpacePointCandidateType { eBottom, eTop };

/// Iterates over dublets and tests the compatibility by applying a series of
/// cuts that can be tested with only two SPs.
///
/// @tparam candidate_type Type of space point candidate (e.g. Bottom or Top)
/// @tparam interaction_point_cut Whether to apply the interaction point cut
/// @tparam sorted_in_r Whether the space points are sorted in radius
///
/// @param cuts Doublet cuts that define the compatibility of space points
/// @param spacePoints Space point container to be used
/// @param middleSp Space point candidate to be used as middle SP in a seed
/// @param middleSpInfo Information about the middle space point
/// @param candidateSps Group of space points to be used as candidates for
///                     middle SP in a seed
/// @param candidateOffset Offset in the candidateSps to start from
/// @param compatibleDoublets Output container for compatible doublets
template <SpacePointCandidateType candidateType, bool interactionPointCut,
          bool sortedInR>
void createDoublets(
    const DoubletSeedFinder::DerivedCuts& cuts,
    const SpacePointContainer2& spacePoints,
    const ConstSpacePointProxy2& middleSp,
    const DoubletSeedFinder::MiddleSpInfo& middleSpInfo,
    std::span<const SpacePointIndex2> candidateSps,
    std::size_t& candidateOffset,
    DoubletSeedFinder::DoubletsForMiddleSp& compatibleDoublets) {
  constexpr bool isBottomCandidate =
      candidateType == SpacePointCandidateType::eBottom;

  const float impactMax = isBottomCandidate ? -cuts.impactMax : cuts.impactMax;

  const float xM = middleSp.x();
  const float yM = middleSp.y();
  const float zM = middleSp.z();
  const float rM = middleSp.r();

  // equivalent to impactMax / (rM * rM);
  float vIPAbs = impactMax * middleSpInfo.uIP2;

  float deltaR = 0.;
  float deltaZ = 0.;

  const auto outsideRangeCheck = [](const float value, const float min,
                                    const float max) {
    // intentionally using `|` after profiling. faster due to better branch
    // prediction
    return static_cast<bool>(static_cast<int>(value < min) |
                             static_cast<int>(value > max));
  };

  const auto calculateError = [&](const ConstSpacePointProxy2& otherSp,
                                  float iDeltaR2, float cotTheta) {
    using enum SpacePointKnownExtraColumn;

    // TOD use some reasonable defaults
    float varianceZM =
        spacePoints.hasExtraColumns(VarianceZ) ? middleSp.varianceZ() : 0;
    float varianceZO =
        spacePoints.hasExtraColumns(VarianceZ) ? otherSp.varianceZ() : 0;
    float varianceRM =
        spacePoints.hasExtraColumns(VarianceR) ? middleSp.varianceR() : 0;
    float varianceRO =
        spacePoints.hasExtraColumns(VarianceR) ? otherSp.varianceR() : 0;

    return iDeltaR2 * ((varianceZM + varianceZO) +
                       (cotTheta * cotTheta) * (varianceRM + varianceRO));
  };

  if constexpr (sortedInR) {
    // find the first SP inside the radius region of interest and update
    // the iterator so we don't need to look at the other SPs again
    for (; candidateOffset < candidateSps.size(); ++candidateOffset) {
      ConstSpacePointProxy2 otherSp =
          spacePoints.at(candidateSps[candidateOffset]);

      if constexpr (isBottomCandidate) {
        // if r-distance is too big, try next SP in bin
        if (rM - otherSp.r() <= cuts.deltaRMax) {
          break;
        }
      } else {
        // if r-distance is too small, try next SP in bin
        if (otherSp.r() - rM >= cuts.deltaRMin) {
          break;
        }
      }
    }
  }

  for (SpacePointIndex2 otherSpIndex : candidateSps.subspan(candidateOffset)) {
    ConstSpacePointProxy2 otherSp = spacePoints.at(otherSpIndex);

    if constexpr (isBottomCandidate) {
      deltaR = rM - otherSp.r();

      if constexpr (sortedInR) {
        // if r-distance is too small we are done
        if (deltaR < cuts.deltaRMin) {
          break;
        }
      }
    } else {
      deltaR = otherSp.r() - rM;

      if constexpr (sortedInR) {
        // if r-distance is too big we are done
        if (deltaR > cuts.deltaRMax) {
          break;
        }
      }
    }

    if constexpr (!sortedInR) {
      if (outsideRangeCheck(deltaR, cuts.deltaRMin, cuts.deltaRMax)) {
        continue;
      }
    }

    if constexpr (isBottomCandidate) {
      deltaZ = zM - otherSp.z();
    } else {
      deltaZ = otherSp.z() - zM;
    }

    if (outsideRangeCheck(deltaZ, cuts.deltaZMin, cuts.deltaZMax)) {
      continue;
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
    if constexpr (!interactionPointCut) {
      // check if duplet cotTheta is within the region of interest
      // cotTheta is defined as (deltaZ / deltaR) but instead we multiply
      // cotThetaMax by deltaR to avoid division
      if (outsideRangeCheck(deltaZ, -cuts.cotThetaMax * deltaR,
                            cuts.cotThetaMax * deltaR)) {
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
      const float iDeltaR2 = 1 / deltaR2;

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
    const float iDeltaR2 = 1 / deltaR2;

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
            !cuts.experimentCuts(otherSp.r(), cotTheta)) {
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
          !cuts.experimentCuts(otherSp.r(), cotTheta)) {
        continue;
      }
    }

    const float er = calculateError(otherSp, iDeltaR2, cotTheta);

    // fill output vectors
    compatibleDoublets.emplace_back(
        otherSp.index(), {cotTheta, iDeltaR, er, uT, vT, xNewFrame, yNewFrame});
  }
}

}  // namespace

DoubletSeedFinder::DerivedCuts DoubletSeedFinder::Cuts::derive(
    float bFieldInZ) const {
  DerivedCuts result;

  static_cast<Cuts&>(result) = *this;

  // bFieldInZ is in (pT/radius) natively, no need for conversion
  const float pTPerHelixRadius = bFieldInZ;
  result.minHelixDiameter2 = std::powf(result.minPt * 2 / pTPerHelixRadius, 2) *
                             result.helixCutTolerance;

  return result;
}

DoubletSeedFinder::MiddleSpInfo DoubletSeedFinder::computeMiddleSpInfo(
    const ConstSpacePointProxy2& spM) {
  const float rM = spM.r();
  const float uIP = -1 / rM;
  const float cosPhiM = -spM.x() * uIP;
  const float sinPhiM = -spM.y() * uIP;
  const float uIP2 = uIP * uIP;

  return {uIP, uIP2, cosPhiM, sinPhiM};
}

void DoubletSeedFinder::createBottomDoublets(
    const DerivedCuts& cuts, const SpacePointContainer2& spacePoints,
    const ConstSpacePointProxy2& middleSp, const MiddleSpInfo& middleSpInfo,
    std::span<const SpacePointIndex2> candidateSps,
    DoubletsForMiddleSp& compatibleDoublets) {
  std::size_t candidateOffset = 0;
  if (cuts.interactionPointCut) {
    return createDoublets<SpacePointCandidateType::eBottom, true, false>(
        cuts, spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  } else {
    return createDoublets<SpacePointCandidateType::eBottom, false, false>(
        cuts, spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  }
}

void DoubletSeedFinder::createTopDoublets(
    const DerivedCuts& cuts, const SpacePointContainer2& spacePoints,
    const ConstSpacePointProxy2& middleSp, const MiddleSpInfo& middleSpInfo,
    std::span<const SpacePointIndex2> candidateSps,
    DoubletsForMiddleSp& compatibleDoublets) {
  std::size_t candidateOffset = 0;
  if (cuts.interactionPointCut) {
    return createDoublets<SpacePointCandidateType::eTop, true, false>(
        cuts, spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  } else {
    return createDoublets<SpacePointCandidateType::eTop, false, false>(
        cuts, spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  }
}

void DoubletSeedFinder::createSortedBottomDoublets(
    const DerivedCuts& cuts, const SpacePointContainer2& spacePoints,
    const ConstSpacePointProxy2& middleSp, const MiddleSpInfo& middleSpInfo,
    std::span<const SpacePointIndex2> candidateSps,
    std::size_t& candidateOffset, DoubletsForMiddleSp& compatibleDoublets) {
  if (cuts.interactionPointCut) {
    return createDoublets<SpacePointCandidateType::eBottom, true, true>(
        cuts, spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  } else {
    return createDoublets<SpacePointCandidateType::eBottom, false, true>(
        cuts, spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  }
}

void DoubletSeedFinder::createSortedTopDoublets(
    const DerivedCuts& cuts, const SpacePointContainer2& spacePoints,
    const ConstSpacePointProxy2& middleSp, const MiddleSpInfo& middleSpInfo,
    std::span<const SpacePointIndex2> candidateSps,
    std::size_t& candidateOffset, DoubletsForMiddleSp& compatibleDoublets) {
  if (cuts.interactionPointCut) {
    return createDoublets<SpacePointCandidateType::eTop, true, true>(
        cuts, spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  } else {
    return createDoublets<SpacePointCandidateType::eTop, false, true>(
        cuts, spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  }
}

}  // namespace Acts::Experimental
