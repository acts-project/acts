// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/DoubletFinder2.hpp"

namespace Acts {

DoubletFinder2::DerivedCuts DoubletFinder2::Cuts::derive(
    float bFieldInZ) const {
  DerivedCuts result;

  static_cast<Cuts&>(result) = *this;

  // bFieldInZ is in (pT/radius) natively, no need for conversion
  const float pTPerHelixRadius = bFieldInZ;
  result.minHelixDiameter2 = std::pow(result.minPt * 2 / pTPerHelixRadius, 2) *
                             result.helixCutTolerance;

  return result;
}

DoubletFinder2::MiddleSpacePointInfo
DoubletFinder2::computeMiddleSpacePointInfo(
    const ConstSpacePointProxy2& spM,
    const SpacePointContainer2::DenseColumn<float>& rColumn) {
  const float rM = spM.extra(rColumn);
  const float uIP = -1. / rM;
  const float cosPhiM = -spM.x() * uIP;
  const float sinPhiM = -spM.y() * uIP;
  const float uIP2 = uIP * uIP;

  return {uIP, uIP2, cosPhiM, sinPhiM};
}

template <DoubletFinder2::SpacePointCandidateType candidate_type>
void DoubletFinder2::createDoublets(
    const DerivedCuts& cuts,
    const SpacePointContainerPointers2& containerPointers,
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
    if (!cuts.interactionPointCut) {
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

// instantiate the template for both candidate types
template void DoubletFinder2::createDoublets<DoubletFinder2::eBottom>(
    const DerivedCuts& cuts,
    const SpacePointContainerPointers2& containerPointers,
    const ConstSpacePointProxy2& middleSp,
    const MiddleSpacePointInfo& middleSpInfo,
    std::span<const SpacePointIndex2> candidateSps,
    Doublets& compatibleDoublets);
template void DoubletFinder2::createDoublets<DoubletFinder2::eTop>(
    const DerivedCuts& cuts,
    const SpacePointContainerPointers2& containerPointers,
    const ConstSpacePointProxy2& middleSp,
    const MiddleSpacePointInfo& middleSpInfo,
    std::span<const SpacePointIndex2> candidateSps,
    Doublets& compatibleDoublets);

}  // namespace Acts
