// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/DoubletSeedFinder.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <stdexcept>

namespace Acts::Experimental {

namespace {

enum class SpacePointCandidateType { Bottom, Top };

/// Iterates over dublets and tests the compatibility by applying a series of
/// cuts that can be tested with only two SPs.
///
/// @tparam candidate_type Type of space point candidate (e.g. Bottom or Top)
/// @tparam interaction_point_cut Whether to apply the interaction point cut
/// @tparam sorted_in_r Whether the space points are sorted in radius
///
/// @param config Doublet cuts that define the compatibility of space points
/// @param spacePoints Space point container to be used
/// @param middleSp Space point candidate to be used as middle SP in a seed
/// @param middleSpInfo Information about the middle space point
/// @param candidateSps Group of space points to be used as candidates for
///                     middle SP in a seed
/// @param candidateOffset Offset in the candidateSps to start from
/// @param compatibleDoublets Output container for compatible doublets
template <SpacePointCandidateType candidateType, bool interactionPointCut,
          bool sortedInR>
void createDoubletsImpl(
    const DoubletSeedFinder::DerivedConfig& config,
    const SpacePointContainer2& spacePoints,
    const ConstSpacePointProxy2& middleSp,
    const DoubletSeedFinder::MiddleSpInfo& middleSpInfo,
    std::span<const SpacePointIndex2> candidateSps,
    std::size_t& candidateOffset,
    DoubletSeedFinder::DoubletsForMiddleSp& compatibleDoublets) {
  constexpr bool isBottomCandidate =
      candidateType == SpacePointCandidateType::Bottom;

  const float impactMax =
      isBottomCandidate ? -config.impactMax : config.impactMax;

  const float xM = middleSp.x();
  const float yM = middleSp.y();
  const float zM = middleSp.z();
  const float rM = middleSp.r();
  const float varianceZM = middleSp.varianceZ();
  const float varianceRM = middleSp.varianceR();

  // equivalent to impactMax / (rM * rM);
  const float vIPAbs = impactMax * middleSpInfo.uIP2;

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
    const float varianceZO = otherSp.varianceZ();
    const float varianceRO = otherSp.varianceR();

    return iDeltaR2 * ((varianceZM + varianceZO) +
                       (cotTheta * cotTheta) * (varianceRM + varianceRO));
  };

  if constexpr (sortedInR) {
    // find the first SP inside the radius region of interest and update
    // the iterator so we don't need to look at the other SPs again
    for (; candidateOffset < candidateSps.size(); ++candidateOffset) {
      ConstSpacePointProxy2 otherSp =
          spacePoints[candidateSps[candidateOffset]];

      if constexpr (isBottomCandidate) {
        // if r-distance is too big, try next SP in bin
        if (rM - otherSp.r() <= config.deltaRMax) {
          break;
        }
      } else {
        // if r-distance is too small, try next SP in bin
        if (otherSp.r() - rM >= config.deltaRMin) {
          break;
        }
      }
    }
  }

  for (SpacePointIndex2 otherSpIndex : candidateSps.subspan(candidateOffset)) {
    ConstSpacePointProxy2 otherSp = spacePoints[otherSpIndex];

    const float xO = otherSp.x();
    const float yO = otherSp.y();
    const float zO = otherSp.z();
    const float rO = otherSp.r();

    if constexpr (isBottomCandidate) {
      deltaR = rM - rO;

      if constexpr (sortedInR) {
        // if r-distance is too small we are done
        if (deltaR < config.deltaRMin) {
          break;
        }
      }
    } else {
      deltaR = rO - rM;

      if constexpr (sortedInR) {
        // if r-distance is too big we are done
        if (deltaR > config.deltaRMax) {
          break;
        }
      }
    }

    if constexpr (!sortedInR) {
      if (outsideRangeCheck(deltaR, config.deltaRMin, config.deltaRMax)) {
        continue;
      }
    }

    if constexpr (isBottomCandidate) {
      deltaZ = zM - zO;
    } else {
      deltaZ = zO - zM;
    }

    if (outsideRangeCheck(deltaZ, config.deltaZMin, config.deltaZMax)) {
      continue;
    }

    // the longitudinal impact parameter zOrigin is defined as (zM - rM *
    // cotTheta) where cotTheta is the ratio Z/R (forward angle) of space
    // point duplet but instead we calculate (zOrigin * deltaR) and multiply
    // collisionRegion by deltaR to avoid divisions
    const float zOriginTimesDeltaR = zM * deltaR - rM * deltaZ;
    // check if duplet origin on z axis within collision region
    if (outsideRangeCheck(zOriginTimesDeltaR,
                          config.collisionRegionMin * deltaR,
                          config.collisionRegionMax * deltaR)) {
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
      if (outsideRangeCheck(deltaZ, -config.cotThetaMax * deltaR,
                            config.cotThetaMax * deltaR)) {
        continue;
      }

      // transform SP coordinates to the u-v reference frame
      const float deltaX = xO - xM;
      const float deltaY = yO - yM;

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
    const float deltaX = xO - xM;
    const float deltaY = yO - yM;

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
      if (outsideRangeCheck(deltaZ, -config.cotThetaMax * deltaR,
                            config.cotThetaMax * deltaR)) {
        continue;
      }

      const float iDeltaR = std::sqrt(iDeltaR2);
      const float cotTheta = deltaZ * iDeltaR;

      // discard bottom-middle doublets in a certain (r, eta) region according
      // to detector specific cuts
      if constexpr (isBottomCandidate) {
        if (config.experimentCuts.connected() &&
            !config.experimentCuts(rO, cotTheta)) {
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
    if ((bCoef * bCoef) * config.minHelixDiameter2 > 1 + aCoef * aCoef) {
      continue;
    }

    // check if duplet cotTheta is within the region of interest
    // cotTheta is defined as (deltaZ / deltaR) but instead we multiply
    // cotThetaMax by deltaR to avoid division
    if (outsideRangeCheck(deltaZ, -config.cotThetaMax * deltaR,
                          config.cotThetaMax * deltaR)) {
      continue;
    }

    const float iDeltaR = std::sqrt(iDeltaR2);
    const float cotTheta = deltaZ * iDeltaR;

    // discard bottom-middle doublets in a certain (r, eta) region according
    // to detector specific cuts
    if constexpr (isBottomCandidate) {
      if (config.experimentCuts.connected() &&
          !config.experimentCuts(rO, cotTheta)) {
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

DoubletSeedFinder::DerivedConfig::DerivedConfig(const Config& cfg,
                                                float bFieldInZ_)
    : Config(cfg), bFieldInZ(bFieldInZ_) {
  // bFieldInZ is in (pT/radius) natively, no need for conversion
  const float pTPerHelixRadius = bFieldInZ;
  minHelixDiameter2 = square(minPt * 2 / pTPerHelixRadius) * helixCutTolerance;
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

DoubletSeedFinder::DoubletSeedFinder(const DerivedConfig& cfg) : m_cfg(cfg) {}

void DoubletSeedFinder::createDoublets(
    const SpacePointContainer2& spacePoints,
    const ConstSpacePointProxy2& middleSp, const MiddleSpInfo& middleSpInfo,
    std::span<const SpacePointIndex2> candidateSps,
    DoubletsForMiddleSp& compatibleDoublets) const {
  std::size_t candidateOffset = 0;

  if (m_cfg.candidateDirection == Direction::Backward() &&
      m_cfg.interactionPointCut) {
    return createDoubletsImpl<SpacePointCandidateType::Bottom, true, false>(
        config(), spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  }

  if (m_cfg.candidateDirection == Direction::Backward() &&
      !m_cfg.interactionPointCut) {
    return createDoubletsImpl<SpacePointCandidateType::Bottom, false, false>(
        config(), spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  }

  if (m_cfg.candidateDirection == Direction::Forward() &&
      m_cfg.interactionPointCut) {
    return createDoubletsImpl<SpacePointCandidateType::Top, true, false>(
        config(), spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  }

  if (m_cfg.candidateDirection == Direction::Forward() &&
      !m_cfg.interactionPointCut) {
    return createDoubletsImpl<SpacePointCandidateType::Top, false, false>(
        config(), spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  }

  throw std::logic_error("DoubletSeedFinder: unhandled configuration");
}

void DoubletSeedFinder::createSortedDoublets(
    const SpacePointContainer2& spacePoints,
    const ConstSpacePointProxy2& middleSp, const MiddleSpInfo& middleSpInfo,
    std::span<const SpacePointIndex2> candidateSps,
    std::size_t& candidateOffset,
    DoubletsForMiddleSp& compatibleDoublets) const {
  if (m_cfg.candidateDirection == Direction::Backward() &&
      m_cfg.interactionPointCut) {
    return createDoubletsImpl<SpacePointCandidateType::Bottom, true, true>(
        config(), spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  }

  if (m_cfg.candidateDirection == Direction::Backward() &&
      !m_cfg.interactionPointCut) {
    return createDoubletsImpl<SpacePointCandidateType::Bottom, false, true>(
        config(), spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  }

  if (m_cfg.candidateDirection == Direction::Forward() &&
      m_cfg.interactionPointCut) {
    return createDoubletsImpl<SpacePointCandidateType::Top, true, true>(
        config(), spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  }

  if (m_cfg.candidateDirection == Direction::Forward() &&
      !m_cfg.interactionPointCut) {
    return createDoubletsImpl<SpacePointCandidateType::Top, false, true>(
        config(), spacePoints, middleSp, middleSpInfo, candidateSps,
        candidateOffset, compatibleDoublets);
  }

  throw std::logic_error("DoubletSeedFinder: unhandled configuration");
}

}  // namespace Acts::Experimental
