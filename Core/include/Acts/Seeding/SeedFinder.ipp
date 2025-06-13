// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/SeedFinder.hpp"

#include <algorithm>
#include <cmath>

namespace Acts {

template <typename external_spacepoint_t, typename grid_t, typename platform_t>
SeedFinder<external_spacepoint_t, grid_t, platform_t>::SeedFinder(
    const SeedFinderConfig<external_spacepoint_t>& config,
    std::unique_ptr<const Logger> logger)
    : m_config(config), m_logger(std::move(logger)) {
  if (std::isnan(config.deltaRMaxTopSP)) {
    throw std::runtime_error("Value of deltaRMaxTopSP was not initialised");
  }
  if (std::isnan(config.deltaRMinTopSP)) {
    throw std::runtime_error("Value of deltaRMinTopSP was not initialised");
  }
  if (std::isnan(config.deltaRMaxBottomSP)) {
    throw std::runtime_error("Value of deltaRMaxBottomSP was not initialised");
  }
  if (std::isnan(config.deltaRMinBottomSP)) {
    throw std::runtime_error("Value of deltaRMinBottomSP was not initialised");
  }
}

template <typename external_spacepoint_t, typename grid_t, typename platform_t>
template <typename container_t, GridBinCollection sp_range_t>
  requires CollectionStoresSeedsTo<container_t, external_spacepoint_t, 3ul>
void SeedFinder<external_spacepoint_t, grid_t, platform_t>::createSeedsForGroup(
    const SeedFinderOptions& options, SeedingState& state, const grid_t& grid,
    container_t& outputCollection, const sp_range_t& bottomSPsIdx,
    const std::size_t middleSPsIdx, const sp_range_t& topSPsIdx,
    const Range1D<float>& rMiddleSPRange) const {
  // This is used for seed filtering later
  const std::size_t max_num_seeds_per_spm =
      m_config.seedFilter->getSeedFilterConfig().maxSeedsPerSpMConf;
  const std::size_t max_num_quality_seeds_per_spm =
      m_config.seedFilter->getSeedFilterConfig().maxQualitySeedsPerSpMConf;

  state.candidatesCollector.setMaxElements(max_num_seeds_per_spm,
                                           max_num_quality_seeds_per_spm);

  // If there are no bottom or top bins, just return and waste no time
  if (bottomSPsIdx.size() == 0 || topSPsIdx.size() == 0) {
    return;
  }

  // Get the middle space point candidates
  const std::vector<const external_spacepoint_t*>& middleSPs =
      grid.at(middleSPsIdx);
  // Return if somehow there are no middle sp candidates
  if (middleSPs.empty()) {
    return;
  }

  // neighbours
  // clear previous results
  state.bottomNeighbours.clear();
  state.topNeighbours.clear();

  // Fill
  // bottoms
  for (const std::size_t idx : bottomSPsIdx) {
    // Only add an entry if the bin has entries
    if (grid.at(idx).size() == 0) {
      continue;
    }
    state.bottomNeighbours.emplace_back(
        grid, idx, middleSPs.front()->radius() - m_config.deltaRMaxBottomSP);
  }
  // if no bottom candidates, then no need to proceed
  if (state.bottomNeighbours.size() == 0) {
    return;
  }

  // tops
  for (const std::size_t idx : topSPsIdx) {
    // Only add an entry if the bin has entries
    if (grid.at(idx).size() == 0) {
      continue;
    }
    state.topNeighbours.emplace_back(
        grid, idx, middleSPs.front()->radius() + m_config.deltaRMinTopSP);
  }
  // if no top candidates, then no need to proceed
  if (state.topNeighbours.size() == 0) {
    return;
  }

  // we compute this here since all middle space point candidates belong to the
  // same z-bin
  auto [minRadiusRangeForMiddle, maxRadiusRangeForMiddle] =
      retrieveRadiusRangeForMiddle(*middleSPs.front(), rMiddleSPRange);
  ACTS_VERBOSE("Current global bin: " << middleSPsIdx << ", z value of "
                                      << middleSPs.front()->z());
  ACTS_VERBOSE("Validity range (radius) for the middle space point is ["
               << minRadiusRangeForMiddle << ", " << maxRadiusRangeForMiddle
               << "]");

  for (const external_spacepoint_t* spM : middleSPs) {
    const float rM = spM->radius();

    // check if spM is outside our radial region of interest
    if (rM < minRadiusRangeForMiddle) {
      continue;
    }
    if (rM > maxRadiusRangeForMiddle) {
      // break because SPs are sorted in r
      break;
    }

    const float zM = spM->z();
    const float uIP = -1. / rM;
    const float cosPhiM = -spM->x() * uIP;
    const float sinPhiM = -spM->y() * uIP;
    const float uIP2 = uIP * uIP;

    // Iterate over middle-top dublets
    getCompatibleDoublets<SpacePointCandidateType::eTop>(
        options, grid, state.spacePointMutableData, state.topNeighbours, *spM,
        state.linCircleTop, state.compatTopSP, m_config.deltaRMinTopSP,
        m_config.deltaRMaxTopSP, uIP, uIP2, cosPhiM, sinPhiM);

    // no top SP found -> try next spM
    if (state.compatTopSP.empty()) {
      ACTS_VERBOSE("No compatible Tops, moving to next middle candidate");
      continue;
    }

    // apply cut on the number of top SP if seedConfirmation is true
    SeedFilterState seedFilterState;
    if (m_config.seedConfirmation) {
      // check if middle SP is in the central or forward region
      SeedConfirmationRangeConfig seedConfRange =
          (zM > m_config.centralSeedConfirmationRange.zMaxSeedConf ||
           zM < m_config.centralSeedConfirmationRange.zMinSeedConf)
              ? m_config.forwardSeedConfirmationRange
              : m_config.centralSeedConfirmationRange;
      // set the minimum number of top SP depending on whether the middle SP is
      // in the central or forward region
      seedFilterState.nTopSeedConf = rM > seedConfRange.rMaxSeedConf
                                         ? seedConfRange.nTopForLargeR
                                         : seedConfRange.nTopForSmallR;
      // set max bottom radius for seed confirmation
      seedFilterState.rMaxSeedConf = seedConfRange.rMaxSeedConf;
      // continue if number of top SPs is smaller than minimum
      if (state.compatTopSP.size() < seedFilterState.nTopSeedConf) {
        ACTS_VERBOSE(
            "Number of top SPs is "
            << state.compatTopSP.size()
            << " and is smaller than minimum, moving to next middle candidate");
        continue;
      }
    }

    // Iterate over middle-bottom dublets
    getCompatibleDoublets<SpacePointCandidateType::eBottom>(
        options, grid, state.spacePointMutableData, state.bottomNeighbours,
        *spM, state.linCircleBottom, state.compatBottomSP,
        m_config.deltaRMinBottomSP, m_config.deltaRMaxBottomSP, uIP, uIP2,
        cosPhiM, sinPhiM);

    // no bottom SP found -> try next spM
    if (state.compatBottomSP.empty()) {
      ACTS_VERBOSE("No compatible Bottoms, moving to next middle candidate");
      continue;
    }

    ACTS_VERBOSE("Candidates: " << state.compatBottomSP.size()
                                << " bottoms and " << state.compatTopSP.size()
                                << " tops for middle candidate indexed "
                                << spM->index());
    // filter candidates
    if (m_config.useDetailedDoubleMeasurementInfo) {
      filterCandidates<DetectorMeasurementInfo::eDetailed>(
          *spM, options, seedFilterState, state);
    } else {
      filterCandidates<DetectorMeasurementInfo::eDefault>(
          *spM, options, seedFilterState, state);
    }

    m_config.seedFilter->filterSeeds_1SpFixed(state.spacePointMutableData,
                                              state.candidatesCollector,
                                              outputCollection);

  }  // loop on mediums
}

template <typename external_spacepoint_t, typename grid_t, typename platform_t>
template <SpacePointCandidateType candidateType, typename out_range_t>
inline void
SeedFinder<external_spacepoint_t, grid_t, platform_t>::getCompatibleDoublets(
    const SeedFinderOptions& options, const grid_t& grid,
    SpacePointMutableData& mutableData,
    boost::container::small_vector<
        Neighbour<grid_t>, detail::ipow(3, grid_t::DIM)>& otherSPsNeighbours,
    const external_spacepoint_t& mediumSP, std::vector<LinCircle>& linCircleVec,
    out_range_t& outVec, const float deltaRMinSP, const float deltaRMaxSP,
    const float uIP, const float uIP2, const float cosPhiM,
    const float sinPhiM) const {
  float impactMax = m_config.impactMax;

  constexpr bool isBottomCandidate =
      candidateType == SpacePointCandidateType::eBottom;

  if constexpr (isBottomCandidate) {
    impactMax = -impactMax;
  }

  outVec.clear();
  linCircleVec.clear();

  // get number of neighbour SPs
  std::size_t nsp = 0;
  for (const auto& otherSPCol : otherSPsNeighbours) {
    nsp += grid.at(otherSPCol.index).size();
  }

  linCircleVec.reserve(nsp);
  outVec.reserve(nsp);

  const float rM = mediumSP.radius();
  const float xM = mediumSP.x();
  const float yM = mediumSP.y();
  const float zM = mediumSP.z();
  const float varianceRM = mediumSP.varianceR();
  const float varianceZM = mediumSP.varianceZ();

  float vIPAbs = 0;
  if (m_config.interactionPointCut) {
    // equivalent to m_config.impactMax / (rM * rM);
    vIPAbs = impactMax * uIP2;
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

  for (auto& otherSPCol : otherSPsNeighbours) {
    const std::vector<const external_spacepoint_t*>& otherSPs =
        grid.at(otherSPCol.index);
    if (otherSPs.empty()) {
      continue;
    }

    // we make a copy of the iterator here since we need it to remain
    // the same in the Neighbour object
    auto min_itr = otherSPCol.itr;

    // find the first SP inside the radius region of interest and update
    // the iterator so we don't need to look at the other SPs again
    for (; min_itr != otherSPs.end(); ++min_itr) {
      const external_spacepoint_t* otherSP = *min_itr;
      if constexpr (candidateType == SpacePointCandidateType::eBottom) {
        // if r-distance is too big, try next SP in bin
        if ((rM - otherSP->radius()) <= deltaRMaxSP) {
          break;
        }
      } else {
        // if r-distance is too small, try next SP in bin
        if ((otherSP->radius() - rM) >= deltaRMinSP) {
          break;
        }
      }
    }
    // We update the iterator in the Neighbour object
    // that mean that we have changed the middle space point
    // and the lower bound has moved accordingly
    otherSPCol.itr = min_itr;

    for (; min_itr != otherSPs.end(); ++min_itr) {
      const external_spacepoint_t* otherSP = *min_itr;

      if constexpr (isBottomCandidate) {
        deltaR = (rM - otherSP->radius());

        // if r-distance is too small, try next SP in bin
        if (deltaR < deltaRMinSP) {
          break;
        }
      } else {
        deltaR = (otherSP->radius() - rM);

        // if r-distance is too big, try next SP in bin
        if (deltaR > deltaRMaxSP) {
          break;
        }
      }

      if constexpr (isBottomCandidate) {
        deltaZ = (zM - otherSP->z());
      } else {
        deltaZ = (otherSP->z() - zM);
      }

      // the longitudinal impact parameter zOrigin is defined as (zM - rM *
      // cotTheta) where cotTheta is the ratio Z/R (forward angle) of space
      // point duplet but instead we calculate (zOrigin * deltaR) and multiply
      // collisionRegion by deltaR to avoid divisions
      const float zOriginTimesDeltaR = (zM * deltaR - rM * deltaZ);
      // check if duplet origin on z axis within collision region
      if (outsideRangeCheck(zOriginTimesDeltaR,
                            m_config.collisionRegionMin * deltaR,
                            m_config.collisionRegionMax * deltaR)) {
        continue;
      }

      // if interactionPointCut is false we apply z cuts before coordinate
      // transformation to avoid unnecessary calculations. If
      // interactionPointCut is true we apply the curvature cut first because it
      // is more frequent but requires the coordinate transformation
      if (!m_config.interactionPointCut) {
        // check if duplet cotTheta is within the region of interest
        // cotTheta is defined as (deltaZ / deltaR) but instead we multiply
        // cotThetaMax by deltaR to avoid division
        if (outsideRangeCheck(deltaZ, -m_config.cotThetaMax * deltaR,
                              m_config.cotThetaMax * deltaR)) {
          continue;
        }
        // if z-distance between SPs is within max and min values
        if (outsideRangeCheck(deltaZ, -m_config.deltaZMax,
                              m_config.deltaZMax)) {
          continue;
        }

        // transform SP coordinates to the u-v reference frame
        const float deltaX = otherSP->x() - xM;
        const float deltaY = otherSP->y() - yM;

        const float xNewFrame = deltaX * cosPhiM + deltaY * sinPhiM;
        const float yNewFrame = deltaY * cosPhiM - deltaX * sinPhiM;

        const float deltaR2 = (deltaX * deltaX + deltaY * deltaY);
        const float iDeltaR2 = 1. / deltaR2;

        const float uT = xNewFrame * iDeltaR2;
        const float vT = yNewFrame * iDeltaR2;

        const float iDeltaR = std::sqrt(iDeltaR2);
        const float cotTheta = deltaZ * iDeltaR;

        const float Er =
            ((varianceZM + otherSP->varianceZ()) +
             (cotTheta * cotTheta) * (varianceRM + otherSP->varianceR())) *
            iDeltaR2;

        // fill output vectors
        linCircleVec.emplace_back(cotTheta, iDeltaR, Er, uT, vT, xNewFrame,
                                  yNewFrame);

        mutableData.setDeltaR(otherSP->index(),
                              std::sqrt(deltaR2 + (deltaZ * deltaZ)));
        outVec.push_back(otherSP);
        continue;
      }

      // transform SP coordinates to the u-v reference frame
      const float deltaX = otherSP->x() - xM;
      const float deltaY = otherSP->y() - yM;

      const float xNewFrame = deltaX * cosPhiM + deltaY * sinPhiM;
      const float yNewFrame = deltaY * cosPhiM - deltaX * sinPhiM;

      const float deltaR2 = (deltaX * deltaX + deltaY * deltaY);
      const float iDeltaR2 = 1. / deltaR2;

      const float uT = xNewFrame * iDeltaR2;
      const float vT = yNewFrame * iDeltaR2;

      // We check the interaction point by evaluating the minimal distance
      // between the origin and the straight line connecting the two points in
      // the doublets. Using a geometric similarity, the Im is given by
      // yNewFrame * rM / deltaR <= m_config.impactMax
      // However, we make here an approximation of the impact parameter
      // which is valid under the assumption yNewFrame / xNewFrame is small
      // The correct computation would be:
      // yNewFrame * yNewFrame * rM * rM <= m_config.impactMax *
      // m_config.impactMax * deltaR2
      if (std::abs(rM * yNewFrame) <= impactMax * xNewFrame) {
        // check if duplet cotTheta is within the region of interest
        // cotTheta is defined as (deltaZ / deltaR) but instead we multiply
        // cotThetaMax by deltaR to avoid division
        if (outsideRangeCheck(deltaZ, -m_config.cotThetaMax * deltaR,
                              m_config.cotThetaMax * deltaR)) {
          continue;
        }

        const float iDeltaR = std::sqrt(iDeltaR2);
        const float cotTheta = deltaZ * iDeltaR;

        // discard bottom-middle dublets in a certain (r, eta) region according
        // to detector specific cuts
        if constexpr (isBottomCandidate) {
          if (!m_config.experimentCuts(otherSP->radius(), cotTheta)) {
            continue;
          }
        }

        const float Er =
            ((varianceZM + otherSP->varianceZ()) +
             (cotTheta * cotTheta) * (varianceRM + otherSP->varianceR())) *
            iDeltaR2;

        // fill output vectors
        linCircleVec.emplace_back(cotTheta, iDeltaR, Er, uT, vT, xNewFrame,
                                  yNewFrame);
        mutableData.setDeltaR(otherSP->index(),
                              std::sqrt(deltaR2 + (deltaZ * deltaZ)));
        outVec.emplace_back(otherSP);
        continue;
      }

      // in the rotated frame the interaction point is positioned at x = -rM
      // and y ~= impactParam
      const float vIP = (yNewFrame > 0.) ? -vIPAbs : vIPAbs;

      // we can obtain aCoef as the slope dv/du of the linear function,
      // estimated using du and dv between the two SP bCoef is obtained by
      // inserting aCoef into the linear equation
      const float aCoef = (vT - vIP) / (uT - uIP);
      const float bCoef = vIP - aCoef * uIP;
      // the distance of the straight line from the origin (radius of the
      // circle) is related to aCoef and bCoef by d^2 = bCoef^2 / (1 +
      // aCoef^2) = 1 / (radius^2) and we can apply the cut on the curvature
      if ((bCoef * bCoef) * options.minHelixDiameter2 > (1 + aCoef * aCoef)) {
        continue;
      }

      // check if duplet cotTheta is within the region of interest
      // cotTheta is defined as (deltaZ / deltaR) but instead we multiply
      // cotThetaMax by deltaR to avoid division
      if (outsideRangeCheck(deltaZ, -m_config.cotThetaMax * deltaR,
                            m_config.cotThetaMax * deltaR)) {
        continue;
      }

      const float iDeltaR = std::sqrt(iDeltaR2);
      const float cotTheta = deltaZ * iDeltaR;

      // discard bottom-middle dublets in a certain (r, eta) region according
      // to detector specific cuts
      if constexpr (isBottomCandidate) {
        if (!m_config.experimentCuts(otherSP->radius(), cotTheta)) {
          continue;
        }
      }

      const float Er =
          ((varianceZM + otherSP->varianceZ()) +
           (cotTheta * cotTheta) * (varianceRM + otherSP->varianceR())) *
          iDeltaR2;

      // fill output vectors
      linCircleVec.emplace_back(cotTheta, iDeltaR, Er, uT, vT, xNewFrame,
                                yNewFrame);

      mutableData.setDeltaR(otherSP->index(),
                            std::sqrt(deltaR2 + (deltaZ * deltaZ)));
      outVec.emplace_back(otherSP);
    }
  }
}

template <typename external_spacepoint_t, typename grid_t, typename platform_t>
template <DetectorMeasurementInfo detailedMeasurement>
inline void
SeedFinder<external_spacepoint_t, grid_t, platform_t>::filterCandidates(
    const external_spacepoint_t& spM, const SeedFinderOptions& options,
    SeedFilterState& seedFilterState, SeedingState& state) const {
  const float rM = spM.radius();
  const float cosPhiM = spM.x() / rM;
  const float sinPhiM = spM.y() / rM;
  const float varianceRM = spM.varianceR();
  const float varianceZM = spM.varianceZ();

  std::size_t numTopSp = state.compatTopSP.size();

  // sort: make index vector
  std::vector<std::uint32_t> sortedBottoms(state.compatBottomSP.size());
  for (std::uint32_t i = 0; i < sortedBottoms.size(); ++i) {
    sortedBottoms[i] = i;
  }
  std::vector<std::uint32_t> sortedTops(state.linCircleTop.size());
  for (std::uint32_t i = 0; i < sortedTops.size(); ++i) {
    sortedTops[i] = i;
  }

  if constexpr (detailedMeasurement == DetectorMeasurementInfo::eDefault) {
    std::vector<float> cotThetaBottoms(state.compatBottomSP.size());
    for (std::uint32_t i = 0; i < sortedBottoms.size(); ++i) {
      cotThetaBottoms[i] = state.linCircleBottom[i].cotTheta;
    }
    std::ranges::sort(sortedBottoms, {}, [&](const std::uint32_t s) {
      return cotThetaBottoms[s];
    });

    std::vector<float> cotThetaTops(state.linCircleTop.size());
    for (std::uint32_t i = 0; i < sortedTops.size(); ++i) {
      cotThetaTops[i] = state.linCircleTop[i].cotTheta;
    }
    std::ranges::sort(sortedTops, {},
                      [&](const std::uint32_t s) { return cotThetaTops[s]; });
  }

  // Reserve enough space, in case current capacity is too little
  state.topSpVec.reserve(numTopSp);
  state.curvatures.reserve(numTopSp);
  state.impactParameters.reserve(numTopSp);

  std::size_t t0 = 0;

  // clear previous results and then loop on bottoms and tops
  state.candidatesCollector.clear();

  for (const std::size_t b : sortedBottoms) {
    // break if we reached the last top SP
    if (t0 == numTopSp) {
      break;
    }

    auto lb = state.linCircleBottom[b];
    float cotThetaB = lb.cotTheta;
    float Vb = lb.V;
    float Ub = lb.U;
    float ErB = lb.Er;
    float iDeltaRB = lb.iDeltaR;

    // 1+(cot^2(theta)) = 1/sin^2(theta)
    float iSinTheta2 = (1. + cotThetaB * cotThetaB);
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
    state.topSpVec.clear();
    state.curvatures.clear();
    state.impactParameters.clear();

    // coordinate transformation and checks for middle spacepoint
    // x and y terms for the rotation from UV to XY plane
    float rotationTermsUVtoXY[2] = {0, 0};
    if constexpr (detailedMeasurement == DetectorMeasurementInfo::eDetailed) {
      rotationTermsUVtoXY[0] = cosPhiM * sinTheta;
      rotationTermsUVtoXY[1] = sinPhiM * sinTheta;
    }

    // minimum number of compatible top SPs to trigger the filter for a certain
    // middle bottom pair if seedConfirmation is false we always ask for at
    // least one compatible top to trigger the filter
    std::size_t minCompatibleTopSPs = 2;
    if (!m_config.seedConfirmation ||
        state.compatBottomSP[b]->radius() > seedFilterState.rMaxSeedConf) {
      minCompatibleTopSPs = 1;
    }
    if (m_config.seedConfirmation &&
        state.candidatesCollector.nHighQualityCandidates()) {
      minCompatibleTopSPs++;
    }

    for (std::size_t index_t = t0; index_t < numTopSp; index_t++) {
      const std::size_t t = sortedTops[index_t];

      auto lt = state.linCircleTop[t];

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

      if constexpr (detailedMeasurement == DetectorMeasurementInfo::eDetailed) {
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

        if (!xyzCoordinateCheck(m_config, spM, positionMiddle, rMTransf)) {
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

        auto spB = state.compatBottomSP[b];
        double rBTransf[3];
        if (!xyzCoordinateCheck(m_config, *spB, positionBottom, rBTransf)) {
          continue;
        }

        // coordinate transformation and checks for top spacepoint
        float Ct = 1. - B0 * lt.y;
        float St = A0 + B0 * lt.x;
        double positionTop[3] = {
            rotationTermsUVtoXY[0] * Ct - rotationTermsUVtoXY[1] * St,
            rotationTermsUVtoXY[0] * St + rotationTermsUVtoXY[1] * Ct,
            zPositionMiddle};

        auto spT = state.compatTopSP[t];
        double rTTransf[3];
        if (!xyzCoordinateCheck(m_config, *spT, positionTop, rTTransf)) {
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
      if constexpr (detailedMeasurement == DetectorMeasurementInfo::eDetailed) {
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
      if (deltaCotTheta2 > (error2 + scatteringInRegion2)) {
        // skip top SPs based on cotTheta sorting when producing triplets
        if constexpr (detailedMeasurement ==
                      DetectorMeasurementInfo::eDetailed) {
          continue;
        }
        // break if cotTheta from bottom SP < cotTheta from top SP because
        // the SP are sorted by cotTheta
        if (cotThetaB - cotThetaT < 0) {
          break;
        }
        t0 = index_t + 1;
        continue;
      }

      if constexpr (detailedMeasurement == DetectorMeasurementInfo::eDetailed) {
        rMxy = std::sqrt(rMTransf[0] * rMTransf[0] + rMTransf[1] * rMTransf[1]);
        float irMxy = 1 / rMxy;
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

      if constexpr (detailedMeasurement == DetectorMeasurementInfo::eDetailed) {
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
      if (!std::isinf(m_config.maxPtScattering)) {
        // if pT > maxPtScattering, calculate allowed scattering angle using
        // maxPtScattering instead of pt.
        // To avoid 0-divison the pT check is skipped in case of B2==0, and
        // p2scatterSigma is calculated directly from maxPtScattering
        if (B2 == 0 || options.pTPerHelixRadius * std::sqrt(S2 / B2) >
                           2. * m_config.maxPtScattering) {
          float pTscatterSigma =
              (m_config.highland / m_config.maxPtScattering) *
              m_config.sigmaScattering;
          p2scatterSigma = pTscatterSigma * pTscatterSigma * iSinTheta2;
        }
      }

      // if deltaTheta larger than allowed scattering for calculated pT, skip
      if (deltaCotTheta2 > (error2 + p2scatterSigma)) {
        if constexpr (detailedMeasurement ==
                      DetectorMeasurementInfo::eDetailed) {
          continue;
        }
        if (cotThetaB - cotThetaT < 0) {
          break;
        }
        t0 = index_t;
        continue;
      }
      // A and B allow calculation of impact params in U/V plane with linear
      // function
      // (in contrast to having to solve a quadratic function in x/y plane)
      float Im = 0;
      if constexpr (detailedMeasurement == DetectorMeasurementInfo::eDetailed) {
        Im = std::abs((A - B * rMxy) * rMxy);
      } else {
        Im = std::abs((A - B * rM) * rM);
      }

      if (Im > m_config.impactMax) {
        continue;
      }

      state.topSpVec.push_back(state.compatTopSP[t]);
      // inverse diameter is signed depending on if the curvature is
      // positive/negative in phi
      state.curvatures.push_back(B / std::sqrt(S2));
      state.impactParameters.push_back(Im);
    }  // loop on tops

    // continue if number of top SPs is smaller than minimum required for filter
    if (state.topSpVec.size() < minCompatibleTopSPs) {
      continue;
    }

    seedFilterState.zOrigin = spM.z() - rM * lb.cotTheta;

    m_config.seedFilter->filterSeeds_2SpFixed(
        state.spacePointMutableData, *state.compatBottomSP[b], spM,
        state.topSpVec, state.curvatures, state.impactParameters,
        seedFilterState, state.candidatesCollector);
  }  // loop on bottoms
}

template <typename external_spacepoint_t, typename grid_t, typename platform_t>
std::pair<float, float> SeedFinder<external_spacepoint_t, grid_t, platform_t>::
    retrieveRadiusRangeForMiddle(const external_spacepoint_t& spM,
                                 const Range1D<float>& rMiddleSPRange) const {
  if (m_config.useVariableMiddleSPRange) {
    return {rMiddleSPRange.min(), rMiddleSPRange.max()};
  }
  if (!m_config.rRangeMiddleSP.empty()) {
    /// get zBin position of the middle SP
    auto pVal = std::lower_bound(m_config.zBinEdges.begin(),
                                 m_config.zBinEdges.end(), spM.z());
    int zBin = std::distance(m_config.zBinEdges.begin(), pVal);
    /// protects against zM at the limit of zBinEdges
    zBin == 0 ? zBin : --zBin;
    return {m_config.rRangeMiddleSP[zBin][0], m_config.rRangeMiddleSP[zBin][1]};
  }
  return {m_config.rMinMiddle, m_config.rMaxMiddle};
}

}  // namespace Acts
