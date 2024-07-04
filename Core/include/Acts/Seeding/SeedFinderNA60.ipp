// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <type_traits>

namespace Acts {

template <typename external_spacepoint_t, typename grid_t, typename platform_t>
SeedFinderNA60<external_spacepoint_t, grid_t, platform_t>::SeedFinderNA60(
    const Acts::SeedFinderConfigNA60<external_spacepoint_t>& config)
    : m_config(config) {
  if (not config.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFinderConfigNA60 not in ACTS internal units in SeedFinderNA60");
  }
  if (std::isnan(config.deltaYMaxTopSP)) {
    throw std::runtime_error("Value of deltaYMaxTopSP was not initialised");
  }
  if (std::isnan(config.deltaYMinTopSP)) {
    throw std::runtime_error("Value of deltaYMinTopSP was not initialised");
  }
  if (std::isnan(config.deltaYMaxBottomSP)) {
    throw std::runtime_error("Value of deltaYMaxBottomSP was not initialised");
  }
  if (std::isnan(config.deltaYMinBottomSP)) {
    throw std::runtime_error("Value of deltaYMinBottomSP was not initialised");
  }
}

template <typename external_spacepoint_t, typename grid_t, typename platform_t>
template <template <typename...> typename container_t, typename sp_range_t>
void SeedFinderNA60<external_spacepoint_t, grid_t, platform_t>::createSeedsForGroup(
    const Acts::SeedFinderOptionsNA60& options, SeedingState& state,
    const grid_t& grid,
    std::back_insert_iterator<container_t<Seed<external_spacepoint_t>>> outIt,
    const sp_range_t& bottomSPsIdx, const std::size_t middleSPsIdx,
    const sp_range_t& topSPsIdx,
    const Acts::Range1D<float>& yMiddleSPRange,
    float zTarget) const {
  if (not options.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFinderOptionsNA60 not in ACTS internal units in SeedFinderNA60");
  }

  // This is used for seed filtering later
  const std::size_t max_num_seeds_per_spm =
      m_config.seedFilter->getSeedFilterConfig().maxSeedsPerSpMConf;
  const std::size_t max_num_quality_seeds_per_spm =
      m_config.seedFilter->getSeedFilterConfig().maxQualitySeedsPerSpMConf;

  state.candidates_collector.setMaxElements(max_num_seeds_per_spm,
                                            max_num_quality_seeds_per_spm);

  // If there are no bottom or top bins, just return and waste no time
  if (bottomSPsIdx.size() == 0 or topSPsIdx.size() == 0) {
    return;
  }

  // Get the middle space point candidates
  const auto& middleSPs = grid.at(middleSPsIdx);

  // neighbours
  // clear previous results
  state.bottomNeighbours.clear();
  state.topNeighbours.clear();

  // Fill
  // bottoms
  for (const std::size_t idx : bottomSPsIdx) {
    state.bottomNeighbours.emplace_back(
        grid, idx, middleSPs.front()->y() - m_config.deltaYMaxBottomSP);
  }
  // tops
  for (const std::size_t idx : topSPsIdx) {
    state.topNeighbours.emplace_back(
        grid, idx, middleSPs.front()->y() + m_config.deltaYMinTopSP);
  }

  for (const auto& spM : middleSPs) {
    float yM = spM->y()+zTarget;
    float rM = spM->radius();

    // check if spM is outside our radial region of interest
    if (m_config.useVariableMiddleSPRange) {
      if (yM < yMiddleSPRange.min()) {
        continue;
      }
      if (yM > yMiddleSPRange.max()) {
        // break because SPs are sorted in r
        break;
      }
    } else if (not m_config.yRangeMiddleSP.empty()) {
      /// get zBin position of the middle SP
      auto pVal = std::lower_bound(m_config.zBinEdges.begin(),
                                   m_config.zBinEdges.end(), spM->z());
      int zBin = std::distance(m_config.zBinEdges.begin(), pVal);
      /// protects against zM at the limit of zBinEdges
      zBin == 0 ? zBin : --zBin;
      if (yM < m_config.yRangeMiddleSP[zBin][0]) {
        continue;
      }
      if (yM > m_config.yRangeMiddleSP[zBin][1]) {
        // break because SPs are sorted in r
        break;
      }
    } else {
      if (m_config.verbose)
        std::cout << "\nNA60+_SeedFinderNA60_middleSP_yM,yMinMiddle,yMaxMiddle: "
                  << yM << " " << m_config.yMinMiddle << " "
                  << m_config.yMaxMiddle << std::endl;
      if (yM < m_config.yMinMiddle) {
        continue;
      }
      if (yM > m_config.yMaxMiddle) {
        // break because SPs are sorted in y
        break;
      }
    }

    // remove middle SPs on the last layer since there would be no outer SPs to
    // complete a seed
    float zM = spM->z();
    if (zM < m_config.zOutermostLayers.first or
        zM > m_config.zOutermostLayers.second) {
      continue;
    }

    const float uIP = -1. / rM;
    const float cosPhiM = -spM->x() * uIP;
    const float sinPhiM = -spM->y() * uIP;
    const float uIP2 = uIP * uIP;

    // Iterate over middle-top dublets
    getCompatibleDoublets<Acts::SpacePointCandidateTypeNA60::eTop>(
        state.spacePointData, options, grid, state.topNeighbours, *spM.get(),
        state.linCircleTop, state.compatTopSP, m_config.deltaYMinTopSP,
        m_config.deltaYMaxTopSP, uIP, uIP2, cosPhiM, sinPhiM);

    // no top SP found -> try next spM
    if (state.compatTopSP.empty()) {
      continue;
    }

    // apply cut on the number of top SP if seedConfirmation is true
    SeedFilterStateNA60 seedFilterState;
    if (m_config.seedConfirmation) {
      // check if middle SP is in the central or forward region
      SeedConfirmationRangeConfig seedConfRange = m_config.seedConfirmationRange;

      // set the minimum number of top SP depending on whether the middle SP is
      // in the central or forward region
      seedFilterState.nTopSeedConf = zM > seedConfRange.zMaxSeedConf
                                         ? seedConfRange.nTopForLargeR
                                         : seedConfRange.nTopForSmallR;
      // We plane 2,3,(4?) for the Middle SP
      // if 1,2,3 we can require 2
      // if 2,3,4 we can require 1
      // however we risk to loose something due to inefficiencies or decays

      if (m_config.verbose)
        std::cout
            << "NA60+_in seedfinder confirmation seedFilterState.nTopSeedConf "
            << seedFilterState.nTopSeedConf << std::endl;

      // set max bottom radius for seed confirmation
      seedFilterState.yMaxSeedConf = seedConfRange.rMaxSeedConf;

      if (m_config.verbose)
        std::cout
            << "NA60+_in seedfinder confirmation seedFilterState.yMaxSeedConf "
            << seedFilterState.yMaxSeedConf << std::endl;

      // continue if number of top SPs is smaller than minimum
      if (state.compatTopSP.size() < seedFilterState.nTopSeedConf) {
        continue;
      }
    }

    // Iterate over middle-bottom dublets
    getCompatibleDoublets<Acts::SpacePointCandidateTypeNA60::eBottom>(
        state.spacePointData, options, grid, state.bottomNeighbours, *spM.get(),
        state.linCircleBottom, state.compatBottomSP, m_config.deltaYMinBottomSP,
        m_config.deltaYMaxBottomSP, uIP, uIP2, cosPhiM, sinPhiM);

    // no bottom SP found -> try next spM
    if (state.compatBottomSP.empty()) {
      continue;
    }

    // filter candidates
    if (m_config.useDetailedDoubleMeasurementInfo) {
      filterCandidates<Acts::DetectorMeasurementInfoNA60::eDetailed>(
          state.spacePointData, *spM.get(), options, seedFilterState, state);
    } else {
      filterCandidates<Acts::DetectorMeasurementInfoNA60::eDefault>(
          state.spacePointData, *spM.get(), options, seedFilterState, state);
    }

    m_config.seedFilter->filterSeeds_1SpFixed(
        state.spacePointData, state.candidates_collector,
        seedFilterState.numQualitySeeds, outIt);

  }  // loop on mediums
}

template <typename external_spacepoint_t, typename grid_t, typename platform_t>
template <Acts::SpacePointCandidateTypeNA60 candidateType, typename out_range_t>
inline void
SeedFinderNA60<external_spacepoint_t, grid_t, platform_t>::getCompatibleDoublets(
    Acts::SpacePointData& spacePointData,
    const Acts::SeedFinderOptionsNA60& options,
    const grid_t& grid,
    boost::container::small_vector<Acts::Neighbour<grid_t>,
                                Acts::detail::ipow(3, grid_t::DIM)>&
    otherSPsNeighbours,
    const InternalSpacePoint<external_spacepoint_t>& mediumSP,
    std::vector<LinCircle>& linCircleVec, out_range_t& outVec,
    const float deltaYMinSP, const float deltaYMaxSP, const float uIP,
    const float uIP2, const float cosPhiM, const float sinPhiM) const {
  float impactMax = m_config.impactMax;
  if constexpr (candidateType == Acts::SpacePointCandidateTypeNA60::eBottom) {
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

  float deltaY1 = 0.;
  float deltaZ = 0.;
  float deltaR = 0.;

  for (auto& otherSPCol : otherSPsNeighbours) {
    const auto& otherSPs = grid.at(otherSPCol.index);
    if (otherSPs.size() == 0) {
      continue;
    }

    // we make a copy of the iterator here since we need it to remain
    // the same in the Neighbour object
    auto min_itr = otherSPCol.itr;

    // find the first SP inside the radius region of interest and update
    // the iterator so we don't need to look at the other SPs again
    for (; min_itr != otherSPs.end(); ++min_itr) {
      const auto& otherSP = *min_itr;
      if (m_config.verbose)
        std::cout << "NA60+_MIDDLE " << xM << " " << yM << " " << zM
                  << std::endl;
      if constexpr (candidateType == Acts::SpacePointCandidateTypeNA60::eBottom)
        if (m_config.verbose)
          std::cout << "BOTTOM " << otherSP->x() << " " << otherSP->y() << " "
                    << otherSP->z() << std::endl;
      if constexpr (candidateType == Acts::SpacePointCandidateTypeNA60::eTop)
        if (m_config.verbose)
          std::cout << "TOP " << otherSP->x() << " " << otherSP->y() << " "
                    << otherSP->z() << std::endl;

      if constexpr (candidateType == Acts::SpacePointCandidateTypeNA60::eBottom) {
        if (m_config.verbose)
          std::cout << "NA60+_SeedFinderNA60_getCompatibleDoublets_Bottom_DeltaR_"
                       "deltaRMin,deltRMax= "
                    << rM - otherSP->radius() << " " << deltaYMinSP << " "
                    << deltaYMaxSP << std::endl;
        // if r-distance is too big, try next SP in bin
        if ((rM - otherSP->radius()) <= deltaYMaxSP) {
          break;
        }
      } else {
        if (m_config.verbose)
          std::cout << "NA60+_SeedFinderNA60_getCompatibleDoublets_Top_DeltaR_"
                       "deltaRMin,deltRMax= "
                    << otherSP->y() - yM << " " << deltaYMinSP << " "
                    << deltaYMaxSP << std::endl;
        // if r-distance is too small, try next SP in bin
        if ((otherSP->y() - yM) >= deltaYMinSP) {
          break;
        }
      }
    }
    // We update the iterator in the Neighbour object
    // that mean that we have changed the middle space point
    // and the lower bound has moved accordingly
    otherSPCol.itr = min_itr;

    for (; min_itr != otherSPs.end(); ++min_itr) {
      const auto& otherSP = *min_itr;

      if constexpr (candidateType == Acts::SpacePointCandidateTypeNA60::eBottom) {
        if (m_config.verbose)
          std::cout << "NA60+_---BOTTOM " << otherSP->x() << " " << otherSP->y()
                    << " " << otherSP->z() << std::endl;
        deltaY1 = (yM - otherSP->y());
        deltaR = (rM - otherSP->radius());
        
        if (m_config.verbose)
          std::cout << "NA60+_SeedFinderNA60_getCompatibleDoublets_Bottom_DeltaR_"
                       "deltRMin= "
                    << (yM - otherSP->y()) << " " << deltaYMinSP
                    << std::endl;

        // if r-distance is too small, try next SP in bin
        if (deltaY1 < deltaYMinSP) {
          break;
        }
      } else {
        deltaY1 = (otherSP->radius() - rM);
        deltaR = (otherSP->radius() - rM);
        if (m_config.verbose) {
          std::cout << "NA60+_---TOP " << otherSP->x() << " " << otherSP->y()
                    << " " << otherSP->z() << std::endl;
          std::cout
              << "NA60+_SeedFinderNA60_getCompatibleDoublets_Top_DeltaR_deltRMax= "
              << otherSP->radius() - rM << " " << deltaYMaxSP << std::endl;
        }
        // if r-distance is too big, try next SP in bin
        if (deltaY1 > deltaYMaxSP) {
          break;
        }
      }

      if constexpr (candidateType == Acts::SpacePointCandidateTypeNA60::eBottom) {
        deltaZ = (zM - otherSP->z());
      } else {
        deltaZ = (otherSP->z() - zM);
      }

      // the longitudinal impact parameter zOrigin is defined as (zM - rM *
      // cotTheta) where cotTheta is the ratio Z/R (forward angle) of space
      // point duplet but instead we calculate (zOrigin * deltaY) and multiply
      // collisionRegion by deltaY to avoid divisions
      const float zOriginTimesDeltaR = (zM * deltaR - rM * deltaZ);

      if (m_config.verbose) {
        if constexpr (candidateType == Acts::SpacePointCandidateTypeNA60::eBottom)
          std::cout << "NA60+_SeedFinderNA60_getCompatibleDoublets_BOTTOM_"
                       "zOriginTimesDeltaR/DeltaR_collRegionMin,collRegionMax= "
                    << zOriginTimesDeltaR / deltaR << " "
                    << m_config.collisionRegionMin << " "
                    << m_config.collisionRegionMax << std::endl;
        if constexpr (candidateType == Acts::SpacePointCandidateTypeNA60::eTop)
          std::cout << "NA60+_SeedFinderNA60_getCompatibleDoublets_TOP_"
                       "zOriginTimesDeltaR/DeltaR_collRegionMin,collRegionMax= "
                    << zOriginTimesDeltaR / deltaR << " "
                    << m_config.collisionRegionMin << " "
                    << m_config.collisionRegionMax << std::endl;
      }
      // check if duplet origin on z axis within collision region
      if (zOriginTimesDeltaR < m_config.collisionRegionMin * deltaR or
          zOriginTimesDeltaR > m_config.collisionRegionMax * deltaR) {
        continue;
      }

      // if interactionPointCut is false we apply z cuts before coordinate
      // transformation to avoid unnecessary calculations. If
      // interactionPointCut is true we apply the curvature cut first because it
      // is more frequent but requires the coordinate transformation
      if (not m_config.interactionPointCut) {
        if (m_config.verbose) {
          if constexpr (candidateType == Acts::SpacePointCandidateTypeNA60::eBottom)
            std::cout << "NA60+_SeedFinderNA60_getCompatibleDoublets_BOTTOM_deltaZ/"
                         "deltaY,cotThetaMin,cotThetaMax: "
                      << deltaZ / deltaR << " " << -m_config.cotThetaMax << " "
                      << m_config.cotThetaMax << std::endl;
          if constexpr (candidateType == Acts::SpacePointCandidateTypeNA60::eTop)
            std::cout << "NA60+_SeedFinderNA60_getCompatibleDoublets_TOP_deltaZ/"
                         "deltaY,cotThetaMin,cotThetaMax: "
                      << deltaZ / deltaR << " " << -m_config.cotThetaMax << " "
                      << m_config.cotThetaMax << std::endl;
        }
        // check if duplet cotTheta is within the region of interest
        // cotTheta is defined as (deltaZ / deltaY) but instead we multiply
        // cotThetaMax by deltaY to avoid division
        if (deltaZ > m_config.cotThetaMax * deltaR or
            deltaZ < -m_config.cotThetaMax * deltaR) {
          continue;
        }

        if (m_config.verbose) {
          if constexpr (candidateType == Acts::SpacePointCandidateTypeNA60::eBottom)
            std::cout << "NA60+_SeedFinderNA60_getCompatibleDoublets_BOTTOM_deltaZ,"
                         "deltaZMin,deltaZMax: "
                      << deltaZ << " " << -m_config.deltaZMax << " "
                      << m_config.deltaZMax << std::endl;
          if constexpr (candidateType == Acts::SpacePointCandidateTypeNA60::eTop)
            std::cout << "NA60+_SeedFinderNA60_getCompatibleDoublets_TOP_deltaZ,"
                         "deltaZMin,deltaZMax: "
                      << deltaZ << " " << -m_config.deltaZMax << " "
                      << m_config.deltaZMax << std::endl;
        }
        // if z-distance between SPs is within max and min values
        if (deltaZ > m_config.deltaZMax or deltaZ < -m_config.deltaZMax) {
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
        spacePointData.setDeltaR(otherSP->index(),
                                 std::sqrt(deltaR2 + (deltaZ * deltaZ)));
        outVec.push_back(otherSP.get());
        if (m_config.verbose) {
          std::cout << "NA60+_cotTheta in getCompatible = " << cotTheta
                    << std::endl;
        }
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

      // interactionPointCut == true we apply this cut first cuts before
      // coordinate transformation to avoid unnecessary calculations
      if (std::abs(rM * yNewFrame) <= impactMax * xNewFrame) {
        const float iDeltaR = std::sqrt(iDeltaR2);
        const float cotTheta = deltaZ * iDeltaR;

        const float Er =
            ((varianceZM + otherSP->varianceZ()) +
             (cotTheta * cotTheta) * (varianceRM + otherSP->varianceR())) *
            iDeltaR2;

        // fill output vectors
        linCircleVec.emplace_back(cotTheta, iDeltaR, Er, uT, vT, xNewFrame,
                                  yNewFrame);
        spacePointData.setDeltaR(otherSP->index(),
                                 std::sqrt(deltaR2 + (deltaZ * deltaZ)));
        outVec.emplace_back(otherSP.get());
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


      const float iDeltaR = std::sqrt(iDeltaR2);
      const float cotTheta = deltaZ * iDeltaR;

      const float Er =
          ((varianceZM + otherSP->varianceZ()) +
           (cotTheta * cotTheta) * (varianceRM + otherSP->varianceR())) *
          iDeltaR2;

      // fill output vectors
      linCircleVec.emplace_back(cotTheta, iDeltaR, Er, uT, vT, xNewFrame,
                                yNewFrame);
      spacePointData.setDeltaR(otherSP->index(),
                               std::sqrt(deltaR2 + (deltaZ * deltaZ)));
      outVec.emplace_back(otherSP.get());
    }
  }
}

template <typename external_spacepoint_t, typename grid_t, typename platform_t>
template <Acts::DetectorMeasurementInfoNA60 detailedMeasurement>
inline void SeedFinderNA60<external_spacepoint_t, grid_t, platform_t>::filterCandidates(
    Acts::SpacePointData& spacePointData,
    const InternalSpacePoint<external_spacepoint_t>& spM,
    const Acts::SeedFinderOptionsNA60& options, SeedFilterStateNA60& seedFilterState,
    SeedingState& state) const {
  float rM = spM.radius();
  float cosPhiM = spM.x() / rM;
  float sinPhiM = spM.y() / rM;
  float varianceRM = spM.varianceR();
  float varianceZM = spM.varianceZ();

  std::size_t numTopSP = state.compatTopSP.size();

  // sort: make index vector
  std::vector<std::size_t> sorted_bottoms(state.linCircleBottom.size());
  for (std::size_t i(0); i < sorted_bottoms.size(); ++i) {
    sorted_bottoms[i] = i;
  }

  std::vector<std::size_t> sorted_tops(state.linCircleTop.size());
  for (std::size_t i(0); i < sorted_tops.size(); ++i) {
    sorted_tops[i] = i;
  }

  if constexpr (detailedMeasurement ==
                Acts::DetectorMeasurementInfoNA60::eDefault) {
    std::sort(sorted_bottoms.begin(), sorted_bottoms.end(),
              [&state](const std::size_t a, const std::size_t b) -> bool {
                return state.linCircleBottom[a].cotTheta <
                       state.linCircleBottom[b].cotTheta;
              });

    std::sort(sorted_tops.begin(), sorted_tops.end(),
              [&state](const std::size_t a, const std::size_t b) -> bool {
                return state.linCircleTop[a].cotTheta <
                       state.linCircleTop[b].cotTheta;
              });
  }

  // Reserve enough space, in case current capacity is too little

  if (m_config.verbose) {
    std::cout << "NA60+_SeedFinderNA60_filterCandidates numTopSP= " << numTopSP
              << std::endl;
    std::cout << "\nNA60+_SeedFinderNA60_filterCandidates sorted_tops= "
              << sorted_tops.size()
              << " sorted_bottoms= " << sorted_bottoms.size() << std::endl;
  }
  state.topSpVec.reserve(numTopSP);
  state.curvatures.reserve(numTopSP);
  state.impactParameters.reserve(numTopSP);

  size_t t0 = 0;

  // clear previous results and then loop on bottoms and tops
  state.candidates_collector.clear();

  for (const std::size_t b : sorted_bottoms) {
    if (m_config.verbose)
      std::cout << "NA60+_---- filterCandidates - Loop on bottom" << std::endl;
    // break if we reached the last top SP
    if (t0 == numTopSP) {
      break;
    }

    auto lb = state.linCircleBottom[b];
    float cotThetaB = lb.cotTheta;
    float Vb = lb.V;
    float Ub = lb.U;
    float ErB = lb.Er;
    float iDeltaRB = lb.iDeltaR;

    if (m_config.verbose)
      std::cout << "NA60+_cotThetaB in SeedFinderNA60_filterCandidates = "
                << cotThetaB << std::endl;

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
    if constexpr (detailedMeasurement ==
                  Acts::DetectorMeasurementInfoNA60::eDetailed) {
      rotationTermsUVtoXY[0] = cosPhiM * sinTheta;
      rotationTermsUVtoXY[1] = sinPhiM * sinTheta;
    }

    // minimum number of compatible top SPs to trigger the filter for a certain
    // middle bottom pair if seedConfirmation is false we always ask for at
    // least one compatible top to trigger the filter

    size_t minCompatibleTopSPs = 2;
    if (!m_config.seedConfirmation or
        state.compatBottomSP[b]->radius() > seedFilterState.yMaxSeedConf) {
      if (m_config.verbose)
        std::cout << "NA60+_state.compatBottomSP[b]->radius()= "
                  << state.compatBottomSP[b]->radius()
                  << " seedFilterState.rMaxSeedConf= "
                  << seedFilterState.yMaxSeedConf << std::endl;
      minCompatibleTopSPs = 1;
    }
    if (m_config.seedConfirmation and seedFilterState.numQualitySeeds) {
      minCompatibleTopSPs++;
    }

    // added by me - test
    minCompatibleTopSPs = 1;

    for (size_t index_t = t0; index_t < numTopSP; index_t++) {
      if (m_config.verbose)
        std::cout << "NA60+_---- filterCandidates - Loop on top" << std::endl;
      const std::size_t t = sorted_tops[index_t];

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

      if (m_config.verbose)
        std::cout << "NA60+_cotThetaT in SeedFinderNA60_filterCandidates = "
                  << cotThetaT << std::endl;

      if constexpr (detailedMeasurement ==
                    Acts::DetectorMeasurementInfoNA60::eDetailed) {
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

        if (!xyzCoordinateCheck(spacePointData, m_config, spM, positionMiddle,
                                rMTransf)) {
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
        if (!xyzCoordinateCheck(spacePointData, m_config, *spB, positionBottom,
                                rBTransf)) {
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
        if (!xyzCoordinateCheck(spacePointData, m_config, *spT, positionTop,
                                rTTransf)) {
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

      if (m_config.verbose)
        std::cout << "NA60+_SeedFinderNA60_cotThetaAvg2= " << cotThetaAvg2
                  << " cotThetaB= " << cotThetaB << " cotThetaT= " << cotThetaT
                  << std::endl;
      if constexpr (detailedMeasurement ==
                    Acts::DetectorMeasurementInfoNA60::eDetailed) {
        // use arithmetic average
        float averageCotTheta = 0.5 * (cotThetaB + cotThetaT);
        cotThetaAvg2 = averageCotTheta * averageCotTheta;
      } else if (cotThetaAvg2 <= 0) {
        if (m_config.verbose)
          std::cout << "NA60+_Skip SP because cotThetaAvg2 <= 0 [SP have not a "
                       "correct order in z]"
                    << std::endl;
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

      if (m_config.verbose) {
        std::cout << "NA60+_varianceRM = " << varianceRM
                  << " varianceZM= " << varianceZM << std::endl;
        std::cout << "NA60+_cotThetaB= " << cotThetaB
                  << " cotThetaT= " << cotThetaT << std::endl;
        std::cout << "NA60+_SeedFinderNA60_filterCandidates_deltaCotTheta2,error2,"
                     "scatteringInRegion2: "
                  << deltaCotTheta2 << " " << error2 << " "
                  << scatteringInRegion2 << std::endl;
      }
      if (deltaCotTheta2 > (error2 + scatteringInRegion2)) {
        // skip top SPs based on cotTheta sorting when producing triplets
        if constexpr (detailedMeasurement ==
                      Acts::DetectorMeasurementInfoNA60::eDetailed) {
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

      if constexpr (detailedMeasurement ==
                    Acts::DetectorMeasurementInfoNA60::eDetailed) {
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

      if constexpr (detailedMeasurement ==
                    Acts::DetectorMeasurementInfoNA60::eDetailed) {
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

      if (m_config.verbose) {
        std::cout << "NA60+_SeedFinderNA60_filterCandidates_minHelixDiameter2: "
                  << S2 / B2 << " "
                  << options.minHelixDiameter2 * m_config.helixCut << std::endl;
        std::cout << "NA60+_S2/B2= " << S2 / B2
                  << " options.minHelixDiameter2=" << options.minHelixDiameter2
                  << " m_config.helixCut= " << m_config.helixCut << std::endl;
      }
      if (S2 < B2 * options.minHelixDiameter2 * m_config.helixCut) {
        continue;
      }
      // std::cout << "NA60+_Compatibility cut 2" << std::endl;

      // refinement of the cut on the compatibility between the r-z slope of
      // the two seed segments using a scattering term scaled by the actual
      // measured pT (p2scatterSigma)
      float iHelixDiameter2 = B2 / S2;
      // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
      // from rad to deltaCotTheta
      float p2scatterSigma = iHelixDiameter2 * sigmaSquaredPtDependent;

      if (m_config.verbose)
        std::cout << "NA60+_p2scatterSigma= " << p2scatterSigma << std::endl;
      if (!std::isinf(m_config.maxPtScattering)) {
        // if pT > maxPtScattering, calculate allowed scattering angle using
        // maxPtScattering instead of pt.
        // To avoid 0-divison the pT check is skipped in case of B2==0, and
        // p2scatterSigma is calculated directly from maxPtScattering
        if (B2 == 0 or options.pTPerHelixRadius * std::sqrt(S2 / B2) >
                           2. * m_config.maxPtScattering) {
          float pTscatterSigma =
              (m_config.highland / m_config.maxPtScattering) *
              m_config.sigmaScattering;
          p2scatterSigma = pTscatterSigma * pTscatterSigma * iSinTheta2;

          if (m_config.verbose)
            std::cout << "NA60+_inside p2scatterSigma= " << p2scatterSigma
                      << std::endl;
        }
      }

      // if deltaTheta larger than allowed scattering for calculated pT, skip

      if (m_config.verbose)
        std::cout << "NA60+_SeedFinderNA60_filterCandidates_deltaTheta_scattering,"
                     "error2,p2scatterSigma: "
                  << deltaCotTheta2 << " " << error2 << " " << p2scatterSigma
                  << std::endl;
      if (deltaCotTheta2 > (error2 + p2scatterSigma)) {
        if constexpr (detailedMeasurement ==
                      Acts::DetectorMeasurementInfoNA60::eDetailed) {
          continue;
        }
        if (cotThetaB - cotThetaT < 0) {
          break;
        }
        t0 = index_t;
        continue;
      }
      //      std::cout << "NA60+_Compatibility cut 3" << std::endl;

      // A and B allow calculation of impact params in U/V plane with linear
      // function
      // (in contrast to having to solve a quadratic function in x/y plane)
      float Im = 0;
      if constexpr (detailedMeasurement ==
                    Acts::DetectorMeasurementInfoNA60::eDetailed) {
        Im = std::abs((A - B * rMxy) * rMxy);
      } else {
        Im = std::abs((A - B * rM) * rM);
      }

      if (m_config.verbose)
        std::cout << "NA60+_SeedFinderNA60_filterCandidates_Im_impactMax: " << Im
                  << " " << m_config.impactMax << std::endl;

      if (Im > m_config.impactMax) {
        continue;
      }
      //      std::cout << "NA60+_Compatibility cut 4" << std::endl;

      state.topSpVec.push_back(state.compatTopSP[t]);
      // inverse diameter is signed depending on if the curvature is
      // positive/negative in phi
      state.curvatures.push_back(B / std::sqrt(S2));
      state.impactParameters.push_back(Im);

      // ADDED
      // store parameters of seeds surviving cuts (before filtering)

      auto Bzfield = 1.5;  // by hand
      if (m_config.verbose) {
        if ((cotThetaB < 0) && (cotThetaT < 0))
          std::cout << "NA60+_SeedParameter_cotThetaAvg= "
                    << -sqrt(cotThetaAvg2) << std::endl;
        else
          std::cout << "NA60+_SeedParameter_cotThetaAvg= " << sqrt(cotThetaAvg2)
                    << std::endl;
        std::cout << "NA60+_SeedParameter_pT= "
                  << 1. / (B / std::sqrt(S2)) * 300 * Bzfield / 2. / 1e6
                  << std::endl;  // from pdf Louis (why 1e6?
        std::cout << "NA60+_SeedParameter_d0= " << Im << std::endl;
        std::cout << "NA60+_SeedParameter_curvature= " << B / std::sqrt(S2)
                  << std::endl;
      }
    }  // loop on tops

    // continue if number of top SPs is smaller than minimum required for filter

    if (m_config.verbose)
      std::cout << "NA60+_state.topSpVec.size()= " << state.topSpVec.size()
                << " minCompatibleTopSPs= " << minCompatibleTopSPs << std::endl;
    if (state.topSpVec.size() < minCompatibleTopSPs) {
      continue;
    }

    seedFilterState.zOrigin = spM.z() - rM * lb.cotTheta;

    if (m_config.verbose)
      std::cout << "NA60+_SeedFilter_filterSeeds_2SpFixed" << std::endl;

    m_config.seedFilter->filterSeeds_2SpFixed(
        state.spacePointData, *state.compatBottomSP[b], spM, state.topSpVec,
        state.curvatures, state.impactParameters, seedFilterState,
        state.candidates_collector);
  }  // loop on bottoms
}

template <typename external_spacepoint_t, typename grid_t, typename platform_t>
template <typename sp_range_t>
std::vector<Seed<external_spacepoint_t>>
SeedFinderNA60<external_spacepoint_t, grid_t, platform_t>::createSeedsForGroup(
    const Acts::SeedFinderOptionsNA60& options,
    const grid_t& grid,
    const sp_range_t& bottomSPs, const std::size_t middleSPs,
    const sp_range_t& topSPs) const {
  SeedingState state;
  const Acts::Range1D<float> yMiddleSPRange;
  float zTarget;
  std::vector<Seed<external_spacepoint_t>> ret;

  createSeedsForGroup(options, state, grid, std::back_inserter(ret), bottomSPs,
                      middleSPs, topSPs, yMiddleSPRange, zTarget);

  return ret;
}
}  // namespace Acts
