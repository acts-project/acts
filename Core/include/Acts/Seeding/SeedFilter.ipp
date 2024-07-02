// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <numeric>
#include <utility>

namespace Acts {
// constructor
template <typename external_spacepoint_t>
SeedFilter<external_spacepoint_t>::SeedFilter(
    SeedFilterConfig config,
    IExperimentCuts<external_spacepoint_t>* expCuts /* = 0*/)
    : m_cfg(config), m_experimentCuts(expCuts) {
  if (!config.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFilterConfig not in ACTS internal units in SeedFilter");
  }
}
// function to filter seeds based on all seeds with same bottom- and
// middle-spacepoint.
// return vector must contain weight of each seed
template <typename external_spacepoint_t>
void SeedFilter<external_spacepoint_t>::filterSeeds_2SpFixed(
    Acts::SpacePointData& spacePointData,
    const InternalSpacePoint<external_spacepoint_t>& bottomSP,
    const InternalSpacePoint<external_spacepoint_t>& middleSP,
    const std::vector<const InternalSpacePoint<external_spacepoint_t>*>&
        topSpVec,
    const std::vector<float>& invHelixDiameterVec,
    const std::vector<float>& impactParametersVec,
    SeedFilterState& seedFilterState,
    CandidatesForMiddleSp<const InternalSpacePoint<external_spacepoint_t>>&
        candidates_collector) const {
  // seed confirmation
  SeedConfirmationRangeConfig seedConfRange;
  if (m_cfg.seedConfirmation) {
    // check if bottom SP is in the central or forward region
    seedConfRange =
        (bottomSP.z() > m_cfg.centralSeedConfirmationRange.zMaxSeedConf ||
         bottomSP.z() < m_cfg.centralSeedConfirmationRange.zMinSeedConf)
            ? m_cfg.forwardSeedConfirmationRange
            : m_cfg.centralSeedConfirmationRange;
    // set the minimum number of top SP depending on whether the bottom SP is
    // in the central or forward region
    seedFilterState.nTopSeedConf =
        bottomSP.radius() > seedConfRange.rMaxSeedConf
            ? seedConfRange.nTopForLargeR
            : seedConfRange.nTopForSmallR;
  }

  std::size_t maxWeightSeedIndex = 0;
  bool maxWeightSeed = false;
  float weightMax = std::numeric_limits<float>::lowest();
  float zOrigin = seedFilterState.zOrigin;

  // initialize original index locations
  std::vector<std::size_t> topSPIndexVec(topSpVec.size());
  for (std::size_t i(0); i < topSPIndexVec.size(); ++i) {
    topSPIndexVec[i] = i;
  }

  if (topSpVec.size() > 2) {
    // sort indexes based on comparing values in invHelixDiameterVec
    std::sort(
        topSPIndexVec.begin(), topSPIndexVec.end(),
        [&invHelixDiameterVec](const std::size_t i1, const std::size_t i2) {
          return invHelixDiameterVec[i1] < invHelixDiameterVec[i2];
        });
  }

  // vector containing the radius of all compatible seeds
  std::vector<float> compatibleSeedR;
  compatibleSeedR.reserve(m_cfg.compatSeedLimit);

  std::size_t beginCompTopIndex = 0;
  // loop over top SPs and other compatible top SP candidates
  for (const std::size_t topSPIndex : topSPIndexVec) {
    // if two compatible seeds with high distance in r are found, compatible
    // seeds span 5 layers
    compatibleSeedR.clear();

    float invHelixDiameter = invHelixDiameterVec[topSPIndex];
    float lowerLimitCurv = invHelixDiameter - m_cfg.deltaInvHelixDiameter;
    float upperLimitCurv = invHelixDiameter + m_cfg.deltaInvHelixDiameter;
    // use deltaR instead of top radius
    float currentTopR =
        m_cfg.useDeltaRorTopRadius
            ? spacePointData.deltaR(topSpVec[topSPIndex]->index())
            : topSpVec[topSPIndex]->radius();
    float impact = impactParametersVec[topSPIndex];

    float weight = -(impact * m_cfg.impactWeightFactor);

    // loop over compatible top SP candidates
    for (std::size_t variableCompTopIndex = beginCompTopIndex;
         variableCompTopIndex < topSPIndexVec.size(); variableCompTopIndex++) {
      std::size_t compatibleTopSPIndex = topSPIndexVec[variableCompTopIndex];
      if (compatibleTopSPIndex == topSPIndex) {
        continue;
      }

      float otherTopR =
          m_cfg.useDeltaRorTopRadius
              ? spacePointData.deltaR(topSpVec[compatibleTopSPIndex]->index())
              : topSpVec[compatibleTopSPIndex]->radius();

      // curvature difference within limits?
      if (invHelixDiameterVec[compatibleTopSPIndex] < lowerLimitCurv) {
        // the SPs are sorted in curvature so we skip unnecessary iterations
        beginCompTopIndex = variableCompTopIndex + 1;
        continue;
      }
      if (invHelixDiameterVec[compatibleTopSPIndex] > upperLimitCurv) {
        // the SPs are sorted in curvature so we skip unnecessary iterations
        break;
      }
      // compared top SP should have at least deltaRMin distance
      float deltaR = currentTopR - otherTopR;
      if (std::abs(deltaR) < m_cfg.deltaRMin) {
        continue;
      }
      bool newCompSeed = true;
      for (const float previousDiameter : compatibleSeedR) {
        // original ATLAS code uses higher min distance for 2nd found compatible
        // seed (20mm instead of 5mm)
        // add new compatible seed only if distance larger than rmin to all
        // other compatible seeds
        if (std::abs(previousDiameter - otherTopR) < m_cfg.deltaRMin) {
          newCompSeed = false;
          break;
        }
      }
      if (newCompSeed) {
        compatibleSeedR.push_back(otherTopR);
        weight += m_cfg.compatSeedWeight;
      }
      if (compatibleSeedR.size() >= m_cfg.compatSeedLimit) {
        break;
      }
    }

    if (m_experimentCuts != nullptr) {
      // add detector specific considerations on the seed weight
      weight += m_experimentCuts->seedWeight(bottomSP, middleSP,
                                             *topSpVec[topSPIndex]);
      // discard seeds according to detector specific cuts (e.g.: weight)
      if (!m_experimentCuts->singleSeedCut(weight, bottomSP, middleSP,
                                           *topSpVec[topSPIndex])) {
        continue;
      }
    }

    // increment in seed weight if number of compatible seeds is larger than
    // numSeedIncrement
    if (compatibleSeedR.size() > m_cfg.numSeedIncrement) {
      weight += m_cfg.seedWeightIncrement;
    }

    if (m_cfg.seedConfirmation) {
      // seed confirmation cuts - keep seeds if they have specific values of
      // impact parameter, z-origin and number of compatible seeds inside a
      // pre-defined range that also depends on the region of the detector (i.e.
      // forward or central region) defined by SeedConfirmationRange
      int deltaSeedConf =
          compatibleSeedR.size() + 1 - seedFilterState.nTopSeedConf;
      if (deltaSeedConf < 0 ||
          (seedFilterState.numQualitySeeds != 0 && deltaSeedConf == 0)) {
        continue;
      }
      bool seedRangeCuts =
          bottomSP.radius() < seedConfRange.seedConfMinBottomRadius ||
          std::abs(zOrigin) > seedConfRange.seedConfMaxZOrigin;
      if (seedRangeCuts && deltaSeedConf == 0 &&
          impact > seedConfRange.minImpactSeedConf) {
        continue;
      }

      // term on the weight that depends on the value of zOrigin
      weight += -(std::abs(zOrigin) * m_cfg.zOriginWeightFactor) +
                m_cfg.compatSeedWeight;

      // skip a bad quality seed if any of its constituents has a weight larger
      // than the seed weight
      if (weight < spacePointData.quality(bottomSP.index()) &&
          weight < spacePointData.quality(middleSP.index()) &&
          weight < spacePointData.quality(topSpVec[topSPIndex]->index())) {
        continue;
      }

      if (deltaSeedConf > 0) {
        // if we have not yet reached our max number of quality seeds we add the
        // new seed to outCont

        // Internally, "push" will also check the max number of quality seeds
        // for a middle sp.
        // If this is reached, we remove the seed with the lowest weight.
        candidates_collector.push(bottomSP, middleSP, *topSpVec[topSPIndex],
                                  weight, zOrigin, true);
        if (seedFilterState.numQualitySeeds < m_cfg.maxQualitySeedsPerSpMConf) {
          // fill high quality seed
          seedFilterState.numQualitySeeds++;
        }

      } else if (weight > weightMax) {
        // store weight and index of the best "lower quality" seed
        weightMax = weight;
        maxWeightSeedIndex = topSPIndex;
        maxWeightSeed = true;
      }
    } else {
      // keep the normal behavior without seed quality confirmation
      // if we have not yet reached our max number of seeds we add the new seed
      // to outCont

      candidates_collector.push(bottomSP, middleSP, *topSpVec[topSPIndex],
                                weight, zOrigin, false);
      if (seedFilterState.numSeeds < m_cfg.maxSeedsPerSpMConf) {
        // fill seed
        seedFilterState.numSeeds++;
      }
    }
  }  // loop on tops
  // if no high quality seed was found for a certain middle+bottom SP pair,
  // lower quality seeds can be accepted
  if (m_cfg.seedConfirmation && maxWeightSeed &&
      seedFilterState.numQualitySeeds == 0) {
    // if we have not yet reached our max number of seeds we add the new seed to
    // outCont

    candidates_collector.push(bottomSP, middleSP, *topSpVec[maxWeightSeedIndex],
                              weightMax, zOrigin, false);
    if (seedFilterState.numSeeds < m_cfg.maxSeedsPerSpMConf) {
      // fill seed
      seedFilterState.numSeeds++;
    }
  }
}

// after creating all seeds with a common middle space point, filter again

template <typename external_spacepoint_t>
void SeedFilter<external_spacepoint_t>::filterSeeds_1SpFixed(
    Acts::SpacePointData& spacePointData,
    CandidatesForMiddleSp<const InternalSpacePoint<external_spacepoint_t>>&
        candidates_collector,
    const std::size_t numQualitySeeds,
    GenericBackInserter<Seed<external_spacepoint_t>> outIt) const {
  // retrieve all candidates
  // this collection is already sorted
  // higher weights first
  auto extended_collection = candidates_collector.storage();
  filterSeeds_1SpFixed(spacePointData, extended_collection, numQualitySeeds,
                       outIt);
}

template <typename external_spacepoint_t>
void SeedFilter<external_spacepoint_t>::filterSeeds_1SpFixed(
    Acts::SpacePointData& spacePointData,
    std::vector<typename CandidatesForMiddleSp<
        const InternalSpacePoint<external_spacepoint_t>>::value_type>&
        candidates,
    const std::size_t numQualitySeeds,
    GenericBackInserter<Seed<external_spacepoint_t>> outIt) const {
  if (m_experimentCuts != nullptr) {
    candidates = m_experimentCuts->cutPerMiddleSP(std::move(candidates));
  }

  unsigned int maxSeeds = candidates.size();

  if (maxSeeds > m_cfg.maxSeedsPerSpM) {
    maxSeeds = m_cfg.maxSeedsPerSpM + 1;
  }

  // default filter removes the last seeds if maximum amount exceeded
  // ordering by weight by filterSeeds_2SpFixed means these are the lowest
  // weight seeds
  unsigned int numTotalSeeds = 0;
  for (const auto& [bottom, medium, top, bestSeedQuality, zOrigin,
                    qualitySeed] : candidates) {
    // stop if we reach the maximum number of seeds
    if (numTotalSeeds >= maxSeeds) {
      break;
    }

    if (m_cfg.seedConfirmation) {
      // continue if higher-quality seeds were found
      if (numQualitySeeds > 0 && !qualitySeed) {
        continue;
      }
      if (bestSeedQuality < spacePointData.quality(bottom->index()) &&
          bestSeedQuality < spacePointData.quality(medium->index()) &&
          bestSeedQuality < spacePointData.quality(top->index())) {
        continue;
      }
    }

    // set quality of seed components
    spacePointData.setQuality(bottom->index(), bestSeedQuality);
    spacePointData.setQuality(medium->index(), bestSeedQuality);
    spacePointData.setQuality(top->index(), bestSeedQuality);

    outIt = Seed<external_spacepoint_t>{bottom->sp(), medium->sp(), top->sp(),
                                        zOrigin, bestSeedQuality};
    ++numTotalSeeds;
  }
}

}  // namespace Acts
