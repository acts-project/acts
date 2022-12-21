// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
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
  if (not config.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFilterConfig not in ACTS internal units in SeedFilter");
  }
}
// function to filter seeds based on all seeds with same bottom- and
// middle-spacepoint.
// return vector must contain weight of each seed
template <typename external_spacepoint_t>
void SeedFilter<external_spacepoint_t>::filterSeeds_2SpFixed(
    InternalSpacePoint<external_spacepoint_t>& bottomSP,
    InternalSpacePoint<external_spacepoint_t>& middleSP,
    std::vector<InternalSpacePoint<external_spacepoint_t>*>& topSpVec,
    std::vector<float>& invHelixDiameterVec,
    std::vector<float>& impactParametersVec, SeedFilterState& seedFilterState,
    std::vector<std::pair<
        float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>&
        outCont) const {
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

  size_t maxWeightSeedIndex = 0;
  bool maxWeightSeed = false;
  float weightMax = -std::numeric_limits<float>::max();
  float zOrigin = seedFilterState.zOrigin;

  // initialize original index locations
  std::vector<size_t> topSPIndexVec(topSpVec.size());
  std::iota(topSPIndexVec.begin(), topSPIndexVec.end(), 0);

  if (m_cfg.curvatureSortingInFilter and topSpVec.size() > 2) {
    // sort indexes based on comparing values in invHelixDiameterVec
    std::sort(topSPIndexVec.begin(), topSPIndexVec.end(),
              [&invHelixDiameterVec](size_t i1, size_t i2) {
                return invHelixDiameterVec[i1] < invHelixDiameterVec[i2];
              });
  }

  size_t beginCompTopIndex = 0;
  // loop over top SPs and other compatible top SP candidates
  for (auto& topSPIndex : topSPIndexVec) {
    // if two compatible seeds with high distance in r are found, compatible
    // seeds span 5 layers
    // -> weaker requirement for a good seed
    std::vector<float> compatibleSeedR;

    float invHelixDiameter = invHelixDiameterVec[topSPIndex];
    float lowerLimitCurv = invHelixDiameter - m_cfg.deltaInvHelixDiameter;
    float upperLimitCurv = invHelixDiameter + m_cfg.deltaInvHelixDiameter;
    // use deltaR instead of top radius
    float currentTopR = m_cfg.useDeltaRorTopRadius
                            ? topSpVec[topSPIndex]->deltaR()
                            : topSpVec[topSPIndex]->radius();
    float impact = impactParametersVec[topSPIndex];

    float weight = -(impact * m_cfg.impactWeightFactor);

    for (size_t variableCompTopIndex = beginCompTopIndex;
         variableCompTopIndex < topSPIndexVec.size(); variableCompTopIndex++) {
      size_t compatibleTopSPIndex = topSPIndexVec[variableCompTopIndex];
      if (compatibleTopSPIndex == topSPIndex) {
        continue;
      }

      float otherTopR = m_cfg.useDeltaRorTopRadius
                            ? topSpVec[compatibleTopSPIndex]->deltaR()
                            : topSpVec[compatibleTopSPIndex]->radius();

      // curvature difference within limits?
      if (invHelixDiameterVec[compatibleTopSPIndex] < lowerLimitCurv) {
        // if SPs are sorted in curvature we skip unnecessary iterations
        if (m_cfg.curvatureSortingInFilter) {
          beginCompTopIndex = variableCompTopIndex + 1;
        }
        continue;
      }
      if (invHelixDiameterVec[compatibleTopSPIndex] > upperLimitCurv) {
        // if SPs are sorted in curvature we skip unnecessary iterations
        if (m_cfg.curvatureSortingInFilter) {
          break;
        }
        continue;
      }
      // compared top SP should have at least deltaRMin distance
      float deltaR = currentTopR - otherTopR;
      if (std::abs(deltaR) < m_cfg.deltaRMin) {
        continue;
      }
      bool newCompSeed = true;
      for (float previousDiameter : compatibleSeedR) {
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
          (seedFilterState.numQualitySeeds != 0 and deltaSeedConf == 0)) {
        continue;
      }
      bool seedRangeCuts =
          bottomSP.radius() < seedConfRange.seedConfMinBottomRadius ||
          std::abs(zOrigin) > seedConfRange.seedConfMaxZOrigin;
      if (seedRangeCuts and deltaSeedConf == 0 and
          impact > seedConfRange.minImpactSeedConf) {
        continue;
      }

      // term on the weight that depends on the value of zOrigin
      weight += -(std::abs(zOrigin) * m_cfg.zOriginWeightFactor) +
                m_cfg.compatSeedWeight;

      // skip a bad quality seed if any of its constituents has a weight larger
      // than the seed weight
      if (weight < bottomSP.quality() and weight < middleSP.quality() and
          weight < topSpVec[topSPIndex]->quality()) {
        continue;
      }

      if (deltaSeedConf > 0) {
        // if we have not yet reached our max number of quality seeds we add the
        // new seed to outCont
        if (seedFilterState.numQualitySeeds < m_cfg.maxQualitySeedsPerSpMConf) {
          // fill high quality seed
          seedFilterState.numQualitySeeds++;
          outCont.push_back(std::make_pair(
              weight,
              std::make_unique<const InternalSeed<external_spacepoint_t>>(
                  bottomSP, middleSP, *topSpVec[topSPIndex], zOrigin, true)));
        } else {
          // otherwise we check if there is a lower quality seed to remove
          checkReplaceSeeds(bottomSP, middleSP, *topSpVec[topSPIndex], zOrigin,
                            true, weight, outCont);
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
      if (seedFilterState.numSeeds < m_cfg.maxSeedsPerSpMConf) {
        // fill seed
        seedFilterState.numSeeds++;
        outCont.push_back(std::make_pair(
            weight,
            std::make_unique<const InternalSeed<external_spacepoint_t>>(
                bottomSP, middleSP, *topSpVec[topSPIndex], zOrigin, false)));
      } else {
        // otherwise we check if there is a lower quality seed to remove
        checkReplaceSeeds(bottomSP, middleSP, *topSpVec[topSPIndex], zOrigin,
                          false, weight, outCont);
      }
    }
  }
  // if no high quality seed was found for a certain middle+bottom SP pair,
  // lower quality seeds can be accepted
  if (m_cfg.seedConfirmation and maxWeightSeed and
      seedFilterState.numQualitySeeds == 0) {
    // if we have not yet reached our max number of seeds we add the new seed to
    // outCont
    if (seedFilterState.numSeeds < m_cfg.maxSeedsPerSpMConf) {
      // fill seed
      seedFilterState.numSeeds++;
      outCont.push_back(std::make_pair(
          weightMax,
          std::make_unique<const InternalSeed<external_spacepoint_t>>(
              bottomSP, middleSP, *topSpVec[maxWeightSeedIndex], zOrigin,
              false)));
    } else {
      // otherwise we check if there is a lower quality seed to remove
      checkReplaceSeeds(bottomSP, middleSP, *topSpVec[maxWeightSeedIndex],
                        zOrigin, false, weightMax, outCont);
    }
  }
}

// after creating all seeds with a common middle space point, filter again
template <typename external_spacepoint_t>
void SeedFilter<external_spacepoint_t>::filterSeeds_1SpFixed(
    std::vector<std::pair<
        float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>&
        seedsPerSpM,
    int& numQualitySeeds,
    std::back_insert_iterator<std::vector<Seed<external_spacepoint_t>>> outIt)
    const {
  // sort by weight and iterate only up to configured max number of seeds per
  // middle SP
  std::sort((seedsPerSpM.begin()), (seedsPerSpM.end()),
            [](const std::pair<float, std::unique_ptr<const Acts::InternalSeed<
                                          external_spacepoint_t>>>& i1,
               const std::pair<float, std::unique_ptr<const Acts::InternalSeed<
                                          external_spacepoint_t>>>& i2) {
              if (i1.first != i2.first) {
                return i1.first > i2.first;
              } else {
                // This is for the case when the weights from different seeds
                // are same. This makes cpu & cuda results same
                float seed1_sum = 0;
                float seed2_sum = 0;
                for (int i = 0; i < 3; i++) {
                  seed1_sum += pow(i1.second->sp[i]->y(), 2) +
                               pow(i1.second->sp[i]->z(), 2);
                  seed2_sum += pow(i2.second->sp[i]->y(), 2) +
                               pow(i2.second->sp[i]->z(), 2);
                }
                return seed1_sum > seed2_sum;
              }
            });
  if (m_experimentCuts != nullptr) {
    seedsPerSpM = m_experimentCuts->cutPerMiddleSP(std::move(seedsPerSpM));
  }
  unsigned int maxSeeds = seedsPerSpM.size();

  if (maxSeeds > m_cfg.maxSeedsPerSpM) {
    maxSeeds = m_cfg.maxSeedsPerSpM + 1;
  }
  auto itBegin = seedsPerSpM.begin();
  auto it = seedsPerSpM.begin();
  // default filter removes the last seeds if maximum amount exceeded
  // ordering by weight by filterSeeds_2SpFixed means these are the lowest
  // weight seeds
  unsigned int numTotalSeeds = 0;
  for (; it < itBegin + seedsPerSpM.size(); ++it) {
    // stop if we reach the maximum number of seeds
    if (numTotalSeeds >= maxSeeds) {
      break;
    }

    float bestSeedQuality = (*it).first;

    if (m_cfg.seedConfirmation) {
      // continue if higher-quality seeds were found
      if (numQualitySeeds > 0 and (*it).second->qualitySeed() == false) {
        continue;
      }
      if (bestSeedQuality < (*it).second->sp[0]->quality() and
          bestSeedQuality < (*it).second->sp[1]->quality() and
          bestSeedQuality < (*it).second->sp[2]->quality()) {
        continue;
      }
    }

    // set quality of seed components
    (*it).second->sp[0]->setQuality(bestSeedQuality);
    (*it).second->sp[1]->setQuality(bestSeedQuality);
    (*it).second->sp[2]->setQuality(bestSeedQuality);

    outIt = Seed<external_spacepoint_t>{
        (*it).second->sp[0]->sp(), (*it).second->sp[1]->sp(),
        (*it).second->sp[2]->sp(), (*it).second->z(), bestSeedQuality};
    numTotalSeeds += 1;
  }
}

template <typename external_spacepoint_t>
void SeedFilter<external_spacepoint_t>::checkReplaceSeeds(
    InternalSpacePoint<external_spacepoint_t>& bottomSP,
    InternalSpacePoint<external_spacepoint_t>& middleSP,
    InternalSpacePoint<external_spacepoint_t>& topSp, float zOrigin,
    bool isQualitySeed, float weight,
    std::vector<std::pair<
        float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>&
        outCont) const {
  // find the index of the seeds with qualitySeed() == isQualitySeed in outCont
  // and store in seed_indices
  std::vector<size_t> seed_indices;
  seed_indices.reserve(outCont.size());
  auto it = std::find_if(
      outCont.begin(), outCont.end(),
      [&](std::pair<float,
                    std::unique_ptr<const InternalSeed<external_spacepoint_t>>>&
              weight_seed) {
        return weight_seed.second->qualitySeed() == isQualitySeed;
      });
  while (it != outCont.end()) {
    seed_indices.emplace_back(std::distance(std::begin(outCont), it));
    it = std::find_if(
        std::next(it), outCont.end(),
        [&](std::pair<
            float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>&
                outContCheck) {
          return outContCheck.second->qualitySeed() == isQualitySeed;
        });
  }

  // find index of the seed with the minimum weight
  size_t index =
      *std::min_element(std::begin(seed_indices), std::end(seed_indices),
                        [&outCont](const size_t& a, const size_t& b) {
                          return outCont.at(a).first < outCont.at(b).first;
                        });
  // replace that seed with the new one if new one is better
  if (outCont.at(index).first < weight) {
    outCont.at(index) = std::make_pair(
        weight, std::make_unique<const InternalSeed<external_spacepoint_t>>(
                    bottomSP, middleSP, topSp, zOrigin, isQualitySeed));
  }
}

}  // namespace Acts
