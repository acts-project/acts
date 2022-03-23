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
    : m_cfg(config.toInternalUnits()), m_experimentCuts(expCuts) {}

// function to filter seeds based on all seeds with same bottom- and
// middle-spacepoint.
// return vector must contain weight of each seed
template <typename external_spacepoint_t>
void SeedFilter<external_spacepoint_t>::filterSeeds_2SpFixed(
    InternalSpacePoint<external_spacepoint_t>& bottomSP,
    InternalSpacePoint<external_spacepoint_t>& middleSP,
    std::vector<InternalSpacePoint<external_spacepoint_t>*>& topSpVec,
    std::vector<float>& invHelixDiameterVec,
    std::vector<float>& impactParametersVec, float zOrigin,
    int& numQualitySeeds,
    std::back_insert_iterator<std::vector<std::pair<
        float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>>
        outIt) const {
  // seed confirmation
  int nTopSeedConf = 0;
  if (m_cfg.seedConfirmation) {
    float rMaxSeedConfirmation =
        std::abs(bottomSP.z()) < m_cfg.centralSeedConfirmationRange.zMaxSeedConf
            ? m_cfg.centralSeedConfirmationRange.rMaxSeedConf
            : m_cfg.forwardSeedConfirmationRange.rMaxSeedConf;
    nTopSeedConf = rMaxSeedConfirmation.nTopForSmallR;
    if (bottomSP.radius() > rMaxSeedConfirmation)
      nTopSeedConf = rMaxSeedConfirmation.nTopForLargeR;
  }

  size_t minWeightSeedIndex = 0;
  bool minWeightSeed = false;
  float weightMin = std::numeric_limits<float>::min();

  // initialize original index locations
  std::vector<size_t> idx(topSpVec.size());
  std::iota(idx.begin(), idx.end(), 0);

  if (m_cfg.curvatureSortingInFilter) {
    // sort indexes based on comparing values in invHelixDiameterVec
    std::sort(idx.begin(), idx.end(),
              [&invHelixDiameterVec](size_t i1, size_t i2) {
                return invHelixDiameterVec[i1] < invHelixDiameterVec[i2];
              });
  }

  for (auto& i : idx) {
    // if two compatible seeds with high distance in r are found, compatible
    // seeds span 5 layers
    // -> very good seed
    std::vector<float> compatibleSeedR;

    float invHelixDiameter = invHelixDiameterVec[i];
    float lowerLimitCurv = invHelixDiameter - m_cfg.deltaInvHelixDiameter;
    float upperLimitCurv = invHelixDiameter + m_cfg.deltaInvHelixDiameter;
    float currentTop_r;
    float impact = impactParametersVec[i];

    if (m_cfg.useDeltaRTopRadius) {
      // use deltaR instead of top radius
      currentTop_r = topSpVec[i]->deltaR();
    } else {
      currentTop_r = topSpVec[i]->radius();
    }

    float weight = -(impact * m_cfg.impactWeightFactor);
    for (auto& j : idx) {
      if (i == j) {
        continue;
      }
      float otherTop_r;
      if (m_cfg.useDeltaRTopRadius) {
        otherTop_r = topSpVec[i]->deltaR();
      } else {
        otherTop_r = topSpVec[i]->radius();
      }
      // compared top SP should have at least deltaRMin distance
      float deltaR = currentTop_r - otherTop_r;
      if (std::abs(deltaR) < m_cfg.deltaRMin) {
        continue;
      }
      // curvature difference within limits?
      if (invHelixDiameterVec[j] < lowerLimitCurv) {
        continue;
      }
      if (invHelixDiameterVec[j] > upperLimitCurv) {
        if (m_cfg.curvatureSortingInFilter) {
          break;
        }
        continue;
      }
      bool newCompSeed = true;
      for (float previousDiameter : compatibleSeedR) {
        // original ATLAS code uses higher min distance for 2nd found compatible
        // seed (20mm instead of 5mm)
        // add new compatible seed only if distance larger than rmin to all
        // other compatible seeds
        if (std::abs(previousDiameter - otherTop_r) < m_cfg.deltaRMin) {
          newCompSeed = false;
          break;
        }
      }
      if (newCompSeed) {
        compatibleSeedR.push_back(otherTop_r);
        weight += m_cfg.compatSeedWeight;
      }
      if (compatibleSeedR.size() >= m_cfg.compatSeedLimit) {
        break;
      }
    }

    if (m_experimentCuts != nullptr) {
      // add detector specific considerations on the seed weight
      weight += m_experimentCuts->seedWeight(bottomSP, middleSP, *topSpVec[i]);
      // discard seeds according to detector specific cuts (e.g.: weight)
      if (!m_experimentCuts->singleSeedCut(weight, bottomSP, middleSP,
                                           *topSpVec[i])) {
        continue;
      }
    }

    int deltaSeedConf;
    if (m_cfg.seedConfirmation) {
      // seed confirmation cuts
      deltaSeedConf = compatibleSeedR.size() + 1 - nTopSeedConf;
      if (deltaSeedConf < 0 || (numQualitySeeds and !deltaSeedConf)) {
        continue;
      }
      bool seedConfMinRange =
          bottomSP.radius() < m_cfg.seedConfMinBottomRadius ||
          std::abs(zOrigin) > m_cfg.seedConfMaxZOrigin;
      if (seedConfMinRange and !deltaSeedConf and
          impact > m_cfg.minImpactSeedConf) {
        continue;
      }

      // term on the weight that depends on the value of zOrigin
      weight += -std::abs(zOrigin) + m_cfg.compatSeedWeight;

      // skip a bad quality if any of the components has a weight smaller than
      // the seed weight
      if (weight < bottomSP.quality() and weight < middleSP.quality() and
          weight < topSpVec[i]->quality()) {
        continue;
      }

      if (deltaSeedConf) {
        // fill high quality seed
        ++numQualitySeeds;
        outIt = std::make_pair(
            weight, std::make_unique<const InternalSeed<external_spacepoint_t>>(
                        bottomSP, middleSP, *topSpVec[i], zOrigin));
      } else if (weight > weightMin) {
        // store index of lower quality seeds with a weight greater than the
        // minimum
        weightMin = weight;
        minWeightSeedIndex = i;
        minWeightSeed = true;
      }
    } else {
      // keep the normal behavior without seed quality confirmation
      outIt = std::make_pair(
          weight, std::make_unique<const InternalSeed<external_spacepoint_t>>(
                      bottomSP, middleSP, *topSpVec[i], zOrigin));
    }
  }
  // if no high quality seed was found for a certain middle+bottom SP pair,
  // lower quality seeds can be accepted
  if (m_cfg.seedConfirmation and minWeightSeed and !numQualitySeeds) {
    outIt = std::make_pair(
        weightMin,
        std::make_unique<const InternalSeed<external_spacepoint_t>>(
            bottomSP, middleSP, *topSpVec[minWeightSeedIndex], zOrigin));
  }
}

// after creating all seeds with a common middle space point, filter again
template <typename external_spacepoint_t>
void SeedFilter<external_spacepoint_t>::filterSeeds_1SpFixed(
    std::vector<std::pair<
        float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>&
        seedsPerSpM,
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
                  seed1_sum += pow(i1.second->sp[i]->sp().y(), 2) +
                               pow(i1.second->sp[i]->sp().z(), 2);
                  seed2_sum += pow(i2.second->sp[i]->sp().y(), 2) +
                               pow(i2.second->sp[i]->sp().z(), 2);
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
  for (; it < itBegin + maxSeeds; ++it) {
    float bestSeedQuality = (*it).first;

    (*it).second->sp[0]->setQuality(bestSeedQuality);
    (*it).second->sp[1]->setQuality(bestSeedQuality);
    (*it).second->sp[2]->setQuality(bestSeedQuality);

    outIt = Seed<external_spacepoint_t>{
        (*it).second->sp[0]->sp(), (*it).second->sp[1]->sp(),
        (*it).second->sp[2]->sp(), (*it).second->z()};
  }
}

}  // namespace Acts
