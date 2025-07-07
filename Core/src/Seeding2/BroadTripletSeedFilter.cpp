// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/BroadTripletSeedFilter.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding2/ITripletSeedCuts.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <numeric>

namespace Acts::Experimental {

namespace {

float getBestSeedQuality(
    const SpacePointContainer2& spacePoints,
    const std::unordered_map<SpacePointIndex2, float>& bestSeedQualityMap,
    SpacePointIndex2 sp) {
  if (spacePoints.hasExtraColumns(SpacePointKnownExtraColumn::CopyFromIndex)) {
    sp = spacePoints[sp].extra(spacePoints.copyFromIndexColumn());
  }
  auto it = bestSeedQualityMap.find(sp);
  if (it != bestSeedQualityMap.end()) {
    return it->second;
  }
  return std::numeric_limits<float>::lowest();
}

void setBestSeedQuality(
    const SpacePointContainer2& spacePoints,
    std::unordered_map<SpacePointIndex2, float>& bestSeedQualityMap,
    SpacePointIndex2 top, SpacePointIndex2 middle, SpacePointIndex2 bottom,
    float quality) {
  if (spacePoints.hasExtraColumns(SpacePointKnownExtraColumn::CopyFromIndex)) {
    top = spacePoints[top].copyFromIndex();
    middle = spacePoints[middle].copyFromIndex();
    bottom = spacePoints[bottom].copyFromIndex();
  }

  const auto set = [&](SpacePointIndex2 sp) {
    auto it = bestSeedQualityMap.find(sp);
    if (it != bestSeedQualityMap.end()) {
      it->second = std::max(quality, it->second);
    } else {
      bestSeedQualityMap.emplace(sp, quality);
    }
  };

  for (SpacePointIndex2 sp : {top, middle, bottom}) {
    set(sp);
  }
}

}  // namespace

BroadTripletSeedFilter::BroadTripletSeedFilter(
    const Config& config, std::unique_ptr<const Logger> logger)
    : m_cfg(config), m_logger(std::move(logger)) {}

void BroadTripletSeedFilter::filter2SpFixed(
    State& state, Cache& cache, const SpacePointContainer2& spacePoints,
    SpacePointIndex2 bottomSp, SpacePointIndex2 middleSp,
    std::span<const SpacePointIndex2> topSpVec,
    std::span<const float> invHelixDiameterVec,
    std::span<const float> impactParametersVec, float zOrigin,
    CandidatesForMiddleSp2& candidatesCollector) const {
  auto spB = spacePoints[bottomSp];
  auto spM = spacePoints[middleSp];

  // seed confirmation
  SeedConfirmationRangeConfig seedConfRange;
  std::size_t nTopSeedConf = 0;
  if (m_cfg.seedConfirmation) {
    // check if bottom SP is in the central or forward region
    const bool isForwardRegion =
        spB.z() > m_cfg.centralSeedConfirmationRange.zMaxSeedConf ||
        spB.z() < m_cfg.centralSeedConfirmationRange.zMinSeedConf;
    seedConfRange = isForwardRegion ? m_cfg.forwardSeedConfirmationRange
                                    : m_cfg.centralSeedConfirmationRange;
    // set the minimum number of top SP depending on whether the bottom SP is
    // in the central or forward region
    nTopSeedConf = spB.r() > seedConfRange.rMaxSeedConf
                       ? seedConfRange.nTopForLargeR
                       : seedConfRange.nTopForSmallR;
  }

  std::size_t maxWeightTopSp = 0;
  bool maxWeightSeed = false;
  float weightMax = std::numeric_limits<float>::lowest();

  // initialize original index locations
  cache.topSpIndexVec.resize(topSpVec.size());
  std::iota(cache.topSpIndexVec.begin(), cache.topSpIndexVec.end(), 0);
  std::ranges::sort(cache.topSpIndexVec, {},
                    [&invHelixDiameterVec](const std::size_t t) {
                      return invHelixDiameterVec[t];
                    });

  // vector containing the radius of all compatible seeds
  cache.compatibleSeedR.reserve(m_cfg.compatSeedLimit);

  const auto getTopR = [&](ConstSpacePointProxy2 spT) {
    if (m_cfg.useDeltaRinsteadOfTopRadius) {
      return fastHypot(spT.r() - spM.r(), spT.z() - spM.z());
    }
    return spT.r();
  };

  std::size_t beginCompTopIndex = 0;
  // loop over top SPs and other compatible top SP candidates
  for (const std::size_t topSpIndex : cache.topSpIndexVec) {
    auto topSp = topSpVec[topSpIndex];
    auto spT = spacePoints[topSp];

    cache.compatibleSeedR.clear();

    float invHelixDiameter = invHelixDiameterVec[topSpIndex];
    float lowerLimitCurv = invHelixDiameter - m_cfg.deltaInvHelixDiameter;
    float upperLimitCurv = invHelixDiameter + m_cfg.deltaInvHelixDiameter;
    float currentTopR = getTopR(spT);
    float impact = impactParametersVec[topSpIndex];

    float weight = -impact * m_cfg.impactWeightFactor;

    // loop over compatible top SP candidates
    for (std::size_t variableCompTopIndex = beginCompTopIndex;
         variableCompTopIndex < cache.topSpIndexVec.size();
         variableCompTopIndex++) {
      std::size_t compatibletopSpIndex =
          cache.topSpIndexVec[variableCompTopIndex];
      if (compatibletopSpIndex == topSpIndex) {
        continue;
      }
      auto otherSpT = spacePoints[topSpVec[compatibletopSpIndex]];

      float otherTopR = getTopR(otherSpT);

      // curvature difference within limits?
      if (invHelixDiameterVec[compatibletopSpIndex] < lowerLimitCurv) {
        // the SPs are sorted in curvature so we skip unnecessary iterations
        beginCompTopIndex = variableCompTopIndex + 1;
        continue;
      }
      if (invHelixDiameterVec[compatibletopSpIndex] > upperLimitCurv) {
        // the SPs are sorted in curvature so we skip unnecessary iterations
        break;
      }
      // compared top SP should have at least deltaRMin distance
      float deltaR = currentTopR - otherTopR;
      if (std::abs(deltaR) < m_cfg.deltaRMin) {
        continue;
      }
      bool newCompSeed = true;
      for (const float previousDiameter : cache.compatibleSeedR) {
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
        cache.compatibleSeedR.push_back(otherTopR);
        weight += m_cfg.compatSeedWeight;
      }
      if (cache.compatibleSeedR.size() >= m_cfg.compatSeedLimit) {
        break;
      }
    }

    if (m_cfg.experimentCuts != nullptr) {
      // add detector specific considerations on the seed weight
      weight += m_cfg.experimentCuts->seedWeight(spB, spM, spT);
      // discard seeds according to detector specific cuts (e.g.: weight)
      if (!m_cfg.experimentCuts->singleSeedCut(weight, spB, spM, spT)) {
        continue;
      }
    }

    // increment in seed weight if number of compatible seeds is larger than
    // numSeedIncrement
    if (cache.compatibleSeedR.size() > m_cfg.numSeedIncrement) {
      weight += m_cfg.seedWeightIncrement;
    }

    if (m_cfg.seedConfirmation) {
      // seed confirmation cuts - keep seeds if they have specific values of
      // impact parameter, z-origin and number of compatible seeds inside a
      // pre-defined range that also depends on the region of the detector (i.e.
      // forward or central region) defined by SeedConfirmationRange
      int deltaSeedConf = cache.compatibleSeedR.size() + 1 - nTopSeedConf;
      if (deltaSeedConf < 0 ||
          (candidatesCollector.nHighQualityCandidates() != 0 &&
           deltaSeedConf == 0)) {
        continue;
      }
      bool seedRangeCuts = spB.r() < seedConfRange.seedConfMinBottomRadius ||
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
      if (weight < getBestSeedQuality(spacePoints, state.bestSeedQualityMap,
                                      bottomSp) &&
          weight < getBestSeedQuality(spacePoints, state.bestSeedQualityMap,
                                      middleSp) &&
          weight < getBestSeedQuality(spacePoints, state.bestSeedQualityMap,
                                      topSp)) {
        continue;
      }

      if (deltaSeedConf > 0) {
        // if we have not yet reached our max number of quality seeds we add the
        // new seed to outCont

        // Internally, "push" will also check the max number of quality seeds
        // for a middle sp.
        // If this is reached, we remove the seed with the lowest weight.
        candidatesCollector.push(bottomSp, middleSp, topSp, weight, zOrigin,
                                 true);
      } else if (weight > weightMax) {
        // store weight and index of the best "lower quality" seed
        weightMax = weight;
        maxWeightTopSp = topSp;
        maxWeightSeed = true;
      }
    } else {
      // keep the normal behavior without seed quality confirmation
      // if we have not yet reached our max number of seeds we add the new seed
      // to outCont

      candidatesCollector.push(bottomSp, middleSp, topSp, weight, zOrigin,
                               false);
    }
  }  // loop on tops

  // if no high quality seed was found for a certain middle+bottom SP pair,
  // lower quality seeds can be accepted
  if (m_cfg.seedConfirmation && maxWeightSeed &&
      candidatesCollector.nHighQualityCandidates() == 0) {
    // if we have not yet reached our max number of seeds we add the new seed to
    // outCont

    candidatesCollector.push(bottomSp, middleSp, maxWeightTopSp, weightMax,
                             zOrigin, false);
  }
}

void BroadTripletSeedFilter::filter1SpFixed(
    State& state, const SpacePointContainer2& spacePoints,
    std::span<TripletCandidate2> candidates, std::size_t numQualitySeeds,
    SeedContainer2& outputCollection) const {
  if (m_cfg.experimentCuts != nullptr) {
    m_cfg.experimentCuts->cutPerMiddleSp(candidates);
  }

  std::size_t maxSeeds = candidates.size();

  if (maxSeeds > m_cfg.maxSeedsPerSpM) {
    maxSeeds = m_cfg.maxSeedsPerSpM + 1;
  }

  // default filter removes the last seeds if maximum amount exceeded
  // ordering by weight by filterSeeds_2SpFixed means these are the lowest
  // weight seeds
  std::size_t numTotalSeeds = 0;
  for (const auto& [bottom, middle, top, bestSeedQuality, zOrigin,
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
      if (bestSeedQuality < getBestSeedQuality(spacePoints,
                                               state.bestSeedQualityMap,
                                               bottom) &&
          bestSeedQuality < getBestSeedQuality(spacePoints,
                                               state.bestSeedQualityMap,
                                               middle) &&
          bestSeedQuality <
              getBestSeedQuality(spacePoints, state.bestSeedQualityMap, top)) {
        continue;
      }
    }

    // set quality of seed components
    setBestSeedQuality(spacePoints, state.bestSeedQualityMap, top, middle,
                       bottom, bestSeedQuality);

    ACTS_VERBOSE("Adding seed: [b=" << bottom << ", m=" << middle << ", t="
                                    << top << "], quality=" << bestSeedQuality
                                    << ", vertexZ=" << zOrigin);

    auto seed = outputCollection.createSeed(
        std::array<SpacePointIndex2, 3>{bottom, middle, top});
    seed.vertexZ() = zOrigin;
    seed.quality() = bestSeedQuality;

    ++numTotalSeeds;
  }

  ACTS_VERBOSE("Identified " << numTotalSeeds << " seeds");
}

}  // namespace Acts::Experimental
