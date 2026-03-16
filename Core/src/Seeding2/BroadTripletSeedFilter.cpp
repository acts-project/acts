// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/BroadTripletSeedFilter.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding2/TripletSeedFinder.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <numeric>

namespace Acts {

namespace {

float getBestSeedQuality(
    const std::unordered_map<SpacePointIndex2, float>& bestSeedQualityMap,
    SpacePointIndex2 sp) {
  auto it = bestSeedQualityMap.find(sp);
  if (it != bestSeedQualityMap.end()) {
    return it->second;
  }
  return std::numeric_limits<float>::lowest();
}

void setBestSeedQuality(
    std::unordered_map<SpacePointIndex2, float>& bestSeedQualityMap,
    SpacePointIndex2 bottom, SpacePointIndex2 middle, SpacePointIndex2 top,
    float quality) {
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

BroadTripletSeedFilter::BroadTripletSeedFilter(const Config& config,
                                               State& state, Cache& cache,
                                               const Logger& logger)
    : m_cfg(&config), m_state(&state), m_cache(&cache), m_logger(&logger) {
  state.candidatesCollector =
      CandidatesForMiddleSp2(this->config().maxSeedsPerSpMConf,
                             this->config().maxQualitySeedsPerSpMConf);
}

bool BroadTripletSeedFilter::sufficientTopDoublets(
    const SpacePointContainer2& /*spacePoints*/,
    const ConstSpacePointProxy2& spM,
    const DoubletsForMiddleSp& topDoublets) const {
  // apply cut on the number of top SP if seedConfirmation is true
  if (!config().seedConfirmation) {
    return true;
  }

  // check if middle SP is in the central or forward region
  const bool isForwardRegion =
      spM.zr()[0] > config().centralSeedConfirmationRange.zMaxSeedConf ||
      spM.zr()[0] < config().centralSeedConfirmationRange.zMinSeedConf;
  SeedConfirmationRangeConfig seedConfRange =
      isForwardRegion ? config().forwardSeedConfirmationRange
                      : config().centralSeedConfirmationRange;
  // set the minimum number of top SP depending on whether the middle SP is
  // in the central or forward region
  std::size_t nTopSeedConf = spM.zr()[1] > seedConfRange.rMaxSeedConf
                                 ? seedConfRange.nTopForLargeR
                                 : seedConfRange.nTopForSmallR;
  // set max bottom radius for seed confirmation
  state().rMaxSeedConf = seedConfRange.rMaxSeedConf;
  // continue if number of top SPs is smaller than minimum
  if (topDoublets.size() < nTopSeedConf) {
    ACTS_VERBOSE("Number of top SPs is "
                 << topDoublets.size()
                 << " and is smaller than minimum, returning");
    return false;
  }

  return true;
}

void BroadTripletSeedFilter::filterTripletTopCandidates(
    const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
    const DoubletsForMiddleSp::Proxy& bottomLink,
    const TripletTopCandidates& tripletTopCandidates) const {
  auto spB = spacePoints[bottomLink.spacePointIndex()];
  auto bottomSp = spB.index();
  auto middleSp = spM.index();

  // minimum number of compatible top SPs to trigger the filter for a certain
  // middle bottom pair if seedConfirmation is false we always ask for at
  // least one compatible top to trigger the filter
  std::size_t minCompatibleTopSPs = 2;
  if (!config().seedConfirmation || spB.zr()[1] > state().rMaxSeedConf) {
    minCompatibleTopSPs = 1;
  }
  if (config().seedConfirmation &&
      state().candidatesCollector.nHighQualityCandidates() > 0) {
    minCompatibleTopSPs++;
  }
  // continue if number of top SPs is smaller than minimum required for filter
  if (tripletTopCandidates.size() < minCompatibleTopSPs) {
    return;
  }
  float zOrigin = spM.zr()[0] - spM.zr()[1] * bottomLink.cotTheta();

  // seed confirmation
  SeedConfirmationRangeConfig seedConfRange;
  std::size_t nTopSeedConf = 0;
  if (config().seedConfirmation) {
    // check if bottom SP is in the central or forward region
    const bool isForwardRegion =
        spB.zr()[0] > config().centralSeedConfirmationRange.zMaxSeedConf ||
        spB.zr()[0] < config().centralSeedConfirmationRange.zMinSeedConf;
    seedConfRange = isForwardRegion ? config().forwardSeedConfirmationRange
                                    : config().centralSeedConfirmationRange;
    // set the minimum number of top SP depending on whether the bottom SP is
    // in the central or forward region
    nTopSeedConf = spB.zr()[1] > seedConfRange.rMaxSeedConf
                       ? seedConfRange.nTopForLargeR
                       : seedConfRange.nTopForSmallR;
  }

  std::size_t maxWeightTopSp = 0;
  bool maxWeightSeed = false;
  float weightMax = std::numeric_limits<float>::lowest();

  // initialize original index locations
  cache().topSpIndexVec.resize(tripletTopCandidates.size());
  std::iota(cache().topSpIndexVec.begin(), cache().topSpIndexVec.end(), 0);
  std::ranges::sort(cache().topSpIndexVec, {},
                    [&tripletTopCandidates](const std::size_t t) {
                      return tripletTopCandidates.curvatures()[t];
                    });

  // vector containing the radius of all compatible seeds
  cache().compatibleSeedR.reserve(config().compatSeedLimit);

  const auto getTopR = [&](ConstSpacePointProxy2 spT) {
    if (config().useDeltaRinsteadOfTopRadius) {
      return fastHypot(spT.zr()[1] - spM.zr()[1], spT.zr()[0] - spM.zr()[0]);
    }
    return spT.zr()[1];
  };

  std::size_t beginCompTopIndex = 0;
  // loop over top SPs and other compatible top SP candidates
  for (const std::size_t topSpIndex : cache().topSpIndexVec) {
    auto topSp = tripletTopCandidates.topSpacePoints()[topSpIndex];
    auto spT = spacePoints[topSp];

    cache().compatibleSeedR.clear();

    float invHelixDiameter = tripletTopCandidates.curvatures()[topSpIndex];
    float lowerLimitCurv = invHelixDiameter - config().deltaInvHelixDiameter;
    float upperLimitCurv = invHelixDiameter + config().deltaInvHelixDiameter;
    float currentTopR = getTopR(spT);
    float impact = tripletTopCandidates.impactParameters()[topSpIndex];

    float weight = -impact * config().impactWeightFactor;

    // loop over compatible top SP candidates
    for (std::size_t variableCompTopIndex = beginCompTopIndex;
         variableCompTopIndex < cache().topSpIndexVec.size();
         variableCompTopIndex++) {
      std::size_t compatibleTopSpIndex =
          cache().topSpIndexVec[variableCompTopIndex];
      if (compatibleTopSpIndex == topSpIndex) {
        continue;
      }
      auto otherSpT = spacePoints[tripletTopCandidates
                                      .topSpacePoints()[compatibleTopSpIndex]];

      float otherTopR = getTopR(otherSpT);

      // curvature difference within limits?
      if (tripletTopCandidates.curvatures()[compatibleTopSpIndex] <
          lowerLimitCurv) {
        // the SPs are sorted in curvature so we skip unnecessary iterations
        beginCompTopIndex = variableCompTopIndex + 1;
        continue;
      }
      if (tripletTopCandidates.curvatures()[compatibleTopSpIndex] >
          upperLimitCurv) {
        // the SPs are sorted in curvature so we skip unnecessary iterations
        break;
      }
      // compared top SP should have at least deltaRMin distance
      float deltaR = currentTopR - otherTopR;
      if (std::abs(deltaR) < config().deltaRMin) {
        continue;
      }
      bool newCompSeed = true;
      for (const float previousDiameter : cache().compatibleSeedR) {
        // original ATLAS code uses higher min distance for 2nd found compatible
        // seed (20mm instead of 5mm)
        // add new compatible seed only if distance larger than rmin to all
        // other compatible seeds
        if (std::abs(previousDiameter - otherTopR) < config().deltaRMin) {
          newCompSeed = false;
          break;
        }
      }
      if (newCompSeed) {
        cache().compatibleSeedR.push_back(otherTopR);
        weight += config().compatSeedWeight;
      }
      if (cache().compatibleSeedR.size() >= config().compatSeedLimit) {
        break;
      }
    }

    if (config().experimentCuts != nullptr) {
      // add detector specific considerations on the seed weight
      weight += config().experimentCuts->seedWeight(spB, spM, spT);
      // discard seeds according to detector specific cuts (e.g.: weight)
      if (!config().experimentCuts->singleSeedCut(weight, spB, spM, spT)) {
        continue;
      }
    }

    // increment in seed weight if number of compatible seeds is larger than
    // numSeedIncrement
    if (cache().compatibleSeedR.size() > config().numSeedIncrement) {
      weight += config().seedWeightIncrement;
    }

    if (config().seedConfirmation) {
      // seed confirmation cuts - keep seeds if they have specific values of
      // impact parameter, z-origin and number of compatible seeds inside a
      // pre-defined range that also depends on the region of the detector (i.e.
      // forward or central region) defined by SeedConfirmationRange
      int deltaSeedConf = cache().compatibleSeedR.size() + 1 - nTopSeedConf;
      if (deltaSeedConf < 0 ||
          (state().candidatesCollector.nHighQualityCandidates() != 0 &&
           deltaSeedConf == 0)) {
        continue;
      }
      bool seedRangeCuts =
          spB.zr()[1] < seedConfRange.seedConfMinBottomRadius ||
          std::abs(zOrigin) > seedConfRange.seedConfMaxZOrigin;
      if (seedRangeCuts && deltaSeedConf == 0 &&
          impact > seedConfRange.minImpactSeedConf) {
        continue;
      }

      // term on the weight that depends on the value of zOrigin
      weight += -(std::abs(zOrigin) * config().zOriginWeightFactor) +
                config().compatSeedWeight;

      // skip a bad quality seed if any of its constituents has a weight larger
      // than the seed weight
      if (weight <
              getBestSeedQuality(state().bestSeedQualityMap, spB.index()) &&
          weight <
              getBestSeedQuality(state().bestSeedQualityMap, spM.index()) &&
          weight <
              getBestSeedQuality(state().bestSeedQualityMap, spT.index())) {
        continue;
      }

      if (deltaSeedConf > 0) {
        // if we have not yet reached our max number of quality seeds we add the
        // new seed to outCont

        // Internally, "push" will also check the max number of quality seeds
        // for a middle sp.
        // If this is reached, we remove the seed with the lowest weight.
        state().candidatesCollector.push(bottomSp, middleSp, topSp, weight,
                                         zOrigin, true);
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

      state().candidatesCollector.push(bottomSp, middleSp, topSp, weight,
                                       zOrigin, false);
    }
  }  // loop on tops

  // if no high quality seed was found for a certain middle+bottom SP pair,
  // lower quality seeds can be accepted
  if (config().seedConfirmation && maxWeightSeed &&
      state().candidatesCollector.nHighQualityCandidates() == 0) {
    // if we have not yet reached our max number of seeds we add the new seed to
    // outCont

    state().candidatesCollector.push(bottomSp, middleSp, maxWeightTopSp,
                                     weightMax, zOrigin, false);
  }
}

void BroadTripletSeedFilter::filterTripletsMiddleFixed(
    const SpacePointContainer2& spacePoints,
    SeedContainer2& outputCollection) const {
  const std::size_t numQualitySeeds =
      state().candidatesCollector.nHighQualityCandidates();

  cache().sortedCandidates.clear();
  state().candidatesCollector.toSortedCandidates(cache().sortedCandidates);
  std::span<TripletCandidate2> sortedCandidates = cache().sortedCandidates;

  if (config().experimentCuts != nullptr) {
    config().experimentCuts->cutPerMiddleSp(sortedCandidates);
  }

  std::size_t maxSeeds = sortedCandidates.size();

  if (maxSeeds > config().maxSeedsPerSpM) {
    maxSeeds = config().maxSeedsPerSpM + 1;
  }

  // default filter removes the last seeds if maximum amount exceeded
  // ordering by weight by filterSeeds_2SpFixed means these are the lowest
  // weight seeds
  std::size_t numTotalSeeds = 0;
  for (const auto& [bottom, middle, top, bestSeedQuality, zOrigin,
                    qualitySeed] : sortedCandidates) {
    // stop if we reach the maximum number of seeds
    if (numTotalSeeds >= maxSeeds) {
      break;
    }

    std::array<SpacePointIndex2, 3> triplet{spacePoints[bottom].index(),
                                            spacePoints[middle].index(),
                                            spacePoints[top].index()};

    if (config().seedConfirmation) {
      // continue if higher-quality seeds were found
      if (numQualitySeeds > 0 && !qualitySeed) {
        continue;
      }
      if (bestSeedQuality <
              getBestSeedQuality(state().bestSeedQualityMap, triplet[0]) &&
          bestSeedQuality <
              getBestSeedQuality(state().bestSeedQualityMap, triplet[1]) &&
          bestSeedQuality <
              getBestSeedQuality(state().bestSeedQualityMap, triplet[2])) {
        continue;
      }
    }

    // set quality of seed components
    setBestSeedQuality(state().bestSeedQualityMap, triplet[0], triplet[1],
                       triplet[2], bestSeedQuality);

    ACTS_VERBOSE("Adding seed: original indices=["
                 << triplet[0] << ", " << triplet[1] << ", " << triplet[2]
                 << "], internal indices=[" << bottom << ", " << middle << ", "
                 << top << "], quality=" << bestSeedQuality
                 << ", vertexZ=" << zOrigin);

    auto seed = outputCollection.createSeed();
    seed.assignSpacePointIndices(triplet);
    seed.vertexZ() = zOrigin;
    seed.quality() = bestSeedQuality;

    ++numTotalSeeds;
  }

  ACTS_VERBOSE("Identified " << numTotalSeeds << " seeds");
}

}  // namespace Acts
