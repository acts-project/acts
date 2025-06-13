// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/TripletSeedFilter2.hpp"

#include <numeric>

namespace Acts {

using namespace UnitLiterals;

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

}  // namespace

TripletSeedFilter2::DerivedConfig TripletSeedFilter2::Config::derive() const {
  DerivedConfig result;

  static_cast<Config&>(result) = *this;

  // TODO get rid of unit conversions
  {
    result.deltaRMin /= 1_mm;
    result.deltaInvHelixDiameter /= 1. / 1_mm;
  }

  return result;
}

TripletSeedFilter2::TripletSeedFilter2(const DerivedConfig& config,
                                       std::unique_ptr<const Logger> logger)
    : m_cfg(config), m_logger(std::move(logger)) {}

void TripletSeedFilter2::filter2SpFixed(
    const Options& options, State& state,
    const SpacePointContainer2& spacePoints,
    const SpacePointContainer2::DenseColumn<float>& rColumn,
    SpacePointIndex2 bottomSp, SpacePointIndex2 middleSp,
    const std::vector<SpacePointIndex2>& topSpVec,
    const std::vector<float>& invHelixDiameterVec,
    const std::vector<float>& impactParametersVec, float zOrigin,
    CandidatesForMiddleSp2& candidatesCollector) const {
  auto spB = spacePoints.at(bottomSp);
  auto spM = spacePoints.at(middleSp);

  std::size_t maxWeightTopSp = 0;
  bool maxWeightSeed = false;
  float weightMax = std::numeric_limits<float>::lowest();

  // initialize original index locations
  state.topSpIndexVec.resize(topSpVec.size());
  std::iota(state.topSpIndexVec.begin(), state.topSpIndexVec.end(), 0);
  std::ranges::sort(state.topSpIndexVec, {},
                    [&invHelixDiameterVec](const std::size_t t) {
                      return invHelixDiameterVec[t];
                    });

  // vector containing the radius of all compatible seeds
  state.compatibleSeedR.reserve(m_cfg.compatSeedLimit);

  std::size_t beginCompTopIndex = 0;
  // loop over top SPs and other compatible top SP candidates
  for (const std::size_t topSpIndex : state.topSpIndexVec) {
    auto topSp = topSpVec[topSpIndex];
    auto spT = spacePoints.at(topSp);

    // if two compatible seeds with high distance in r are found, compatible
    // seeds span 5 layers
    state.compatibleSeedR.clear();

    float invHelixDiameter = invHelixDiameterVec[topSpIndex];
    float lowerLimitCurv = invHelixDiameter - m_cfg.deltaInvHelixDiameter;
    float upperLimitCurv = invHelixDiameter + m_cfg.deltaInvHelixDiameter;
    // use deltaR instead of top radius
    assert(!m_cfg.useDeltaRorTopRadius && "TODO not implemented");
    float currentTopR = spT.extra(rColumn);
    float impact = impactParametersVec[topSpIndex];

    float weight = -(impact * m_cfg.impactWeightFactor);

    // loop over compatible top SP candidates
    for (std::size_t variableCompTopIndex = beginCompTopIndex;
         variableCompTopIndex < state.topSpIndexVec.size();
         variableCompTopIndex++) {
      std::size_t compatibletopSpIndex =
          state.topSpIndexVec[variableCompTopIndex];
      if (compatibletopSpIndex == topSpIndex) {
        continue;
      }
      auto otherSpT = spacePoints.at(topSpVec[compatibletopSpIndex]);

      assert(!m_cfg.useDeltaRorTopRadius && "TODO not implemented");
      float otherTopR = otherSpT.extra(rColumn);

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
      for (const float previousDiameter : state.compatibleSeedR) {
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
        state.compatibleSeedR.push_back(otherTopR);
        weight += m_cfg.compatSeedWeight;
      }
      if (state.compatibleSeedR.size() >= m_cfg.compatSeedLimit) {
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
    if (state.compatibleSeedR.size() > m_cfg.numSeedIncrement) {
      weight += m_cfg.seedWeightIncrement;
    }

    if (options.seedConfirmation) {
      // seed confirmation cuts - keep seeds if they have specific values of
      // impact parameter, z-origin and number of compatible seeds inside a
      // pre-defined range that also depends on the region of the detector (i.e.
      // forward or central region) defined by SeedConfirmationRange
      int deltaSeedConf =
          state.compatibleSeedR.size() + 1 - options.nTopSeedConf;
      if (deltaSeedConf < 0 ||
          (candidatesCollector.nHighQualityCandidates() != 0 &&
           deltaSeedConf == 0)) {
        continue;
      }
      bool seedRangeCuts =
          spB.extra(rColumn) < options.seedConfRange.seedConfMinBottomRadius ||
          std::abs(zOrigin) > options.seedConfRange.seedConfMaxZOrigin;
      if (seedRangeCuts && deltaSeedConf == 0 &&
          impact > options.seedConfRange.minImpactSeedConf) {
        continue;
      }

      // term on the weight that depends on the value of zOrigin
      weight += -(std::abs(zOrigin) * m_cfg.zOriginWeightFactor) +
                m_cfg.compatSeedWeight;

      // skip a bad quality seed if any of its constituents has a weight larger
      // than the seed weight
      if (weight < getBestSeedQuality(state.bestSeedQualityMap, bottomSp) &&
          weight < getBestSeedQuality(state.bestSeedQualityMap, middleSp) &&
          weight < getBestSeedQuality(state.bestSeedQualityMap, topSp)) {
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
  if (options.seedConfirmation && maxWeightSeed &&
      candidatesCollector.nHighQualityCandidates() == 0) {
    // if we have not yet reached our max number of seeds we add the new seed to
    // outCont

    candidatesCollector.push(bottomSp, middleSp, maxWeightTopSp, weightMax,
                             zOrigin, false);
  }
}

void TripletSeedFilter2::filter1SpFixed(
    const Options& options, State& state,
    std::vector<TripletCandidate2>& candidates, std::size_t numQualitySeeds,
    SeedContainer2& outputCollection) const {
  if (m_cfg.experimentCuts != nullptr) {
    m_cfg.experimentCuts->cutPerMiddleSp(candidates);
  }

  unsigned int maxSeeds = candidates.size();

  if (maxSeeds > m_cfg.maxSeedsPerSpM) {
    maxSeeds = m_cfg.maxSeedsPerSpM + 1;
  }

  // default filter removes the last seeds if maximum amount exceeded
  // ordering by weight by filterSeeds_2SpFixed means these are the lowest
  // weight seeds
  unsigned int numTotalSeeds = 0;
  for (const auto& [bottom, middle, top, bestSeedQuality, zOrigin,
                    qualitySeed] : candidates) {
    // stop if we reach the maximum number of seeds
    if (numTotalSeeds >= maxSeeds) {
      break;
    }

    if (options.seedConfirmation) {
      // continue if higher-quality seeds were found
      if (numQualitySeeds > 0 && !qualitySeed) {
        continue;
      }
      if (bestSeedQuality <
              getBestSeedQuality(state.bestSeedQualityMap, bottom) &&
          bestSeedQuality <
              getBestSeedQuality(state.bestSeedQualityMap, middle) &&
          bestSeedQuality < getBestSeedQuality(state.bestSeedQualityMap, top)) {
        continue;
      }
    }

    // set quality of seed components
    state.bestSeedQualityMap[bottom] = bestSeedQuality;
    state.bestSeedQualityMap[middle] = bestSeedQuality;
    state.bestSeedQualityMap[top] = bestSeedQuality;

    ACTS_VERBOSE("Adding seed: [b=" << bottom << ", m=" << middle << ", t="
                                    << top << "], quality=" << bestSeedQuality
                                    << ", vertexZ=" << zOrigin);

    auto seed = outputCollection.makeSeed(bottom, middle, top);
    seed.vertexZ() = zOrigin;
    seed.quality() = bestSeedQuality;

    ++numTotalSeeds;
  }

  ACTS_VERBOSE("Identified " << numTotalSeeds << " seeds");
}

}  // namespace Acts
