// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"

#include <memory>

namespace Acts {
// forward declaration to avoid cyclic dependence
template <typename T>
class SeedFilter;

template <typename SpacePoint>
struct SeedFinderOrthogonalConfig {
  std::shared_ptr<Acts::SeedFilter<SpacePoint>> seedFilter;

  // Seed Cuts
  // lower cutoff for seeds
  float minPt = 400. * Acts::UnitConstants::MeV;
  // cot of maximum theta angle
  // equivalent to 2.7 eta (pseudorapidity)
  float cotThetaMax = 7.40627;
  // minimum distance in r between middle and top SP in one seed
  float deltaRMinTopSP = 5 * Acts::UnitConstants::mm;
  // maximum distance in r between middle and top SP in one seed
  float deltaRMaxTopSP = 270 * Acts::UnitConstants::mm;
  // minimum distance in r between middle and bottom SP in one seed
  float deltaRMinBottomSP = 5 * Acts::UnitConstants::mm;
  // maximum distance in r between middle and bottom SP in one seed
  float deltaRMaxBottomSP = 270 * Acts::UnitConstants::mm;

  // impact parameter
  float impactMax = 20. * Acts::UnitConstants::mm;

  // how many sigmas of scattering angle should be considered?
  float sigmaScattering = 5;
  // Upper pt limit for scattering calculation
  float maxPtScattering = 10 * Acts::UnitConstants::GeV;

  // for how many seeds can one SpacePoint be the middle SpacePoint?
  unsigned int maxSeedsPerSpM = 5;

  // Geometry Settings
  // Detector ROI
  // limiting location of collision region in z
  float collisionRegionMin = -150 * Acts::UnitConstants::mm;
  float collisionRegionMax = +150 * Acts::UnitConstants::mm;
  float phiMin = -M_PI;
  float phiMax = M_PI;
  // limiting location of measurements
  float zMin = -2800 * Acts::UnitConstants::mm;
  float zMax = 2800 * Acts::UnitConstants::mm;
  float rMax = 600 * Acts::UnitConstants::mm;
  // WARNING: if rMin is smaller than impactMax, the bin size will be 2*pi,
  // which will make seeding very slow!
  float rMin = 33 * Acts::UnitConstants::mm;

  // z of last layers to avoid iterations
  std::pair<float, float> zOutermostLayers{-2700 * Acts::UnitConstants::mm,
                                           2700 * Acts::UnitConstants::mm};

  // radial range for middle SP
  // variable range based on SP radius
  bool useVariableMiddleSPRange = true;
  float deltaRMiddleMinSPRange = 10. * Acts::UnitConstants::mm;
  float deltaRMiddleMaxSPRange = 10. * Acts::UnitConstants::mm;
  // range defined in vector for each z region
  std::vector<std::vector<float>> rRangeMiddleSP;
  // range defined by rMinMiddle and rMaxMiddle
  float rMinMiddle = 60.f * Acts::UnitConstants::mm;
  float rMaxMiddle = 120.f * Acts::UnitConstants::mm;

  float deltaPhiMax = 0.085;

  // cut to the maximum value of delta z between SPs
  float deltaZMax =
      std::numeric_limits<float>::infinity() * Acts::UnitConstants::mm;

  // cut on bottom SPs in a certain (r, eta) region of the detector for fast
  // seeding
  bool fastTrackingCut = false;
  float fastTrackingRMin = 50. * Acts::UnitConstants::mm;
  float fastTrackingCotThetaMax = 2.13;

  // enable cut on the compatibility between interaction point and SPs
  bool interactionPointCut = false;

  // seed confirmation
  bool seedConfirmation = false;
  // parameters for central seed confirmation
  SeedConfirmationRangeConfig centralSeedConfirmationRange;
  // parameters for forward seed confirmation
  SeedConfirmationRangeConfig forwardSeedConfirmationRange;

  // average radiation lengths of material on the length of a seed. used for
  // scattering.
  // default is 5%
  // TODO: necessary to make amount of material dependent on detector region?
  float radLengthPerSeed = 0.05;

  // derived values, set on SeedFinder construction
  float highland = 0;
  float maxScatteringAngle2 = 0;

  bool isInInternalUnits = false;

  SeedFinderOrthogonalConfig calculateDerivedQuantities() const {
    if (not isInInternalUnits) {
      throw std::runtime_error(
          "SeedFinderOrthogonalConfig not in ACTS internal units in "
          "calculateDerivedQuantities");
    }
    SeedFinderOrthogonalConfig config = *this;
    // calculation of scattering using the highland formula
    // convert pT to p once theta angle is known
    config.highland = 13.6 * std::sqrt(radLengthPerSeed) *
                      (1 + 0.038 * std::log(radLengthPerSeed));
    config.maxScatteringAngle2 = std::pow(config.highland / config.minPt, 2);
    return config;
  }

  SeedFinderOrthogonalConfig toInternalUnits() const {
    if (isInInternalUnits) {
      throw std::runtime_error(
          "SeedFinderOrthogonalConfig already in ACTS internal units in "
          "toInternalUnits");
    }
    using namespace Acts::UnitLiterals;
    SeedFinderOrthogonalConfig config = *this;
    config.isInInternalUnits = true;
    config.minPt /= 1_MeV;
    config.deltaRMinTopSP /= 1_mm;
    config.deltaRMaxTopSP /= 1_mm;
    config.deltaRMinBottomSP /= 1_mm;
    config.deltaRMaxBottomSP /= 1_mm;
    config.impactMax /= 1_mm;
    config.maxPtScattering /= 1_MeV;
    config.collisionRegionMin /= 1_mm;
    config.collisionRegionMax /= 1_mm;
    config.zMin /= 1_mm;
    config.zMax /= 1_mm;
    config.rMax /= 1_mm;
    config.rMin /= 1_mm;

    return config;
  }
};

}  // namespace Acts
