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

  float rMinMiddle = 60.f * Acts::UnitConstants::mm;
  float rMaxMiddle = 120.f * Acts::UnitConstants::mm;

  float deltaPhiMax = 0.085;

  float bFieldInZ = 2.08 * Acts::UnitConstants::T;
  // location of beam in x,y plane.
  // used as offset for Space Points
  Acts::Vector2 beamPos{0 * Acts::UnitConstants::mm,
                        0 * Acts::UnitConstants::mm};

  // cut to the maximum value of delta z between SPs
  float deltaZMax =
      std::numeric_limits<float>::infinity() * Acts::UnitConstants::mm;

  // enable cut on the compatibility between interaction point and SPs
  bool interactionPointCut = false;

  // skip top SPs based on cotTheta sorting when producing triplets
  bool skipPreviousTopSP = false;

  // average radiation lengths of material on the length of a seed. used for
  // scattering.
  // default is 5%
  // TODO: necessary to make amount of material dependent on detector region?
  float radLengthPerSeed = 0.05;

  // derived values, set on SeedFinder construction
  float highland = 0;
  float maxScatteringAngle2 = 0;
  float pTPerHelixRadius = 0;
  float minHelixDiameter2 = 0;
  float pT2perRadius = 0;
  float sigmapT2perRadius = 0;

  SeedFinderOrthogonalConfig toInternalUnits() const {
    using namespace Acts::UnitLiterals;
    SeedFinderOrthogonalConfig config = *this;
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
    config.bFieldInZ /= 1000. * 1_T;

    config.beamPos[0] /= 1_mm;
    config.beamPos[1] /= 1_mm;

    return config;
  }
};

}  // namespace Acts
