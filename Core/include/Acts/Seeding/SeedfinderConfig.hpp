// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <memory>
#include <vector>

namespace Acts {

// forward declaration to avoid cyclic dependence
class SeedFilter;

struct SeedfinderConfig {
  std::shared_ptr<Acts::SeedFilter> seedFilter;

  // Seed Cuts
  // lower cutoff for seeds
  float minPt = 400. * Acts::UnitConstants::MeV;
  // cot of maximum theta angle
  // equivalent to 2.7 eta (pseudorapidity)
  float cotThetaMax = 7.40627;
  // minimum distance in r between two measurements within one seed
  float deltaRMin = 5 * Acts::UnitConstants::mm;
  // maximum distance in r between two measurements within one seed
  float deltaRMax = 270 * Acts::UnitConstants::mm;
  // minimum distance in r between middle and top SP
  float deltaRMinTopSP = 5 * Acts::UnitConstants::mm;
  // maximum distance in r between middle and top SP
  float deltaRMaxTopSP = 270 * Acts::UnitConstants::mm;
  // minimum distance in r between middle and bottom SP
  float deltaRMinBottomSP = 5 * Acts::UnitConstants::mm;
  // maximum distance in r between middle and bottom SP
  float deltaRMaxBottomSP = 270 * Acts::UnitConstants::mm;
  // radial bin size for filling space point grid
  float binSizeR = 1. * Acts::UnitConstants::mm;
  // force sorting in R in space point grid bins
  bool forceRadialSorting = false;

  // radial range for middle SP
  std::vector<std::vector<float>> rRangeMiddleSP;
  bool useVariableMiddleSPRange = false;
  float deltaRMiddleMinSPRange = 10. * Acts::UnitConstants::mm;
  float deltaRMiddleMaxSPRange = 10. * Acts::UnitConstants::mm;

  // seed confirmation
  bool seedConfirmation = false;
  // parameters for central seed confirmation
  SeedConfirmationRangeConfig centralSeedConfirmationRange;
  // parameters for forward seed confirmation
  SeedConfirmationRangeConfig forwardSeedConfirmationRange;

  // cut to the maximum value of delta z between SPs
  float deltaZMax =
      std::numeric_limits<float>::infinity() * Acts::UnitConstants::mm;

  // enable cut on the compatibility between interaction point and SPs
  bool interactionPointCut = false;

  // use arithmetic average in the calculation of the squared error on the
  // difference in tan(theta)
  bool arithmeticAverageCotTheta = false;

  // non equidistant binning in z
  std::vector<float> zBinEdges;

  // additional cut to skip top SPs when producing triplets
  bool skipPreviousTopSP = false;

  // FIXME: this is not used yet
  //        float upperPtResolutionPerSeed = 20* Acts::GeV;

  // the delta for inverse helix radius up to which compared seeds
  // are considered to have a compatible radius. delta of inverse radius
  // leads to this value being the cutoff. unit is 1/mm. default value
  // of 0.00003 leads to all helices with radius>33m to be considered compatible

  // impact parameter
  float impactMax = 20. * Acts::UnitConstants::mm;

  // how many standard deviations of scattering angle should be considered?
  float sigmaScattering = 5;
  // how many standard deviations of measurement uncertainty + alignment
  // uncertainty should be considered?
  float sigmaError = 5;
  // Upper pt limit for scattering calculation
  float maxPtScattering = 10 * Acts::UnitConstants::GeV;

  // for how many seeds can one SpacePoint be the middle SpacePoint?
  unsigned int maxSeedsPerSpM = 5;

  // tolerance parameter used to check the compatibility of SPs coordinates in
  // xyz
  float toleranceParam = 1.1 * Acts::UnitConstants::mm;

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

  float bFieldInZ = 2.08 * Acts::UnitConstants::T;
  // location of beam in x,y plane.
  // used as offset for Space Points
  Acts::Vector2 beamPos{0 * Acts::UnitConstants::mm,
                        0 * Acts::UnitConstants::mm};

  boost::container::small_vector<size_t, 20> zBinsCustomLooping = {};

  // average radiation lengths of material on the length of a seed. used for
  // scattering.
  // default is 5%
  // TODO: necessary to make amount of material dependent on detector region?
  float radLengthPerSeed = 0.05;

  // derived values, set on Seedfinder construction
  float highland = 0;
  float maxScatteringAngle2 = 0;
  float pTPerHelixRadius = 0;
  float minHelixDiameter2 = 0;
  float pT2perRadius = 0;
  float sigmapT2perRadius = 0;

  // only for Cuda plugin
  int maxBlockSize = 1024;
  int nTrplPerSpBLimit = 100;
  int nAvgTrplPerSpBLimit = 2;

  // Delegates for accessors to detailed information on double measurement that
  // produced the space point.
  // This is mainly referring to space points produced when combining
  // measurement from strips on back-to-back modules.
  // Enables setting of the following delegates.
  bool useDetailedDoubleMeasurementInfo = false;
  // Returns half of the length of the top strip.
  Delegate<float(const SpacePoint&)> getTopHalfStripLength;
  // Returns half of the length of the bottom strip.
  Delegate<float(const SpacePoint&)> getBottomHalfStripLength;
  // Returns direction of the top strip.
  Delegate<Acts::Vector3(const SpacePoint&)> getTopStripDirection;
  // Returns direction of the bottom strip.
  Delegate<Acts::Vector3(const SpacePoint&)> getBottomStripDirection;
  // Returns distance between the centers of the two strips.
  Delegate<Acts::Vector3(const SpacePoint&)> getStripCenterDistance;
  // Returns position of the center of the top strip.
  Delegate<Acts::Vector3(const SpacePoint&)> getTopStripCenterPosition;

  SeedfinderConfig toInternalUnits() const {
    using namespace Acts::UnitLiterals;
    SeedfinderConfig config = *this;
    config.minPt /= 1_MeV;
    config.deltaRMin /= 1_mm;
    config.deltaRMax /= 1_mm;
    config.binSizeR /= 1_mm;
    config.deltaRMinTopSP /= 1_mm;
    config.deltaRMaxTopSP /= 1_mm;
    config.deltaRMinBottomSP /= 1_mm;
    config.deltaRMaxBottomSP /= 1_mm;
    config.deltaRMiddleMinSPRange /= 1_mm;
    config.deltaRMiddleMaxSPRange /= 1_mm;
    config.impactMax /= 1_mm;
    config.maxPtScattering /= 1_MeV;  // correct?
    config.collisionRegionMin /= 1_mm;
    config.collisionRegionMax /= 1_mm;
    config.zMin /= 1_mm;
    config.zMax /= 1_mm;
    config.rMax /= 1_mm;
    config.rMin /= 1_mm;
    config.bFieldInZ /= 1000. * 1_T;
    config.deltaZMax /= 1_mm;

    config.beamPos[0] /= 1_mm;
    config.beamPos[1] /= 1_mm;

    config.toleranceParam /= 1_mm;

    return config;
  }
};

}  // namespace Acts