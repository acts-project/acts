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
#include <iostream>

namespace Acts {

// forward declaration to avoid cyclic dependence
template <typename T>
class SeedFilterNA60;

template <typename SpacePoint>
struct SeedFinderConfigNA60 {
  std::shared_ptr<Acts::SeedFilterNA60<SpacePoint>> seedFilter;

  // Seed Cuts
  // lower cutoff for seeds
  float minPt = 100. * Acts::UnitConstants::MeV;
  // cot of maximum theta angle
  // equivalent to 2.7 eta (pseudorapidity)
  float cotThetaMax = 7.40627;
  // minimum distance in r between two measurements within one seed
  float deltaYMin = 5 * Acts::UnitConstants::mm;
  // maximum distance in r between two measurements within one seed
  float deltaYMax = 10 * Acts::UnitConstants::mm; //just a large number
  // minimum distance in r between middle and top SP
  float deltaYMinTopSP = 70 * Acts::UnitConstants::mm; //std::numeric_limits<float>::quiet_NaN();
  // maximum distance in r between middle and top SP
  float deltaYMaxTopSP = 129 * Acts::UnitConstants::mm; //std::numeric_limits<float>::quiet_NaN();
  // minimum distance in r between middle and bottom SP
  float deltaYMinBottomSP = 60* Acts::UnitConstants::mm;//std::numeric_limits<float>::quiet_NaN();
  // maximum distance in r between middle and bottom SP
  float deltaYMaxBottomSP = 90* Acts::UnitConstants::mm;//std::numeric_limits<float>::quiet_NaN();
  // radial bin size for filling space point grid
  float binSizeR = 1. * Acts::UnitConstants::mm;

  // radial range for middle SP
  // variable range based on SP radius
  bool useVariableMiddleSPRange = false;
  float deltaYMiddleMinSPRange = 80. * Acts::UnitConstants::mm;
  float deltaYMiddleMaxSPRange = 210. * Acts::UnitConstants::mm;
  // range defined in vector for each z region
  std::vector<std::vector<float>> yRangeMiddleSP;
  // range defined by rMinMiddle and rMaxMiddle
  float yMinMiddle = 40.f * Acts::UnitConstants::mm; // ok for first 3 layers (or 1234)
  float yMaxMiddle = 60.f * Acts::UnitConstants::mm; // ok for first 3 layers

  // z of last layers to avoid iterations
  std::pair<float, float> zOutermostLayers{-150 * Acts::UnitConstants::mm,
                                           150 * Acts::UnitConstants::mm};

  // cut to the maximum value of delta z between SPs
  float deltaZMax =
      std::numeric_limits<float>::infinity() * Acts::UnitConstants::mm;

  // enable cut on the compatibility between interaction point and SPs
  bool interactionPointCut = false;

  // non equidistant binning in z
  std::vector<float> zBinEdges;

  bool verbose=false;
  // seed confirmation
  bool seedConfirmation = false;
  // parameters for central seed confirmation
  SeedConfirmationRangeConfig seedConfirmationRange;

  // FIXME: this is not used yet
  //        float upperPtResolutionPerSeed = 20* Acts::GeV;

  // the delta for inverse helix radius up to which compared seeds
  // are considered to have a compatible radius. delta of inverse radius
  // leads to this value being the cutoff. unit is 1/mm. default value
  // of 0.00003 leads to all helices with radius>33m to be considered compatible

  // impact parameter
  float impactMax = 20. * Acts::UnitConstants::mm;

  // how many sigmas of scattering angle should be considered?
  float sigmaScattering = 5;
  // Upper pt limit for scattering calculation
  float maxPtScattering = 10 * Acts::UnitConstants::GeV;

  // for how many seeds can one SpacePoint be the middle SpacePoint?
  unsigned int maxSeedsPerSpM = 5;

  // tolerance parameter used to check the compatibility of SPs coordinates in
  // xyz
  float toleranceParam = 1.1 * Acts::UnitConstants::mm;

  // Parameter which can loosen the tolerance of the track seed to form to a
  // helix, useful for (e.g.) misaligned seeding
  float helixCut = 1.;

  // Geometry Settings
  // Detector ROI
  // limiting location of collision region in z
  float collisionRegionMin = -2 * Acts::UnitConstants::mm;
  float collisionRegionMax = +2 * Acts::UnitConstants::mm;
  float phiMin = -M_PI;
  float phiMax = M_PI;
  // limiting location of measurements
  float zMin = -2800 * Acts::UnitConstants::mm;
  float zMax = 2800 * Acts::UnitConstants::mm;
  float xMin = -2800 * Acts::UnitConstants::mm;
  float xMax = 2800 * Acts::UnitConstants::mm;
  float rMax = 600 * Acts::UnitConstants::mm;
  // WARNING: if rMin is smaller than impactMax, the bin size will be 2*pi,
  // which will make seeding very slow!
  float rMin = 33 * Acts::UnitConstants::mm;

  // Order of z bins to loop over when searching for SPs
  std::vector<size_t> zBinsCustomLooping = {};
  // Number of Z bins to skip the search for middle SPs
  std::size_t skipZMiddleBinSearch = 0;

  // average radiation lengths of material on the length of a seed. used for
  // scattering.
  // default is 5%
  // TODO: necessary to make amount of material dependent on detector region?
  float radLengthPerSeed = 0.01;
//  float radLengthPerSeed = 0.05;
  // alignment uncertainties, used for uncertainties in the
  // non-measurement-plane of the modules
  // which otherwise would be 0
  // will be added to spacepoint measurement uncertainties (and therefore also
  // multiplied by sigmaError)
  // FIXME: call align1 and align2
  float zAlign = 0 * Acts::UnitConstants::mm;
  float rAlign = 0 * Acts::UnitConstants::mm;
  // used for measurement (+alignment) uncertainties.
  // find seeds within 5sigma error ellipse
  float sigmaError = 5;

  // derived values, set on SeedFinder construction
  float highland = 0;
  float maxScatteringAngle2 = 0;

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

  bool isInInternalUnits = false;

  SeedFinderConfigNA60 toInternalUnits() const {
    if (isInInternalUnits) {
      throw std::runtime_error(
          "Repeated conversion to internal units for SeedFinderConfigNA60");
    }
    // Make sure the shared ptr to the seed filter is not a nullptr
    // And make sure the seed filter config is in internal units as well
    if (not seedFilter) {
      throw std::runtime_error(
          "Invalid values for the seed filter inside the seed filter config: "
          "nullptr");
    }
    if (not seedFilter->getSeedFilterConfig().isInInternalUnits) {
      throw std::runtime_error(
          "The internal Seed Filter configuration, contained in the seed "
          "finder config, is not in internal units.");
    }

    using namespace Acts::UnitLiterals;
    SeedFinderConfigNA60 config = *this;
    config.isInInternalUnits = true;
    config.minPt /= 1_MeV;
    config.deltaYMin /= 1_mm;
    config.deltaYMax /= 1_mm;
    config.binSizeR /= 1_mm;
    config.deltaYMinTopSP /= 1_mm;
    config.deltaYMaxTopSP /= 1_mm;
    config.deltaYMinBottomSP /= 1_mm;
    config.deltaYMaxBottomSP /= 1_mm;
    config.deltaYMiddleMinSPRange /= 1_mm;
    config.deltaYMiddleMaxSPRange /= 1_mm;
    config.impactMax /= 1_mm;
    config.maxPtScattering /= 1_MeV;  // correct?
    config.collisionRegionMin /= 1_mm;
    config.collisionRegionMax /= 1_mm;
    config.zMin /= 1_mm;
    config.zMax /= 1_mm;
    config.rMax /= 1_mm;
    config.rMin /= 1_mm;
    config.deltaZMax /= 1_mm;

    config.zAlign /= 1_mm;
    config.rAlign /= 1_mm;

    config.toleranceParam /= 1_mm;

    return config;
  }
  SeedFinderConfigNA60 calculateDerivedQuantities() const {
    SeedFinderConfigNA60 config = *this;
    // calculation of scattering using the highland formula
    // convert pT to p once theta angle is known
    config.highland = 13.6 * std::sqrt(radLengthPerSeed) *
                      (1 + 0.038 * std::log(radLengthPerSeed));
    const float maxScatteringAngle = config.highland / minPt;
    config.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;
    return config;
  }
};

struct SeedFinderOptionsNA60 {
  // location of beam in x,y plane.
  // used as offset for Space Points
  Acts::Vector2 beamPos{0 * Acts::UnitConstants::mm,
                        0 * Acts::UnitConstants::mm};
  // field induction
  float bFieldInZ = 2.08 * Acts::UnitConstants::T;

  // derived quantities
  float pTPerHelixRadius = std::numeric_limits<float>::quiet_NaN();
  float minHelixDiameter2 = std::numeric_limits<float>::quiet_NaN();
  float pT2perRadius = std::numeric_limits<float>::quiet_NaN();
  float sigmapT2perRadius = std::numeric_limits<float>::quiet_NaN();
  float multipleScattering2 = std::numeric_limits<float>::quiet_NaN();

  bool isInInternalUnits = false;

  SeedFinderOptionsNA60 toInternalUnits() const {
    if (isInInternalUnits) {
      throw std::runtime_error(
          "Repeated conversion to internal units for SeedFinderOptionsNA60");
    }
    using namespace Acts::UnitLiterals;
    SeedFinderOptionsNA60 options = *this;
    options.isInInternalUnits = true;
    options.beamPos[0] /= 1_mm;
    options.beamPos[1] /= 1_mm;

    options.bFieldInZ /= 1000. * 1_T;
    
      std::cout << "NA60+_SeedFinderOptionsNA60 options.bFieldInZ= " << options.bFieldInZ << std::endl;
    return options;
  }

  template <typename Config>
  SeedFinderOptionsNA60 calculateDerivedQuantities(const Config& config) const {
    using namespace Acts::UnitLiterals;

    if (!isInInternalUnits) {
      throw std::runtime_error(
          "Derived quantities in SeedFinderOptionsNA60 can only be calculated from "
          "Acts internal units");
    }
    SeedFinderOptionsNA60 options = *this;
    // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
    // millimeter
    // TODO: change using ACTS units
	  std::cout << "NA60+_options.bFieldInZ= " << options.bFieldInZ << std::endl;
    options.pTPerHelixRadius = 1_T * 1e6 * options.bFieldInZ;
	  std::cout << "NA60+_1_T= " << 1_T << std::endl;
	  std::cout << "NA60+_options.pTPerHelixRadius= " << options.pTPerHelixRadius << std::endl;
	  std::cout << "NA60+_config.minPt= " << config.minPt << std::endl;
    options.minHelixDiameter2 =
        std::pow(config.minPt * 2 / options.pTPerHelixRadius, 2);
	  std::cout << "NA60+_options.minHelixDiameter2= " << options.minHelixDiameter2 << std::endl;
    options.pT2perRadius =
        std::pow(config.highland / options.pTPerHelixRadius, 2);
    options.sigmapT2perRadius =
        options.pT2perRadius * std::pow(2 * config.sigmaScattering, 2);
    options.multipleScattering2 =
        config.maxScatteringAngle2 * std::pow(config.sigmaScattering, 2);
	  std::cout << "NA60+_options.multipleScattering2= " << options.multipleScattering2 << std::endl;
	
    return options;
  }
};

}  // namespace Acts
