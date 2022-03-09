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

#include <memory>

namespace Acts {

struct SeedConfirmationRange {
  float zMinSeedConf;
  float zMaxSeedConf;
  float rMaxSeedConf;
  size_t nTopSeedConf;
  size_t nTopForLargeR;
  size_t nTopForSmallR;
};

// forward declaration to avoid cyclic dependence
template <typename T>
class SeedFilter;

template <typename SpacePoint>
struct RegionalParameters {
  std::shared_ptr<Acts::SeedFilter<SpacePoint>> seedFilter;
  // Seed Cuts
  // lower cutoff for seeds
  float minPt = 400. * Acts::UnitConstants::MeV;
  // minimum distance in r between middle and top SP
  float deltaRMinTopSP = 5 * Acts::UnitConstants::mm;
  // maximum distance in r between middle and top SP
  float deltaRMaxTopSP = 270 * Acts::UnitConstants::mm;
  // minimum distance in r between middle and bottom SP
  float deltaRMinBottomSP = 5 * Acts::UnitConstants::mm;
  // maximum distance in r between middle and bottom SP
  float deltaRMaxBottomSP = 270 * Acts::UnitConstants::mm;
  // average radiation lengths of material on the length of a seed. used for
  // scattering, default is 5%
  float radLengthPerSeed = 0.05;
  // how many sigmas of scattering angle should be considered?
  float sigmaScattering = 5;
  // Upper pt limit for scattering calculation
  float maxPtScattering = 10 * Acts::UnitConstants::GeV;
  // Average B field in Z
  float bFieldInZ = 2.08 * Acts::UnitConstants::T;

  // derived values, set on Seedfinder construction
  float highland = 0;
  float maxScatteringAngle2 = 0;
  float pTPerHelixRadius = 0;
  float minHelixDiameter2 = 0;
  float pT2perRadius = 0;
  float sigmapT2perRadius = 0;
};

template <typename SpacePoint>
struct SeedfinderConfig {
  // Vector countaning the boundary in Z of the different parameters region
  std::vector<float> zboundaries;
  // Vector countaning the boundary in R of the different parameters region
  std::vector<float> rboundaries;

  /// Return for a given middle space point the associated parameters region
  ///
  /// @param rM R position of the middle space point
  /// @param zM Z position of the middle space point
  ///
  /// @return the output as an output track
  size_t regionalBin(float rM, float zM) const {
    size_t bin = 0;
    size_t binZ = zboundaries.size() + 1;
    for (auto rboundary : rboundaries) {
      if (rM > rboundary) {
        bin = bin + binZ;
      }
    }
    for (auto zboundary : zboundaries) {
      if (zM > zboundary) {
        bin++;
      }
    }
    return bin;
  };

  // Vector of regional parameters for the configuration.
  // The bin of the region are calculated as BinR x TotalBinZ + BinZ
  std::vector<RegionalParameters<SpacePoint>> regionalParameters;

  // cot of maximum theta angle
  // equivalent to 2.7 eta (pseudorapidity)
  float cotThetaMax = 7.40627;

  // radial range for middle SP
  std::vector<std::vector<float>> rRangeMiddleSP;
  bool useVariableMiddleSPRange = false;
  float deltaRMiddleSPRange = 10. * Acts::UnitConstants::mm;

  // seed confirmation
  bool seedConfirmation = false;
  // parameters for central seed confirmation
  SeedConfirmationRange centralSeedConfirmationRange;
  // parameters for forward seed confirmation
  SeedConfirmationRange forwardSeedConfirmationRange;

  // non equidistant binning in z
  std::vector<float> zBinEdges;

  // sort the SP in transformCoordinates method and enables compatibility cuts
  // based on the sorting of cotTheta
  bool enableCutsForSortedSP = false;

  // impact parameter
  float impactMax = 20. * Acts::UnitConstants::mm;

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

  // location of beam in x,y plane.
  // used as offset for Space Points
  Acts::Vector2 beamPos{0 * Acts::UnitConstants::mm,
                        0 * Acts::UnitConstants::mm};

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

  // only for Cuda plugin
  int maxBlockSize = 1024;
  int nTrplPerSpBLimit = 100;
  int nAvgTrplPerSpBLimit = 2;

  SeedfinderConfig toInternalUnits() const {
    using namespace Acts::UnitLiterals;
    SeedfinderConfig config = *this;

    for (auto rparams : config.regionalParameters) {
      rparams.minPt /= 1_MeV;
      rparams.deltaRMinTopSP /= 1_mm;
      rparams.deltaRMaxTopSP /= 1_mm;
      rparams.deltaRMinBottomSP /= 1_mm;
      rparams.deltaRMaxBottomSP /= 1_mm;
      rparams.bFieldInZ /= 1000. * 1_T;
      rparams.maxPtScattering /= 1_MeV;  // correct?
    }

    config.deltaRMiddleSPRange /= 1_mm;
    config.impactMax /= 1_mm;
    config.collisionRegionMin /= 1_mm;
    config.collisionRegionMax /= 1_mm;
    config.zMin /= 1_mm;
    config.zMax /= 1_mm;
    config.rMax /= 1_mm;
    config.rMin /= 1_mm;

    config.beamPos[0] /= 1_mm;
    config.beamPos[1] /= 1_mm;

    config.zAlign /= 1_mm;
    config.rAlign /= 1_mm;

    return config;
  }
};

}  // namespace Acts
