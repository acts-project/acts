// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <limits>
#include <memory>
#include <numbers>
#include <vector>

namespace Acts {

// forward declaration to avoid cyclic dependence
template <typename T>
class SeedFilter;

/// @brief Structure that holds configuration parameters for the seed finder algorithm
template <typename SpacePoint>
struct SeedFinderConfig {
  std::shared_ptr<SeedFilter<SpacePoint>> seedFilter;

  /// Seeding parameters used in the space-point grid creation and bin finding

  /// Geometry Settings + Detector ROI
  /// (r, z, phi) range for limiting location of all measurements and grid
  /// creation
  float phiMin = -std::numbers::pi_v<float>;
  float phiMax = std::numbers::pi_v<float>;
  float zMin = -2800 * UnitConstants::mm;
  float zMax = 2800 * UnitConstants::mm;
  float rMax = 600 * UnitConstants::mm;
  /// WARNING: if rMin is smaller than impactMax, the bin size will be 2*pi,
  /// which will make seeding very slow!
  float rMin = 33 * UnitConstants::mm;

  /// Vector containing the z-bin edges for non equidistant binning in z
  std::vector<float> zBinEdges;

  /// Order of z bins to loop over when searching for SPs
  std::vector<std::size_t> zBinsCustomLooping = {};

  /// Radial bin size used in space-point grid
  float binSizeR = 1. * UnitConstants::mm;

  /// Seeding parameters used to define the region of interest for middle
  /// space-point

  /// Radial range for middle space-point
  /// The range can be defined manually with (rMinMiddle, rMaxMiddle). If
  /// useVariableMiddleSPRange is set to false and the vector rRangeMiddleSP is
  /// empty, we use (rMinMiddle, rMaxMiddle) to cut the middle space-points
  float rMinMiddle = 60.f * UnitConstants::mm;
  float rMaxMiddle = 120.f * UnitConstants::mm;
  /// If useVariableMiddleSPRange is set to false, the vector rRangeMiddleSP can
  /// be used to define a fixed r range for each z bin: {{rMin, rMax}, ...}
  bool useVariableMiddleSPRange = false;
  /// Range defined in vector for each z bin
  std::vector<std::vector<float>> rRangeMiddleSP;
  /// If useVariableMiddleSPRange is true, the radial range will be calculated
  /// based on the maximum and minimum r values of the space-points in the event
  /// and a deltaR (deltaRMiddleMinSPRange, deltaRMiddleMaxSPRange)
  float deltaRMiddleMinSPRange = 10. * UnitConstants::mm;
  float deltaRMiddleMaxSPRange = 10. * UnitConstants::mm;

  /// Seeding parameters used to define the cuts on space-point doublets

  /// Minimum radial distance between two doublet components (prefer
  /// deltaRMinTopSP and deltaRMinBottomSP to set separate values for outer and
  /// inner space-points)
  float deltaRMin = 5 * UnitConstants::mm;
  /// Maximum radial distance between two doublet components (prefer
  /// deltaRMaxTopSP and deltaRMacBottomSP to set separate values for outer and
  /// inner space-points)
  float deltaRMax = 270 * UnitConstants::mm;
  /// Minimum radial distance between middle-outer doublet components
  float deltaRMinTopSP = std::numeric_limits<float>::quiet_NaN();
  /// Maximum radial distance between middle-outer doublet components
  float deltaRMaxTopSP = std::numeric_limits<float>::quiet_NaN();
  /// Minimum radial distance between inner-middle doublet components
  float deltaRMinBottomSP = std::numeric_limits<float>::quiet_NaN();
  /// Maximum radial distance between inner-middle doublet components
  float deltaRMaxBottomSP = std::numeric_limits<float>::quiet_NaN();

  /// Maximum value of z-distance between space-points in doublet
  float deltaZMax = std::numeric_limits<float>::infinity() * UnitConstants::mm;

  /// Maximum allowed cotTheta between two space-points in doublet, used to
  /// check if forward angle is within bounds
  float cotThetaMax = 10.01788;  // equivalent to eta = 3 (pseudorapidity)

  /// Limiting location of collision region in z-axis used to check if doublet
  /// origin is within reasonable bounds
  float collisionRegionMin = -150 * UnitConstants::mm;
  float collisionRegionMax = +150 * UnitConstants::mm;

  /// Enable cut on the compatibility between interaction point and doublet,
  /// this is an useful approximation to speed up the seeding
  bool interactionPointCut = false;

  /// Seeding parameters used to define the cuts on space-point triplets

  /// Minimum transverse momentum (pT) used to check the r-z slope compatibility
  /// of triplets with maximum multiple scattering effect (produced by the
  /// minimum allowed pT particle) + a certain uncertainty term. Check the
  /// documentation for more information
  /// https://acts.readthedocs.io/en/latest/core/reconstruction/pattern_recognition/seeding.html
  float minPt = 400. * UnitConstants::MeV;
  /// Number of sigmas of scattering angle to be considered in the minimum pT
  /// scattering term
  float sigmaScattering = 5;
  /// Term that accounts for the thickness of scattering medium in radiation
  /// lengths in the Lynch & Dahl correction to the Highland equation default is
  /// 5%
  /// TODO: necessary to make amount of material dependent on detector region?
  float radLengthPerSeed = 0.05;
  /// Maximum transverse momentum for scattering calculation
  float maxPtScattering = 10 * UnitConstants::GeV;
  /// Maximum value of impact parameter estimation of the seed candidates
  float impactMax = 20. * UnitConstants::mm;
  /// Parameter which can loosen the tolerance of the track seed to form a
  /// helix. This is useful for e.g. misaligned seeding.
  float helixCutTolerance = 1.;

  /// Seeding parameters used for quality seed confirmation

  /// Enable quality seed confirmation, this is different than default seeding
  /// confirmation because it can also be defined for different (r, z) regions
  /// of the detector (e.g. forward or central region) by SeedConfirmationRange.
  /// Seeds are classified as "high-quality" seeds and normal quality seeds.
  /// Normal quality seeds are only selected if no other "high-quality" seeds
  /// has been found for that inner-middle doublet.
  bool seedConfirmation = false;
  /// Contains parameters for central seed confirmation
  SeedConfirmationRangeConfig centralSeedConfirmationRange;
  /// Contains parameters for forward seed confirmation
  SeedConfirmationRangeConfig forwardSeedConfirmationRange;
  /// Maximum number (minus one) of accepted seeds per middle space-point
  unsigned int maxSeedsPerSpM = 5;

  /// Other parameters

  /// Alignment uncertainties, used for uncertainties in the
  /// non-measurement-plane of the modules
  /// which otherwise would be 0
  /// will be added to spacepoint measurement uncertainties (and therefore also
  /// multiplied by sigmaError)
  /// FIXME: call align1 and align2
  float zAlign = 0 * UnitConstants::mm;
  float rAlign = 0 * UnitConstants::mm;
  /// used for measurement (+alignment) uncertainties.
  /// find seeds within 5sigma error ellipse
  float sigmaError = 5;

  /// derived values, set on SeedFinder construction
  float highland = 0;
  float maxScatteringAngle2 = 0;

  /// only for Cuda plugin
  int maxBlockSize = 1024;
  int nTrplPerSpBLimit = 100;
  int nAvgTrplPerSpBLimit = 2;

  /// Delegates for accessors to detailed information on double measurement that
  /// produced the space point.
  /// This is mainly referring to space points produced when combining
  /// measurement from strips on back-to-back modules.
  /// Enables setting of the following delegates.
  bool useDetailedDoubleMeasurementInfo = false;

  Delegate<bool(const SpacePoint&)> spacePointSelector{
      DelegateFuncTag<voidSpacePointSelector>{}};

  static bool voidSpacePointSelector(const SpacePoint& /*sp*/) { return true; }

  /// Tolerance parameter used to check the compatibility of space-point
  /// coordinates in xyz. This is only used in a detector specific check for
  /// strip modules
  float toleranceParam = 1.1 * UnitConstants::mm;

  // Delegate to apply experiment specific cuts during doublet finding
  Delegate<bool(float /*bottomRadius*/, float /*cotTheta*/)> experimentCuts{
      DelegateFuncTag<&noopExperimentCuts>{}};

  bool isInInternalUnits = true;
  //[[deprecated("SeedFinderConfig uses internal units")]]
  SeedFinderConfig toInternalUnits() const { return *this; }

  SeedFinderConfig calculateDerivedQuantities() const {
    SeedFinderConfig config = *this;
    config.highland = approximateHighlandScattering(config.radLengthPerSeed);
    const float maxScatteringAngle = config.highland / minPt;
    config.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;
    return config;
  }
};

struct SeedFinderOptions {
  // location of beam in x,y plane.
  // used as offset for Space Points
  Vector2 beamPos{0 * UnitConstants::mm, 0 * UnitConstants::mm};
  // field induction
  float bFieldInZ = 2 * UnitConstants::T;

  // derived quantities
  float pTPerHelixRadius = std::numeric_limits<float>::quiet_NaN();
  float minHelixDiameter2 = std::numeric_limits<float>::quiet_NaN();
  float pT2perRadius = std::numeric_limits<float>::quiet_NaN();
  float sigmapT2perRadius = std::numeric_limits<float>::quiet_NaN();
  float multipleScattering2 = std::numeric_limits<float>::quiet_NaN();

  bool isInInternalUnits = true;
  //[[deprecated("SeedFinderOptions uses internal units")]]
  SeedFinderOptions toInternalUnits() const { return *this; }

  template <typename Config>
  SeedFinderOptions calculateDerivedQuantities(const Config& config) const {
    using namespace UnitLiterals;

    SeedFinderOptions options = *this;
    // bFieldInZ is in (pT/radius) natively, no need for conversion
    options.pTPerHelixRadius = options.bFieldInZ;
    options.minHelixDiameter2 =
        std::pow(config.minPt * 2 / options.pTPerHelixRadius, 2) *
        config.helixCutTolerance;
    options.pT2perRadius =
        std::pow(config.highland / options.pTPerHelixRadius, 2);
    options.sigmapT2perRadius =
        options.pT2perRadius * std::pow(2 * config.sigmaScattering, 2);
    options.multipleScattering2 =
        config.maxScatteringAngle2 * std::pow(config.sigmaScattering, 2);
    return options;
  }
};

}  // namespace Acts
