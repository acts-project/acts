// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// The triplet seeding cuts the seeding benchmarks share. The grid and the
/// orthogonal seeder build the same doublet finders, triplet finder and filter
/// and differ only in which space points they offer each other, so sharing the
/// cuts makes a comparison between them one of the binning alone.

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding/DoubletSeedFinder.hpp"
#include "Acts/Seeding/TripletSeedFinder.hpp"

#include <cstddef>
#include <limits>
#include <numbers>
#include <utility>
#include <vector>

namespace ActsTests::TripletSeedingConfig {

using namespace Acts::UnitLiterals;

/// The ATLAS ITk pixel seeding configuration: the values of
/// `ActsTrk::GridTripletSeedingTool` in Athena after `ActsPixelSeedingToolCfg`
/// has applied its pixel overrides. Only the properties that reach ACTS are
/// kept.
struct ItkPixelConfig {
  // -- shared between the grid, the doublet finders and the filter
  float minPt = 900_MeV;
  float cotThetaMax = 27.2899f;  // eta = 4
  float impactMax = 5_mm;        // ITk main pass maxPrimaryImpactSeed
  float bFieldInZ = 2_T;

  // -- grid
  float gridRMax = 320_mm;
  float zMin = -3000_mm;
  float zMax = 3000_mm;
  float deltaRMax = 280_mm;
  float phiMin = -std::numbers::pi_v<float>;
  float phiMax = std::numbers::pi_v<float>;
  int phiBinDeflectionCoverage = 3;
  int maxPhiBins = 200;
  int numPhiNeighbors = 1;
  std::vector<float> zBinEdges{-3000., -2700., -2500., -1400., -925.,
                               -500.,  -250.,  250.,   500.,   925.,
                               1400.,  2500.,  2700.,  3000.};
  std::vector<float> rBinEdges{0., 320.};
  std::vector<std::size_t> zBinsCustomLooping{2,  3, 4, 5, 12, 11,
                                              10, 9, 7, 6, 8};
  std::vector<std::size_t> rBinsCustomLooping{1};
  std::vector<std::pair<int, int>> zBinNeighborsTop{
      {0, 0}, {-1, 0}, {-2, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {-1, 1},
      {0, 1}, {0, 1},  {0, 1},  {0, 2},  {0, 1},  {0, 0}};
  std::vector<std::pair<int, int>> zBinNeighborsBottom{
      {0, 0},  {0, 1},  {0, 1},  {0, 1},  {0, 1},  {0, 1}, {0, 0},
      {-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {0, 0}};
  std::vector<std::pair<int, int>> rBinNeighborsTop{{0, 0}};
  std::vector<std::pair<int, int>> rBinNeighborsBottom{{0, 0}};

  // -- doublet finders
  float deltaRMinTopSP = 6_mm;
  float deltaRMaxTopSP = 280_mm;
  float deltaRMinBottomSP = 6_mm;
  float deltaRMaxBottomSP = 150_mm;
  // the pixel configuration switches the longitudinal doublet cut off
  float deltaZMax = std::numeric_limits<float>::infinity();
  float collisionRegionMin = -200_mm;
  float collisionRegionMax = 200_mm;
  bool interactionPointCut = true;
  float helixCutTolerance = 1.f;

  // -- triplet finder
  float sigmaScattering = 2.f;
  float radLengthPerSeed = 0.098045f;
  float toleranceParam = 1.1_mm;

  // -- middle space point range
  bool useVariableMiddleSPRange = false;
  float deltaRMiddleMinSPRange = 10.f;
  float deltaRMiddleMaxSPRange = 10.f;
  std::vector<std::pair<float, float>> rRangeMiddleSP{
      {0, 0},    {140, 260}, {40, 260}, {40, 260}, {40, 260},
      {40, 260}, {70, 260},  {40, 260}, {40, 260}, {40, 260},
      {40, 260}, {140, 260}, {0, 0}};

  // -- filter
  float deltaRMin = 20_mm;
  float deltaInvHelixDiameter = 0.00003f;
  float compatSeedWeight = 100.f;
  float impactWeightFactor = 100.f;
  float zOriginWeightFactor = 1.f;
  std::size_t maxSeedsPerSpM = 4;
  std::size_t compatSeedLimit = 3;
  float seedWeightIncrement = 0.f;
  // the pixel configuration disables the increment by putting it out of reach
  float numSeedIncrement = std::numeric_limits<float>::infinity();
  bool seedConfirmation = true;
  std::size_t maxSeedsPerSpMConf = 5;
  std::size_t maxQualitySeedsPerSpMConf = 5;
  bool useDeltaRinsteadOfTopRadius = true;
  float absDeltaEtaWeightFactor = 0.f;
  float absDeltaEtaMinImpact = 2_mm;
  /// Drop seeds that are the best one for none of their three space points
  bool seedQualitySelection = true;
};

inline Acts::DoubletSeedFinder::Config makeBottomDoubletConfig(
    const ItkPixelConfig& cfg) {
  Acts::DoubletSeedFinder::Config doubletCfg;
  doubletCfg.spacePointsSortedByRadius = true;
  doubletCfg.candidateDirection = Acts::Direction::Backward();
  doubletCfg.deltaRMin = cfg.deltaRMinBottomSP;
  doubletCfg.deltaRMax = cfg.deltaRMaxBottomSP;
  doubletCfg.deltaZMin = -cfg.deltaZMax;
  doubletCfg.deltaZMax = cfg.deltaZMax;
  doubletCfg.impactMax = cfg.impactMax;
  doubletCfg.interactionPointCut = cfg.interactionPointCut;
  doubletCfg.collisionRegionMin = cfg.collisionRegionMin;
  doubletCfg.collisionRegionMax = cfg.collisionRegionMax;
  doubletCfg.cotThetaMax = cfg.cotThetaMax;
  doubletCfg.minPt = cfg.minPt;
  doubletCfg.helixCutTolerance = cfg.helixCutTolerance;
  return doubletCfg;
}

inline Acts::TripletSeedFinder::Config makeTripletConfig(
    const ItkPixelConfig& cfg) {
  Acts::TripletSeedFinder::Config tripletCfg;
  tripletCfg.useStripInfo = false;
  tripletCfg.sortedByCotTheta = true;
  tripletCfg.minPt = cfg.minPt;
  tripletCfg.sigmaScattering = cfg.sigmaScattering;
  tripletCfg.radLengthPerSeed = cfg.radLengthPerSeed;
  tripletCfg.impactMax = cfg.impactMax;
  tripletCfg.helixCutTolerance = cfg.helixCutTolerance;
  tripletCfg.toleranceParam = cfg.toleranceParam;
  return tripletCfg;
}

inline Acts::BroadTripletSeedFilter::Config makeFilterConfig(
    const ItkPixelConfig& cfg) {
  Acts::BroadTripletSeedFilter::Config filterCfg;
  filterCfg.deltaInvHelixDiameter = cfg.deltaInvHelixDiameter;
  filterCfg.deltaRMin = cfg.deltaRMin;
  filterCfg.compatSeedWeight = cfg.compatSeedWeight;
  filterCfg.impactWeightFactor = cfg.impactWeightFactor;
  filterCfg.zOriginWeightFactor = cfg.zOriginWeightFactor;
  filterCfg.maxSeedsPerSpM = cfg.maxSeedsPerSpM;
  filterCfg.compatSeedLimit = cfg.compatSeedLimit;
  filterCfg.seedWeightIncrement = cfg.seedWeightIncrement;
  filterCfg.numSeedIncrement = cfg.numSeedIncrement;
  filterCfg.seedConfirmation = cfg.seedConfirmation;
  filterCfg.maxSeedsPerSpMConf = cfg.maxSeedsPerSpMConf;
  filterCfg.maxQualitySeedsPerSpMConf = cfg.maxQualitySeedsPerSpMConf;
  filterCfg.useDeltaRinsteadOfTopRadius = cfg.useDeltaRinsteadOfTopRadius;
  filterCfg.absDeltaEtaWeightFactor = cfg.absDeltaEtaWeightFactor;
  filterCfg.absDeltaEtaMinImpact = cfg.absDeltaEtaMinImpact;

  // central seed confirmation
  filterCfg.centralSeedConfirmationRange.zMinSeedConf = -250_mm;
  filterCfg.centralSeedConfirmationRange.zMaxSeedConf = 250_mm;
  filterCfg.centralSeedConfirmationRange.rMaxSeedConf = 140_mm;
  filterCfg.centralSeedConfirmationRange.nTopForLargeR = 1;
  filterCfg.centralSeedConfirmationRange.nTopForSmallR = 2;
  filterCfg.centralSeedConfirmationRange.seedConfMinBottomRadius = 60_mm;
  filterCfg.centralSeedConfirmationRange.seedConfMaxZOrigin = 150_mm;
  filterCfg.centralSeedConfirmationRange.minImpactSeedConf = 1_mm;
  // forward seed confirmation, identical but without the z restriction
  filterCfg.forwardSeedConfirmationRange =
      filterCfg.centralSeedConfirmationRange;
  filterCfg.forwardSeedConfirmationRange.zMinSeedConf = -3000_mm;
  filterCfg.forwardSeedConfirmationRange.zMaxSeedConf = 3000_mm;

  return filterCfg;
}

}  // namespace ActsTests::TripletSeedingConfig
