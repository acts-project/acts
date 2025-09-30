// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/HashingPrototypeSeedingAlgorithm.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace ActsPython {

void addHashing(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  auto hashingExampleModule = mex.def_submodule("_hashing");

  ACTS_PYTHON_DECLARE_ALGORITHM(
      HashingPrototypeSeedingAlgorithm, hashingExampleModule,
      "HashingPrototypeSeedingAlgorithm", inputSpacePoints, outputSeeds,
      bFieldInZ, minPt, cotThetaMax, impactMax, deltaRMin, deltaRMax,
      deltaRMinTop, deltaRMaxTop, deltaRMinBottom, deltaRMaxBottom, rMin, rMax,
      zMin, zMax, phiMin, phiMax, phiBinDeflectionCoverage, maxPhiBins,
      zBinNeighborsTop, zBinNeighborsBottom, numPhiNeighbors, zBinEdges,
      zBinsCustomLooping, rMinMiddle, rMaxMiddle, useVariableMiddleSPRange,
      rRangeMiddleSP, deltaRMiddleMinSPRange, deltaRMiddleMaxSPRange, deltaZMin,
      deltaZMax, interactionPointCut, collisionRegionMin, collisionRegionMax,
      helixCutTolerance, sigmaScattering, radLengthPerSeed, toleranceParam,
      deltaInvHelixDiameter, compatSeedWeight, impactWeightFactor,
      zOriginWeightFactor, maxSeedsPerSpM, compatSeedLimit, seedWeightIncrement,
      numSeedIncrement, seedConfirmation, centralSeedConfirmationRange,
      forwardSeedConfirmationRange, maxSeedsPerSpMConf,
      maxQualitySeedsPerSpMConf, useDeltaRinsteadOfTopRadius, useExtraCuts,
      annoySeed, f, bucketSize, zBins, phiBins, layerRMin, layerRMax, layerZMin,
      layerZMax);
}

}  // namespace ActsPython
