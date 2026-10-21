// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/HashingPrototypeSeedingAlgorithm.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsHashing, hashing) {
  ACTS_PYTHON_DECLARE_ALGORITHM(
      HashingPrototypeSeedingAlgorithm, hashing,
      "HashingPrototypeSeedingAlgorithm", inputSpacePoints, outputSeeds,
      outputBuckets, bFieldInZ, minPt, cotThetaMax, impactMax, deltaRMin,
      deltaRMax, deltaRMinTop, deltaRMaxTop, deltaRMinBottom, deltaRMaxBottom,
      deltaZMin, deltaZMax, interactionPointCut, collisionRegionMin,
      collisionRegionMax, helixCutTolerance, sigmaScattering, radLengthPerSeed,
      toleranceParam, deltaInvHelixDiameter, compatSeedWeight,
      impactWeightFactor, zOriginWeightFactor, maxSeedsPerSpM, compatSeedLimit,
      seedWeightIncrement, numSeedIncrement, seedConfirmation,
      centralSeedConfirmationRange, forwardSeedConfirmationRange,
      maxSeedsPerSpMConf, maxQualitySeedsPerSpMConf,
      useDeltaRinsteadOfTopRadius, useExtraCuts, annoySeed, f, bucketSize,
      zBins, phiBins, layerRMin, layerRMax, layerZMin, layerZMax);
}
