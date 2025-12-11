// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingAlgorithmHashing.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsPlugins;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsHashing, hashing) {
  ACTS_PYTHON_DECLARE_ALGORITHM(
      SeedingAlgorithmHashing, hashing, "SeedingAlgorithmHashing",
      inputSpacePoints, outputSeeds, outputBuckets, seedFilterConfig,
      seedFinderConfig, seedFinderOptions, gridConfig, gridOptions,
      allowSeparateRMax, zBinNeighborsTop, zBinNeighborsBottom, numPhiNeighbors,
      hashingConfig, hashingTrainingConfig);
}
