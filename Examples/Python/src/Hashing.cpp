// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingAlgorithmHashing.hpp"
#include "ActsPlugins/Hashing/HashingAlgorithmConfig.hpp"
#include "ActsPlugins/Hashing/HashingTrainingConfig.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsPlugins;
using namespace ActsExamples;

namespace ActsPython {

void addHashing(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  auto hashingModule = m.def_submodule("hashing");
  auto hashingExampleModule = mex.def_submodule("_hashing");

  {
    using Config = HashingAlgorithmConfig;
    auto c = py::class_<Config>(hashingModule, "HashingAlgorithmConfig")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, bucketSize, zBins, phiBins);
    patchKwargsConstructor(c);
  }

  {
    using Config = HashingTrainingConfig;
    auto c = py::class_<Config>(hashingModule, "HashingTrainingConfig")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, annoySeed, f);
    patchKwargsConstructor(c);
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(
      SeedingAlgorithmHashing, hashingExampleModule, "SeedingAlgorithmHashing",
      inputSpacePoints, outputSeeds, outputBuckets, seedFilterConfig,
      seedFinderConfig, seedFinderOptions, gridConfig, gridOptions,
      allowSeparateRMax, zBinNeighborsTop, zBinNeighborsBottom, numPhiNeighbors,
      hashingConfig, hashingTrainingConfig);
}

}  // namespace ActsPython
