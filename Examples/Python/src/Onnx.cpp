// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingML/AmbiguityResolutionMLAlgorithm.hpp"
#include "ActsExamples/TrackFindingML/SeedFilterMLAlgorithm.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace ActsPython {

void addOnnx(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");
  auto onnx = mex.def_submodule("_onnx");
  ctx.modules["onnx"] = onnx;

  ACTS_PYTHON_DECLARE_ALGORITHM(
      AmbiguityResolutionMLAlgorithm, onnx, "AmbiguityResolutionMLAlgorithm",
      inputTracks, inputDuplicateNN, outputTracks, nMeasurementsMin);

  ACTS_PYTHON_DECLARE_ALGORITHM(SeedFilterMLAlgorithm, onnx,
                                "SeedFilterMLAlgorithm", inputTrackParameters,
                                inputSimSeeds, inputSeedFilterNN,
                                outputTrackParameters, outputSimSeeds,
                                epsilonDBScan, minPointsDBScan, minSeedScore);
}
}  // namespace ActsPython
