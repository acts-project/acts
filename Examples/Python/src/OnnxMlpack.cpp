// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/TrackFindingML/AmbiguityResolutionMLDBScanAlgorithm.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addOnnxMlpack(Context& ctx) {
  auto [m, mex, onnx] = ctx.get("main", "examples", "onnx");
  auto mlpack = mex.def_submodule("_mlpack");

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::AmbiguityResolutionMLDBScanAlgorithm, mlpack,
      "AmbiguityResolutionMLDBScanAlgorithm", inputTracks, inputDuplicateNN,
      outputTracks, nMeasurementsMin, epsilonDBScan, minPointsDBScan);
}
}  // namespace Acts::Python
