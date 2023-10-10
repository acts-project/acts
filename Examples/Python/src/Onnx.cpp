// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/TrackFindingML/AmbiguityResolutionMLAlgorithm.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addOnnx(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");
  auto onnx = mex.def_submodule("_onnx");
  ctx.modules["onnx"] = onnx;

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::AmbiguityResolutionMLAlgorithm,
                                onnx, "AmbiguityResolutionMLAlgorithm",
                                inputTracks, inputDuplicateNN, outputTracks,
                                nMeasurementsMin);
}
}  // namespace Acts::Python
