// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include <ActsExamples/EventData/NeuralCalibrator.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addOnnxNeuralCalibrator(Context &ctx) {
  auto [m, mex, onnx] = ctx.get("main", "examples", "onnx");

  onnx.def(
      "makeNeuralCalibrator",
      [](const char *modelPath, size_t nComp, std::vector<size_t> volumeIds)
          -> std::shared_ptr<MeasurementCalibrator> {
        return std::make_shared<NeuralCalibrator>(modelPath, nComp, volumeIds);
      },
      py::arg("modelPath"), py::arg("nComp") = 1,
      py::arg("volumeIds") = std::vector<size_t>({7, 8, 9}));
}
}  // namespace Acts::Python
