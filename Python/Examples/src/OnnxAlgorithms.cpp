// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingML/AmbiguityResolutionMLAlgorithm.hpp"
#include "ActsExamples/TrackFindingML/SeedFilterMLAlgorithm.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace ActsPython {

/// This adds the ONNX algorithms to the examples module
/// @param mex the context containing the module
void addOnnxAlgorithms(py::module_& mex) {
  auto onnx = mex.def_submodule("_onnx");

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::AmbiguityResolutionMLAlgorithm,
                                onnx, "AmbiguityResolutionMLAlgorithm",
                                inputTracks, inputDuplicateNN, outputTracks,
                                nMeasurementsMin);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::SeedFilterMLAlgorithm, onnx,
                                "SeedFilterMLAlgorithm", inputTrackParameters,
                                inputSimSeeds, inputSeedFilterNN,
                                outputTrackParameters, outputSimSeeds,
                                epsilonDBScan, minPointsDBScan, minSeedScore);

  onnx.def(
      "makeNeuralCalibrator",
      [](const char* modelPath, std::size_t nComp,
         std::vector<std::size_t> volumeIds)
          -> std::shared_ptr<ActsExamples::MeasurementCalibrator> {
        return std::make_shared<ActsExamples::NeuralCalibrator>(
            modelPath, nComp, volumeIds);
      },
      py::arg("modelPath"), py::arg("nComp") = 1,
      py::arg("volumeIds") = std::vector<std::size_t>({7, 8, 9}));
}
}  // namespace ActsPython
