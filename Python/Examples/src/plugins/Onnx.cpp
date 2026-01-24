// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/NeuralCalibrator.hpp"
#include "ActsExamples/TrackFindingML/AmbiguityResolutionMLAlgorithm.hpp"
#include "ActsExamples/TrackFindingML/SeedFilterMLAlgorithm.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsOnnx, onnx) {
  onnx.def(
      "makeNeuralCalibrator",
      [](const char *modelPath, std::size_t nComp,
         std::vector<std::size_t> volumeIds)
          -> std::shared_ptr<MeasurementCalibrator> {
        return std::make_shared<NeuralCalibrator>(modelPath, nComp, volumeIds);
      },
      py::arg("modelPath"), py::arg("nComp") = 1,
      py::arg("volumeIds") = std::vector<std::size_t>({7, 8, 9}));

  ACTS_PYTHON_DECLARE_ALGORITHM(
      AmbiguityResolutionMLAlgorithm, onnx, "AmbiguityResolutionMLAlgorithm",
      inputTracks, inputDuplicateNN, outputTracks, nMeasurementsMin);

  ACTS_PYTHON_DECLARE_ALGORITHM(SeedFilterMLAlgorithm, onnx,
                                "SeedFilterMLAlgorithm", inputTrackParameters,
                                inputSimSeeds, inputSeedFilterNN,
                                outputTrackParameters, outputSimSeeds,
                                epsilonDBScan, minPointsDBScan, minSeedScore);
}
