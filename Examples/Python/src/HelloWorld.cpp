// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/HelloWorld/HelloLoggerAlgorithm.hpp"
#include "ActsExamples/HelloWorld/HelloRandomAlgorithm.hpp"
#include "ActsExamples/HelloWorld/HelloWhiteBoardAlgorithm.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

using namespace ActsExamples;

namespace Acts::Python {

/// This adds the digitization algorithms to the examples module
/// @param mex the examples module
void addHelloWorld(Context& ctx) {
  auto& mex = ctx.get("examples");
  // HelloLoggerAlgorithm: this class has no Config struct, declare directly
  py::class_<HelloLoggerAlgorithm, std::shared_ptr<HelloLoggerAlgorithm>,
             IAlgorithm>(mex, "HelloLoggerAlgorithm")
      .def(py::init<Acts::Logging::Level>(),
           py::arg("level") = Acts::Logging::INFO);

  // Use the marco to define the algorithm, input/output are properties of the
  // nested Config{} struct;
  ACTS_PYTHON_DECLARE_ALGORITHM(HelloWhiteBoardAlgorithm, mex,
                                "HelloWhiteBoardAlgorithm", input, output);

  // More complex configuration parameters for this one
  ACTS_PYTHON_DECLARE_ALGORITHM(
      HelloRandomAlgorithm, mex, "HelloRandomAlgorithm", randomNumbers,
      gaussParameters, uniformParameters, gammaParameters, poissonParameter,
      drawsPerEvent, output);
}

}  // namespace Acts::Python
