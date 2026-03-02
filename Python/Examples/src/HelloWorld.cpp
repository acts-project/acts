// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/HelloWorld/HelloLoggerAlgorithm.hpp"
#include "ActsExamples/HelloWorld/HelloRandomAlgorithm.hpp"
#include "ActsExamples/HelloWorld/HelloWhiteBoardAlgorithm.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace ActsPython {

void addHelloWorld(py::module& mex) {
  py::class_<HelloLoggerAlgorithm, IAlgorithm, std::shared_ptr<HelloLoggerAlgorithm>>(
      mex, "HelloLoggerAlgorithm")
      .def(py::init([](Logging::Level level) {
             return std::make_shared<HelloLoggerAlgorithm>(
                 getDefaultLogger("HelloLoggerAlgorithm", level));
           }),
           py::arg("level"))
      .def(py::init([](std::unique_ptr<const Logger> logger) {
             return std::make_shared<HelloLoggerAlgorithm>(std::move(logger));
           }),
           py::arg("logger"));

  ACTS_PYTHON_DECLARE_ALGORITHM(
      HelloRandomAlgorithm, mex, "HelloRandomAlgorithm", randomNumbers,
      gaussParameters, uniformParameters, gammaParameters, poissonParameter,
      drawsPerEvent, output);

  ACTS_PYTHON_DECLARE_ALGORITHM(HelloWhiteBoardAlgorithm, mex,
                                "HelloWhiteBoardAlgorithm", input, output);
}

}  // namespace ActsPython
