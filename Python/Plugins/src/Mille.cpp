// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Mille/MillePedeResultReader.hpp"
#include "ActsPlugins/Mille/MillePedeSolver.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

using ActsPlugins::MillePedeResultReader;
using ActsPlugins::MillePedeSolver;
using ParameterResult = ActsPlugins::MillePedeResultReader::ParameterResult;

PYBIND11_MODULE(ActsPluginsPythonBindingsMille, mille) {
  {
    auto ps = py::class_<MillePedeSolver, std::shared_ptr<MillePedeSolver>>(
                  mille, "MillePedeSolver")
                  .def(py::init<std::unique_ptr<Acts::Logger>>())
                  .def("solve", &MillePedeSolver::solve);

    auto mr =
        py::class_<MillePedeSolver::Result>(ps, "Result").def(py::init<>());
    ACTS_PYTHON_STRUCT(mr, exitCode, exitStatus, exitMessage, resultsFile,
                       logFile, histoFile, evFile);

    auto sc =
        py::class_<MillePedeSolver::Config>(ps, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(sc, steeringFile, workDir, extraOpts, resFileName,
                       logFileName, histoFileName, evFileName);
  }

  {
    auto ms =
        py::class_<MillePedeResultReader,
                   std::shared_ptr<MillePedeResultReader>>(
            mille, "MillePedeResultReader")
            .def(py::init<std::unique_ptr<Acts::Logger>>())
            .def("readParameters", &MillePedeResultReader::readParameters);

    auto c =
        py::class_<ParameterResult>(ms, "ParameterResult").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, label, val, start, delta, sigma, nRecords);
  }
}
