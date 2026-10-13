// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Mille/MillePedeResultReader.hpp"
#include "ActsPlugins/Mille/MillePedeSolver.hpp"
#include "ActsPlugins/Mille/MillePedeSteering.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

using ActsPlugins::MillePedeSolver;

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
                       redirectStdout, logFileName, histoFileName, evFileName);
  }

  {
    // unwrap Acts::Result - will throw in case of an error-result
    mille.def("readMillePedeResult", [](const std::filesystem::path& mpFile,
                                        const Acts::Logger& logger) {
      return ActsPlugins::readMillePedeResult(mpFile, logger).value();
    });

    auto c = py::class_<ActsPlugins::MillePedeParameterResult>(
                 mille, "MillePedeParameterResult")
                 .def(py::init<>());

    ACTS_PYTHON_STRUCT(c, label, val, start, delta, sigma, nRecords);
  }
  {
    using enum ActsPlugins::MillePedeSolutionStrategy;
    mille.def("generateMillePedeSteeringFile",
              &ActsPlugins::generateMillePedeSteeringFile);

    auto e = py::class_<ActsPlugins::MillePedeEqualityConstraint>(
                 mille, "MillePedeEqualityConstraint")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(e, labelsAndWeights, constraint);

    auto s = py::enum_<ActsPlugins::MillePedeSolutionStrategy>(
                 mille, "MillePedeSolutionStrategy")
                 .value("Inversion", Inversion)
                 .value("Diagonalization", Diagonalization)
                 .value("Decomposition", Decomposition)
                 .value("FullMinRes", FullMinRes)
                 .value("SparseMinRes", SparseMinRes)
                 .value("FullMinResQlp", FullMinResQlp)
                 .value("SparseMinResQlp", SparseMinResQlp)
                 .value("FullLapack", FullLapack)
                 .value("UnpackedLapack", UnpackedLapack)
                 .value("SparsePardiso", SparsePardiso);

    auto c = py::class_<ActsPlugins::MillePedeSteeringConfig>(
                 mille, "MillePedeSteeringConfig")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, strategy, minIterations, convergenceLimit, entriesCut,
                       outlierDownweighting, downweightFractionCut, nOMPthreads,
                       nIOthreads, matIter, printCounts, chi2Cut,
                       monitorResiduals, monitorPulls, skipEmptyCons,
                       countRecords, extraLines, constraints, inputFiles);
  }
}
