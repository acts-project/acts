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

using ActsPlugins::ActsToMille::MillePedeResultReader;
using ActsPlugins::ActsToMille::MillePedeSolver;
using ActsPlugins::ActsToMille::MillePedeSteering;
using ActsPlugins::ActsToMille::mpParameterResult;

PYBIND11_MODULE(ActsPluginsPythonBindingsMille, mille) {
  {
    auto ps = py::class_<MillePedeSolver, std::shared_ptr<MillePedeSolver>>(
                  mille, "MillePedeSolver")
                  .def(py::init<Acts::Logging::Level>())
                  .def("solve", &MillePedeSolver::solve);

    auto mr =
        py::class_<MillePedeSolver::mpResult>(ps, "mpResult").def(py::init<>());
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
            .def(py::init<Acts::Logging::Level>())
            .def("readParameters", &MillePedeResultReader::readParameters);

    auto c = py::class_<mpParameterResult>(mille, "mpParameterResult")
                 .def(py::init<>());

    ACTS_PYTHON_STRUCT(c, label, val, start, delta, sigma, nRecords);
  }
  {
    auto ms = py::class_<MillePedeSteering, std::shared_ptr<MillePedeSteering>>(
                  mille, "MillePedeSteering")
                  .def(py::init<Acts::Logging::Level>())
                  .def("generateSteeringFile",
                       &MillePedeSteering::generateSteeringFile);

    auto e = py::class_<MillePedeSteering::equalityConstraint>(
                 ms, "equalityConstraint")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(e, labelsAndWeights, constraint);

    auto s =
        py::enum_<MillePedeSteering::Strategy>(ms, "Strategy")
            .value("inversion", MillePedeSteering::Strategy::inversion)
            .value("diagonalization",
                   MillePedeSteering::Strategy::diagonalization)
            .value("decomposition", MillePedeSteering::Strategy::decomposition)
            .value("fullMINRES", MillePedeSteering::Strategy::fullMINRES)
            .value("sparseMINRES", MillePedeSteering::Strategy::sparseMINRES)
            .value("fullMINRES_QLP",
                   MillePedeSteering::Strategy::fullMINRES_QLP)
            .value("sparseMINRES_QLP",
                   MillePedeSteering::Strategy::sparseMINRES_QLP)
            .value("fullLAPACK", MillePedeSteering::Strategy::fullLAPACK)
            .value("unpackedLAPACK",
                   MillePedeSteering::Strategy::unpackedLAPACK)
            .value("sparsePARDISO", MillePedeSteering::Strategy::sparsePARDISO);

    auto c =
        py::class_<MillePedeSteering::Config>(ms, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, strategy, minIterations, convergenceLimit, entriesCut,
                       outlierDownweighting, downweightFractionCut, nOMPthreads,
                       nIOthreads, matIter, printCounts, chi2Cut,
                       monitorResiduals, monitorPulls, skipEmptyCons,
                       countRecords, extraLines, constraints, inputFiles);
  }
}
