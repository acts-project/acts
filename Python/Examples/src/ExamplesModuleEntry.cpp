// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPython/Utilities/Helpers.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

/// This adds the examples module entries to the Python module
namespace ActsPython {

void addFramework(py::module& mex);
void addAmbiguityResolution(py::module& mex);
void addGenerators(py::module& mex);
void addExampleAlgorithms(py::module& mex);
void addDispatchAlgorithms(py::module& mex);
void addDetector(py::module& mex);
void addDigitization(py::module& mex);
void addMaterialMapping(py::module& mex);
void addPropagation(py::module& mex);
void addTrackFitting(py::module& mex);
void addTrackFinding(py::module& mex);
void addTruthTracking(py::module& mex);
void addVertexing(py::module& mex);
void addObj(py::module& mex);
void addInput(py::module& mex);
void addOutput(py::module& mex);
void addUtilities(py::module& mex);
}  // namespace ActsPython

PYBIND11_MODULE(ActsExamplesPythonBindings, mex) {
  using namespace ActsPython;

  mex.doc() = "Acts Examples";

  addFramework(mex);
  addAmbiguityResolution(mex);
  addGenerators(mex);
  addExampleAlgorithms(mex);
  addDispatchAlgorithms(mex);
  addDetector(mex);
  addDigitization(mex);
  addMaterialMapping(mex);
  addPropagation(mex);
  addTrackFitting(mex);
  addTrackFinding(mex);
  addTruthTracking(mex);
  addVertexing(mex);
  addObj(mex);
  addInput(mex);
  addOutput(mex);
  addUtilities(mex);
}
