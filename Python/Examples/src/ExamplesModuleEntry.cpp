// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/ActsVersion.hpp"
#include "ActsPython/Utilities/Helpers.hpp"

#include <tuple>
#include <unordered_map>

#include <pybind11/detail/common.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pyerrors.h>

namespace py = pybind11;
using namespace ActsPython;

namespace ActsPython {
void addFramework(py::module_& mex);
void addConcretePropagation(py::module_& mex);
void addDigitization(py::module_& mex);
void addAlignment(py::module_& mex);
void addMaterialMapping(py::module_& mex);
void addOutput(py::module_& mex);
void addDetector(py::module_& mex);
void addExampleAlgorithms(py::module_& mex);
void addInput(py::module_& mex);
void addGenerators(py::module_& mex);
void addTruthTracking(py::module_& mex);
void addTrackFitting(py::module_& mex);
void addTrackFinding(py::module_& mex);
void addVertexing(py::module_& mex);
void addAmbiguityResolution(py::module_& mex);
void addUtilities(py::module_& mex);
void addObj(py::module_& mex);

PYBIND11_MODULE(ActsExamplesPythonBindings, mex) {
  mex.doc() = "Acts Examples Module";

  addFramework(mex);
  addOutput(mex);
  addConcretePropagation(mex);
  addDigitization(mex);
  addAlignment(mex);
  addMaterialMapping(mex);
  addDetector(mex);
  addExampleAlgorithms(mex);
  addInput(mex);
  addGenerators(mex);
  addTruthTracking(mex);
  addTrackFitting(mex);
  addTrackFinding(mex);
  addVertexing(mex);
  addAmbiguityResolution(mex);
  addUtilities(mex);
  addObj(mex);
}

}  // namespace ActsPython
