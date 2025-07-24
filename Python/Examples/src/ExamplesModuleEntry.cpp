// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <pybind11/pybind11.h>

namespace py = pybind11;

/// This adds the python module entries for the Examples
namespace ActsPython {
// These will be added to `acts.examples` module
void addFramework(py::module_& mex);
void addDetector(py::module_& mex);
void addGenerators(py::module_& mex);
void addAmbiguityResolution(py::module_& mex);
void addMagneticFieldMaps(py::module_& mex);
void addDigitizationAlgorithms(py::module_& mex);
void addMaterialMappingAlgorithms(py::module_& mex);
void addPropagationAlgorithms(py::module_& mex);
void addTrackFindingAlgorithms(py::module_& mex);
void addTrackFittingAlgorithms(py::module_& mex);
void addTruthTrackingAlgorithms(py::module_& mex);
void addVertexingAlgorithms(py::module_& mex);
void addPrinterAlgorithms(py::module_& mex);
void addUtilityAlgorithms(py::module_& mex);
void addInput(py::module_& mex);
void addOutput(py::module_& mex);
// Plugin dependent components
void addFatrasAlgorithms(py::module_& mex);
void addGeoModelDetector(py::module_& mex);
void addTGeoDetector(py::module_& mex);
void addDD4hepDetector(py::module_& mex);
void addRootInput(py::module_& mex);
void addRootOutput(py::module_& mex);
void addJsonInputOutput(py::module_& mex);
void addSvgOutput(py::module_& mex);
void addTracccAlgorithms(py::module_& mex);
void addPythia8Generator(py::module_& mex);
void addTruthJetAlgorithms(py::module_& mex);
void addOnnxAlgorithms(py::module_& mex);

}  // namespace ActsPython

PYBIND11_MODULE(ActsExamplesPythonBindings, mex) {
  using namespace ActsPython;

  addFramework(mex);
  addDetector(mex);
  addGenerators(mex);
  addAmbiguityResolution(mex);
  addMagneticFieldMaps(mex);
  addDigitizationAlgorithms(mex);
  addMaterialMappingAlgorithms(mex);
  addPropagationAlgorithms(mex);
  addTrackFindingAlgorithms(mex);
  addTrackFittingAlgorithms(mex);
  addTruthTrackingAlgorithms(mex);
  addVertexingAlgorithms(mex);
  addPrinterAlgorithms(mex);
  addUtilityAlgorithms(mex);
  addInput(mex);
  addOutput(mex);
  // Plugin dependent components
  addFatrasAlgorithms(mex);
  addGeoModelDetector(mex);
  addTGeoDetector(mex);
  addDD4hepDetector(mex);
  addRootInput(mex);
  addRootOutput(mex);
  addJsonInputOutput(mex);
  addSvgOutput(mex);
  addTracccAlgorithms(mex);
  addPythia8Generator(mex);
  addTruthJetAlgorithms(mex);
  addOnnxAlgorithms(mex);
}
