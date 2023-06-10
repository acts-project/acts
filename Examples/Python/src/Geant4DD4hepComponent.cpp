// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/DDG4/DDG4DetectorConstruction.hpp"

#include <functional>
#include <memory>

#include <G4VUserDetectorConstruction.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;
using namespace Acts::Python;

PYBIND11_MODULE(ActsPythonBindingsDDG4, m) {
  py::module_::import("acts.ActsPythonBindingsGeant4");

  py::class_<DDG4DetectorConstructionFactory, DetectorConstructionFactory,
             std::shared_ptr<DDG4DetectorConstructionFactory>>(
      m, "DDG4DetectorConstructionFactory")
      .def(py::init<std::shared_ptr<DD4hep::DD4hepDetector>>());
}
