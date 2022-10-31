// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetectorConstruction.hpp"

#include <functional>
#include <memory>

#include <G4VUserDetectorConstruction.hh>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

PYBIND11_MODULE(ActsPythonBindingsTelescopeDetectorG4, m) {
  py::module_::import("acts.ActsPythonBindingsGeant4");

  // This is the actual class we're binding
  py::class_<TelescopeDetectorConstruction, G4VUserDetectorConstruction>(
      m, "TelescopeDetectorConstructionImpl");

  // This is a python-only factory method that returns the above class.
  // We can apply a return value policy here so that python does NOT assume
  // ownership of the returned pointer, and it is safe to pass to G4
  m.def(
      "TelescopeDetectorConstruction",
      [](DD4hep::DD4hepGeometryService& geometrySvc) {
        return new TelescopeDetectorConstruction(*geometrySvc.lcdd());
      },
      py::return_value_policy::reference);
}
