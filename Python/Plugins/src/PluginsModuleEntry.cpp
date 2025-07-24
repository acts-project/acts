// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <pybind11/pybind11.h>

namespace py = pybind11;

/// This adds the plugins module entries to the Python module
/// There is a certain order necessary as py::class_ definitions
/// need to be registered before they can be used in other modules.
namespace ActsPython {

void addCovfie(py::module_& p);
void addDetray(py::module_& p);
void addHashing(py::module_& p);
void addObj(py::module_& p);
void addJson(py::module_& p);
void addSvg(py::module_& p);
void addFpeMonitoring(py::module_& p);
// Geometry extensions
void addGeoModel(py::module_& p);
void addTGeo(py::module_& p);

}  // namespace ActsPython

PYBIND11_MODULE(ActsPluginsPythonBindings, p) {
  using namespace ActsPython;
  addCovfie(p);
  addDetray(p);
  addHashing(p);
  addObj(p);
  addJson(p);
  addSvg(p);
  addFpeMonitoring(p);
  // Geometry extensions
  addGeoModel(p);
  addTGeo(p);
}
