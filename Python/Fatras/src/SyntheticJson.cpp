// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// @file
/// Reading and writing a synthetic detector, which needs `ActsFatrasJson` and
/// so only exists with the JSON plugin built. A detector as a file is how one
/// is meant to be had, so this is where the bindings a caller reaches for live.

#include "ActsFatras/Json/DataDirectory.hpp"
#include "ActsFatras/Json/DetectorDescriptionJsonConverter.hpp"
#include "ActsFatras/Json/EventConfigJsonConverter.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/EventConfig.hpp"

#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace ActsPython {

/// Bind them onto the synthetic module.
/// @param m the module to add them to
void addFatrasSyntheticJson(py::module_& m) {
  using namespace ActsFatras::Synthetic;

  m.def(
      "readDetectorDescription",
      [](const std::string& path) { return readDetectorDescription(path); },
      py::arg("path"));
  m.def(
      "writeDetectorDescription",
      [](const std::string& path, const DetectorDescription& description) {
        writeDetectorDescription(path, description);
      },
      py::arg("path"), py::arg("description"));
  m.def(
      "readMaterialDecoration",
      [](const std::string& path) { return readMaterialDecoration(path); },
      py::arg("path"));
  m.def(
      "writeMaterialDecoration",
      [](const std::string& path, const MaterialDecoration& decoration) {
        writeMaterialDecoration(path, decoration);
      },
      py::arg("path"), py::arg("decoration"));
  m.def(
      "readEventConfig",
      [](const std::string& path) { return readEventConfig(path); },
      py::arg("path"));
  m.def(
      "writeEventConfig",
      [](const std::string& path, const EventConfig& config) {
        writeEventConfig(path, config);
      },
      py::arg("path"), py::arg("config"));
  m.def(
      "dataPath",
      [](const std::string& name) {
        return ActsFatras::dataPath(name).string();
      },
      py::arg("name"));
}

}  // namespace ActsPython
