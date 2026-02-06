// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/ActsVersion.hpp"
#include "ActsPython/Utilities/Helpers.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

/// This adds the core module entries to the Python module
/// There is a certain order necessary as py::class_ definitions
/// need to be registered before they can be used in other modules.
namespace ActsPython {

void addDefinitions(py::module_& m);
void addMagneticField(py::module_& m);
void addUtilities(py::module_& m);
void addVisualization(py::module_& m);

void addMaterial(py::module_& m);
void addSurfaces(py::module_& m);
void addGeometry(py::module_& m);
void addGeometryGen1(py::module_& m);
void addGeometryGen3(py::module_& m);
void addNavigation(py::module_& m);
void addPropagation(py::module_& m);
void addSeeding(py::module_& mt);
void addTrackFinding(py::module_& m);

}  // namespace ActsPython

namespace py = pybind11;

PYBIND11_MODULE(ActsPythonBindings, m) {
  using namespace ActsPython;

  m.doc() = "Acts";
  m.attr("__version__") =
      std::tuple{Acts::VersionMajor, Acts::VersionMinor, Acts::VersionPatch};

  {
    auto mv = m.def_submodule("version");

    mv.attr("major") = Acts::VersionMajor;
    mv.attr("minor") = Acts::VersionMinor;
    mv.attr("patch") = Acts::VersionPatch;

    mv.attr("commit_hash") = std::string{Acts::CommitHash.value_or("UNKNOWN")};
    mv.attr("commit_hash_short") =
        std::string{Acts::CommitHashShort.value_or("UNKNOWN")};
  }

  addDefinitions(m);
  addMagneticField(m);
  addMaterial(m);
  addUtilities(m);
  addVisualization(m);

  addSurfaces(m);
  addGeometry(m);
  addGeometryGen1(m);
  addGeometryGen3(m);
  addNavigation(m);
  addPropagation(m);
  addSeeding(m);
  addTrackFinding(m);
}
