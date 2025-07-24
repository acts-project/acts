// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"

#include <pybind11/eval.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Acts;

namespace ActsPython {
// This adds the definitions from Core/Surfaces to the python module
/// @param m is the pybind11 core module
void addSurfaces(py::module_& m) {
  {
    py::class_<Surface, std::shared_ptr<Surface>>(m, "Surface")
        // Can't bind directly because GeometryObject is virtual base of Surface
        .def_property_readonly(
            "geometryId", [](const Surface& self) { return self.geometryId(); })
        .def("center", &Surface::center)
        .def_property_readonly("type", &Surface::type)
        .def("visualize", &Surface::visualize)
        .def_property_readonly("surfaceMaterial",
                               &Surface::surfaceMaterialSharedPtr);
  }

  {
    py::enum_<Surface::SurfaceType>(m, "SurfaceType")
        .value("Cone", Surface::SurfaceType::Cone)
        .value("Cylinder", Surface::SurfaceType::Cylinder)
        .value("Disc", Surface::SurfaceType::Disc)
        .value("Perigee", Surface::SurfaceType::Perigee)
        .value("Plane", Surface::SurfaceType::Plane)
        .value("Straw", Surface::SurfaceType::Straw)
        .value("Curvilinear", Surface::SurfaceType::Curvilinear)
        .value("Other", Surface::SurfaceType::Other);
  }

  // Add the surface hierarchy map
  using SurfaceHierarchyMap =
      Acts::GeometryHierarchyMap<std::shared_ptr<Acts::Surface>>;

  py::class_<SurfaceHierarchyMap, std::shared_ptr<SurfaceHierarchyMap>>(
      m, "SurfaceHierarchyMap");
}

}  // namespace ActsPython
