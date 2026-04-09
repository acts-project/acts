// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/ProjectedVisualization.hpp"
#include "Acts/Visualization/ViewConfig.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace Acts;

namespace ActsPython {

/// @brief Add visualization bindings to the given module.
/// @param m    The module to add the bindings to.
void addVisualization(py::module& m) {
  {
    auto c = py::class_<ViewConfig>(m, "ViewConfig").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, visible, color, offset, lineThickness,
                       surfaceThickness, quarterSegments, triangulate,
                       outputName);

    patchKwargsConstructor(c);

    py::class_<Color>(m, "Color")
        .def(py::init<>())
        .def(py::init<int, int, int>())
        .def(py::init<double, double, double>())
        .def(py::init<std::string_view>())
        .def_readonly("rgb", &Color::rgb);
  }

  py::class_<IVisualization3D>(m, "IVisualization3D")
      .def("write", py::overload_cast<const std::filesystem::path&>(
                        &IVisualization3D::write, py::const_));

  py::class_<ObjVisualization3D, IVisualization3D>(m, "ObjVisualization3D")
      .def(py::init<unsigned int, double>(), py::arg("prec") = 4u,
           py::arg("scale") = 1.)
      .def("write",
           py::overload_cast<const std::filesystem::path&>(
               &ObjVisualization3D::write, py::const_),
           py::arg("path"))
      .def("clear", &ObjVisualization3D::clear)
      .def(
          "object",
          [](ObjVisualization3D& self, const std::string& name) {
            self.object(name);
          },
          py::arg("name"));

  py::class_<ProjectedVisualization, IVisualization3D>(m,
                                                       "ProjectedVisualization")
      .def(py::init(
               [](const std::vector<std::pair<
                      std::string, ProjectedVisualization::ProjectionFunction>>&
                      projections) {
                 std::vector<ProjectedVisualization::Projection> converted{
                     projections.begin(), projections.end()};
                 return ProjectedVisualization(std::move(converted));
               }),
           py::arg("projections"))
      .def("write",
           py::overload_cast<const std::filesystem::path&>(
               &ProjectedVisualization::write, py::const_),
           py::arg("path"))
      .def("vertex", &ProjectedVisualization::vertex, py::arg("vtx"),
           py::arg("color") = IVisualization3D::s_defaultColor)
      .def("face", &ProjectedVisualization::face, py::arg("vtxs"),
           py::arg("color") = IVisualization3D::s_defaultColor)
      .def("faces", &ProjectedVisualization::faces, py::arg("vtxs"),
           py::arg("faces"),
           py::arg("color") = IVisualization3D::s_defaultColor)
      .def("line", &ProjectedVisualization::line, py::arg("a"), py::arg("b"),
           py::arg("color") = IVisualization3D::s_defaultColor)
      .def("clear", &ProjectedVisualization::clear)
      .def("object", &ProjectedVisualization::object, py::arg("name"))
      .def_property_readonly("projections",
                             &ProjectedVisualization::projections)
      .def_property_readonly("projectedFaces",
                             &ProjectedVisualization::projectedFaces)
      .def_property_readonly("projectedLines",
                             &ProjectedVisualization::projectedLines)
      .def_property_readonly("projectedVertices",
                             &ProjectedVisualization::projectedVertices);

  m.def("projectToXY", &projectToXY, py::arg("position"));
  m.def("projectToZPhi", &projectToZPhi, py::arg("position"));
  m.def("projectToZR", &projectToZR, py::arg("position"));
}
}  // namespace ActsPython
