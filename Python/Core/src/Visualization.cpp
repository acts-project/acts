// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/AnyTrackProxy.hpp"
#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Visualization/EventDataView3D.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"
#include "Acts/Visualization/VisualizationBuffer.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <array>
#include <span>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace Acts;

namespace ActsPython {

py::array_t<int> colorsToNumpy(std::span<const Color> colors,
                               const py::object& base) {
  std::vector<std::size_t> shape{static_cast<std::size_t>(colors.size()), 3};

  std::vector<std::size_t> strides{static_cast<std::size_t>(sizeof(Color)),
                                   static_cast<std::size_t>(sizeof(int))};

  return py::array_t<int>(shape, strides, colors.data()->rgb.data(), base);
}

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

  py::class_<VisualizationBuffer, IVisualization3D>(m, "VisualizationBuffer")
      .def(py::init<unsigned int, double>(), py::arg("prec") = 4u,
           py::arg("scale") = 1)
      .def_property_readonly("surfaces", &VisualizationBuffer::surfaces)
      .def_property_readonly("segments", &VisualizationBuffer::segments)
      .def_property_readonly("vertices", &VisualizationBuffer::vertices3D)
      .def_property_readonly(
          "faceColors",
          [](const VisualizationBuffer& self) {
            return colorsToNumpy(std::span(self.faceColors()), py::cast(&self));
          })
      .def_property_readonly(
          "lineColor",
          [](const VisualizationBuffer& self) {
            return colorsToNumpy(std::span(self.lineColors()), py::cast(&self));
          })
      .def_property_readonly("lineThickness",
                             &VisualizationBuffer::lineThickness);

  py::class_<EventDataView3D>(m, "EventDataView3D")
      .def_static("drawTrack",
                  [](IVisualization3D& helper, const AnyConstTrackProxy& track,
                     const GeometryContext& gctx) {
                    EventDataView3D::drawTrack(helper, track, gctx);
                  });

  py::class_<GeometryObject>(m, "GeometryObject");
  m.def("default_geometry_coloring", &defaultGeometryColoring);
}
}  // namespace ActsPython
