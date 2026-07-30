// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/PyVisualization.hpp"
#include "Acts/Visualization/EventDataView3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"
#include "Acts/Geometry/GeometryObject.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>
#include <pybind11/functional.h>

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

  py::class_<PyVisualization, IVisualization3D>(m, "PyVisualization")
      .def(py::init<unsigned int, double>(), py::arg("prec") = 4u,
           py::arg("scale") = 1)
      .def_property_readonly("surfaces", &PyVisualization::surfaces)
      .def_property_readonly("segments", &PyVisualization::segments)
      .def_property_readonly("faceColors", [](const PyVisualization &self) {
        const auto &colors = self.faceColors();

        py::array_t<double> arr({colors.size(), size_t(3)});
        // one element in the array
        auto r = arr.mutable_unchecked<2>(); // direct access of elements without internal checking of dimensions

        for (size_t i = 0; i < colors.size(); ++i) {
            r(i, 0) = colors[i].rgb[0];
            r(i, 1) = colors[i].rgb[1];
            r(i, 2) = colors[i].rgb[2];
        }
        return arr;
      }
    )
    .def_property_readonly("lineColor", [](const PyVisualization &self){
      const auto &colors = self.lineColors();

      py::array_t<double> arr({colors.size(), size_t(3)});
      auto r = arr.mutable_unchecked<2>();

      for (size_t i = 0; i < colors.size(); i++) {
          r(i, 0) = colors[i].rgb[0];
          r(i, 1) = colors[i].rgb[1];
          r(i, 2) = colors[i].rgb[2];
        }
      return arr;
      }
    )
    .def_property_readonly("lineThickness", &PyVisualization::lineThickness);

  py::class_<EventDataView3D>(m, "EventDataView3D")
    .def_static("drawTrack",
      [] (IVisualization3D& helper, const Acts::EventDataView3D::ConstTrackProxy& track){
        EventDataView3D::drawTrack(helper, track,  GeometryContext::dangerouslyDefaultConstruct());
      });

  py::class_<GeometryObject>(m, "GeometryObject");
  m.def("makeDefaultColoringFunction", &makeDefaultColoringFunction);
  //m.def("defaultColoring", &defaultColoring);

}
}  // namespace ActsPython
