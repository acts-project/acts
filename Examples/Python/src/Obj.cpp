// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/IVisualization3D.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Detector/DetectorVolume.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Plugins/Python/Utilities.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Visualization/GeometryView3D.hpp>
#include <Acts/Visualization/ObjVisualization3D.hpp>
#include <Acts/Visualization/ViewConfig.hpp>

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;

namespace Acts::Python {
void addObj(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    py::class_<ViewConfig>(m, "ViewConfig")
        .def_readwrite("visible", &ViewConfig::visible)
        .def_readwrite("color", &ViewConfig::color)
        .def_readwrite("offset", &ViewConfig::offset)
        .def_readwrite("lineThickness", &ViewConfig::lineThickness)
        .def_readwrite("surfaceThickness", &ViewConfig::surfaceThickness)
        .def_readwrite("nSegments", &ViewConfig::nSegments)
        .def_readwrite("triangulate", &ViewConfig::triangulate)
        .def_readwrite("outputName", &ViewConfig::outputName);

    py::class_<Color>(m, "Color")
        .def(py::init<>())
        .def(py::init<int, int, int>())
        .def(py::init<double, double, double>())
        .def(py::init<std::string>())
        .def_readonly("rgb", &Color::rgb);
  }

  {
    /// Write a collection of surfaces to an '.obj' file
    ///
    /// @param surfaces is the collection of surfaces
    /// @param viewContext is the geometry context
    /// @param viewRgb is the color of the surfaces
    /// @param viewSegements is the number of segments to approximate a full circle
    /// @param fileName is the path to the output file
    ///
    mex.def("writeSurfacesObj",
            [](const std::vector<std::shared_ptr<Surface>>& surfaces,
               const GeometryContext& viewContext, const ViewConfig& viewConfig,
               const std::string& fileName) {
              Acts::GeometryView3D view3D;
              Acts::ObjVisualization3D obj;

              for (const auto& surface : surfaces) {
                view3D.drawSurface(obj, *surface, viewContext,
                                   Acts::Transform3::Identity(), viewConfig);
              }
              obj.write(fileName);
            });
    mex.def("writeVolumesObj",
            [](const std::vector<std::shared_ptr<Experimental::DetectorVolume>>&
                   Volumes,
               const GeometryContext& viewContext, const ViewConfig& viewConfig,
               const std::string& fileName) {
              Acts::GeometryView3D view3D;
              Acts::ObjVisualization3D obj;

              for (const auto& volume : Volumes) {
                view3D.drawDetectorVolume(obj, *volume, viewContext,
                                          Acts::Transform3::Identity(),
                                          viewConfig);
              }
              obj.write(fileName);
            });
    mex.def("writeVolumesSurfacesObj",
            [](const std::vector<std::shared_ptr<Surface>>& surfaces,
               const std::vector<std::shared_ptr<Experimental::DetectorVolume>>&
                   Volumes,
               const GeometryContext& viewContext, const ViewConfig& viewConfig,
               const std::string& fileName) {
              Acts::GeometryView3D view3D;
              Acts::ObjVisualization3D obj;

              for (const auto& volume : Volumes) {
                view3D.drawDetectorVolume(obj, *volume, viewContext,
                                          Acts::Transform3::Identity(),
                                          viewConfig);
              }
              for (const auto& surface : surfaces) {
                view3D.drawSurface(obj, *surface, viewContext,
                                   Acts::Transform3::Identity(), viewConfig);
              }
              obj.write(fileName);
            });
  }

  py::class_<IVisualization3D>(m, "IVisualization3D")
      .def("write", py::overload_cast<const std::filesystem::path&>(
                        &IVisualization3D::write, py::const_));

  py::class_<ObjVisualization3D<double>, IVisualization3D>(m,
                                                           "ObjVisualization3D")
      .def(py::init<>())
      .def("write", py::overload_cast<const std::filesystem::path&>(
                        &ObjVisualization3D<double>::write, py::const_));
}
}  // namespace Acts::Python
