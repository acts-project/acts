// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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

  // OBJ based visualiation - no dependency
  {
    py::class_<ObjVisualization3D, IVisualization3D>(m, "ObjVisualization3D")
        .def(py::init<>());

    /// Write a collection of surfaces to an '.obj' file
    ///
    /// @param surfaces is the collection of surfaces
    /// @param viewContext is the geometry context
    /// @param viewRgb is the color of the surfaces
    /// @param viewSegments is the number of segments to approximate a quarter of a circle
    /// @param fileName is the path to the output file
    ///
    m.def("writeSurfacesObj",
          [](const std::vector<std::shared_ptr<Surface>>& surfaces,
             const GeometryContext& viewContext, const ViewConfig& viewConfig,
             const std::string& fileName) {
            GeometryView3D view3D;
            ObjVisualization3D obj;

            for (const auto& surface : surfaces) {
              view3D.drawSurface(obj, *surface, viewContext,
                                 Transform3::Identity(), viewConfig);
            }
            obj.write(fileName);
          });
    m.def("writeVolumesObj",
          [](const std::vector<std::shared_ptr<Experimental::DetectorVolume>>&
                 Volumes,
             const GeometryContext& viewContext, const ViewConfig& viewConfig,
             const std::string& fileName) {
            GeometryView3D view3D;
            ObjVisualization3D obj;

            for (const auto& volume : Volumes) {
              view3D.drawDetectorVolume(obj, *volume, viewContext,
                                        Transform3::Identity(), viewConfig);
            }
            obj.write(fileName);
          });
    m.def("writeVolumesSurfacesObj",
          [](const std::vector<std::shared_ptr<Surface>>& surfaces,
             const std::vector<std::shared_ptr<Experimental::DetectorVolume>>&
                 Volumes,
             const GeometryContext& viewContext, const ViewConfig& viewConfig,
             const std::string& fileName) {
            GeometryView3D view3D;
            ObjVisualization3D obj;

            for (const auto& volume : Volumes) {
              view3D.drawDetectorVolume(obj, *volume, viewContext,
                                        Transform3::Identity(), viewConfig);
            }
            for (const auto& surface : surfaces) {
              view3D.drawSurface(obj, *surface, viewContext,
                                 Transform3::Identity(), viewConfig);
            }
            obj.write(fileName);
          });
  }
}
}  // namespace ActsPython
