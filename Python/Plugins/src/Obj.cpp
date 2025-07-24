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
#include "ActsPython/Utilities/Context.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Acts;

namespace ActsPython {
void addObj(Context& ctx) {
  // Register the obj module
  auto& p = ctx.get("plugins");
  auto obj = p.def_submodule("obj");

  {
    /// Write a collection of surfaces to an '.obj' file
    ///
    /// @param surfaces is the collection of surfaces
    /// @param viewContext is the geometry context
    /// @param viewRgb is the color of the surfaces
    /// @param viewSegments is the number of segments to approximate a quarter of a circle
    /// @param fileName is the path to the output file
    ///
    obj.def("writeSurfacesObj",
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
    obj.def("writeVolumesObj",
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
    obj.def("writeVolumesSurfacesObj",
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

  py::class_<ObjVisualization3D, IVisualization3D>(obj, "ObjVisualization3D")
      .def(py::init<>());
}
}  // namespace ActsPython
