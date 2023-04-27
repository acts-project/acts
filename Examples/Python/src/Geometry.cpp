// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/ProtoDetector.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace {
struct GeometryIdentifierHookBinding : public Acts::GeometryIdentifierHook {
  py::object callable;

  Acts::GeometryIdentifier decorateIdentifier(
      Acts::GeometryIdentifier identifier,
      const Acts::Surface& surface) const override {
    return callable(identifier, surface.getSharedPtr())
        .cast<Acts::GeometryIdentifier>();
  }
};
}  // namespace

namespace Acts::Python {
void addGeometry(Context& ctx) {
  auto m = ctx.get("main");
  {
    py::class_<Acts::Surface, std::shared_ptr<Acts::Surface>>(m, "Surface")
        .def("geometryId",
             [](Acts::Surface& self) { return self.geometryId(); })
        .def("center",
             [](Acts::Surface& self) {
               return self.center(Acts::GeometryContext{});
             })
        .def("type", [](Acts::Surface& self) { return self.type(); });
  }

  {
    py::enum_<Acts::Surface::SurfaceType>(m, "SurfaceType")
        .value("Cone", Acts::Surface::SurfaceType::Cone)
        .value("Cylinder", Acts::Surface::SurfaceType::Cylinder)
        .value("Disc", Acts::Surface::SurfaceType::Disc)
        .value("Perigee", Acts::Surface::SurfaceType::Perigee)
        .value("Plane", Acts::Surface::SurfaceType::Plane)
        .value("Straw", Acts::Surface::SurfaceType::Straw)
        .value("Curvilinear", Acts::Surface::SurfaceType::Curvilinear)
        .value("Other", Acts::Surface::SurfaceType::Other);
  }

  {
    py::class_<Acts::TrackingGeometry, std::shared_ptr<Acts::TrackingGeometry>>(
        m, "TrackingGeometry")
        .def("visitSurfaces",
             [](Acts::TrackingGeometry& self, py::function& func) {
               self.visitSurfaces(func);
             });
  }

  {
    py::class_<Acts::Volume, std::shared_ptr<Acts::Volume>>(m, "Volume")
        .def_static(
            "makeCylinderVolume",
            [](double r, double halfZ) {
              auto bounds =
                  std::make_shared<Acts::CylinderVolumeBounds>(0, r, halfZ);
              return std::make_shared<Acts::Volume>(Transform3::Identity(),
                                                    bounds);
            },
            "r"_a, "halfZ"_a);
  }

  {
    py::class_<Acts::GeometryIdentifierHook,
               std::shared_ptr<Acts::GeometryIdentifierHook>>(
        m, "GeometryIdentifierHook")
        .def(py::init([](py::object callable) {
          auto hook = std::make_shared<GeometryIdentifierHookBinding>();
          hook->callable = callable;
          return hook;
        }));
  }

  {
    py::class_<Acts::Experimental::Detector,
               std::shared_ptr<Acts::Experimental::Detector>>(m, "Detector");
  }

  {
    py::class_<Acts::Experimental::DetectorVolume,
               std::shared_ptr<Acts::Experimental::DetectorVolume>>(
        m, "DetectorVolume");
  }
}

}  // namespace Acts::Python
