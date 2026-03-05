// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
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
    py::class_<SurfaceBounds, std::shared_ptr<SurfaceBounds>>(m,
                                                              "SurfaceBounds");
  }

  {
    py::class_<CylinderBounds, SurfaceBounds, std::shared_ptr<CylinderBounds>>(
        m, "CylinderBounds")
        .def(py::init<double, double>())
        .def(py::init(
            [](const std::array<double, CylinderBounds::eSize>& values) {
              return CylinderBounds(values);
            }))
        .def("__getitem__",
             [](const CylinderBounds& self, int index) {
               return self.get(CylinderBounds::BoundValues(index));
             })
        .def("__str__", [](const CylinderBounds& self) {
          std::ostringstream oss;
          oss << self;
          return oss.str();
        });
  }

  {
    py::class_<DiscBounds, SurfaceBounds, std::shared_ptr<DiscBounds>>(
        m, "DiscBounds");

    py::class_<AnnulusBounds, DiscBounds, std::shared_ptr<AnnulusBounds>>(
        m, "AnnulusBounds")
        .def(py::init<double, double, double, double>())
        .def(py::init(
            [](const std::array<double, AnnulusBounds::eSize>& values) {
              return AnnulusBounds(values);
            }))
        .def("__getitem__",
             [](const AnnulusBounds& self, int index) {
               return self.get(AnnulusBounds::BoundValues(index));
             })
        .def("__str__", [](const AnnulusBounds& self) {
          std::ostringstream oss;
          oss << self;
          return oss.str();
        });

    py::class_<RadialBounds, DiscBounds, std::shared_ptr<RadialBounds>>(
        m, "RadialBounds")
        .def(py::init<double, double>())
        .def(
            py::init([](const std::array<double, RadialBounds::eSize>& values) {
              return RadialBounds(values);
            }))
        .def("__getitem__",
             [](const RadialBounds& self, int index) {
               return self.get(RadialBounds::BoundValues(index));
             })
        .def("__str__", [](const RadialBounds& self) {
          std::ostringstream oss;
          oss << self;
          return oss.str();
        });
  }

  {
    py::class_<LineBounds, SurfaceBounds, std::shared_ptr<LineBounds>>(
        m, "LineBounds")
        .def(py::init<double, double>())
        .def(py::init([](const std::array<double, LineBounds::eSize>& values) {
          return LineBounds(values);
        }))
        .def("__getitem__",
             [](const LineBounds& self, int index) {
               return self.get(LineBounds::BoundValues(index));
             })
        .def("__str__", [](const LineBounds& self) {
          std::ostringstream oss;
          oss << self;
          return oss.str();
        });
  }

  {
    py::class_<PlanarBounds, SurfaceBounds, std::shared_ptr<PlanarBounds>>(
        m, "PlanarBounds");

    py::class_<RectangleBounds, PlanarBounds, std::shared_ptr<RectangleBounds>>(
        m, "RectangleBounds")
        .def(py::init<double, double>())
        .def(py::init(
            [](const std::array<double, RectangleBounds::eSize>& values) {
              return RectangleBounds(values);
            }))
        .def("__getitem__",
             [](const RectangleBounds& self, int index) {
               return self.get(RectangleBounds::BoundValues(index));
             })
        .def("__str__", [](const RectangleBounds& self) {
          std::ostringstream oss;
          oss << self;
          return oss.str();
        });

    py::class_<TrapezoidBounds, PlanarBounds, std::shared_ptr<TrapezoidBounds>>(
        m, "TrapezoidBounds")
        .def(py::init<double, double, double, double>())
        .def(py::init(
            [](const std::array<double, TrapezoidBounds::eSize>& values) {
              return TrapezoidBounds(values);
            }))
        .def("__getitem__",
             [](const TrapezoidBounds& self, int index) {
               return self.get(TrapezoidBounds::BoundValues(index));
             })
        .def("__str__", [](const TrapezoidBounds& self) {
          std::ostringstream oss;
          oss << self;
          return oss.str();
        });
  }

  {
    py::class_<Surface, std::shared_ptr<Surface>>(m, "Surface")
        // Can't bind directly because GeometryObject is virtual base of Surface
        .def_property_readonly(
            "geometryId", [](const Surface& self) { return self.geometryId(); })
        .def("assignGeometryId",
             [](Surface& self, const GeometryIdentifier& id) {
               self.assignGeometryId(id);
             })
        .def("center", &Surface::center)
        .def_property_readonly("type", &Surface::type)
        .def("visualize", &Surface::visualize)
        .def_property_readonly("surfaceMaterial",
                               &Surface::surfaceMaterialSharedPtr)
        .def("createCylinder",
             [](const Transform3& transform,
                const std::shared_ptr<const CylinderBounds>& bounds) {
               return Surface::makeShared<CylinderSurface>(transform, bounds);
             })
        .def("createDisc",
             [](const Transform3& transform,
                const std::shared_ptr<const DiscBounds>& bounds) {
               return Surface::makeShared<DiscSurface>(transform, bounds);
             })
        .def("createDisc",
             [](const std::shared_ptr<const DiscBounds>& bounds,
                const SurfacePlacementBase& detelement) {
               return Surface::makeShared<DiscSurface>(bounds, detelement);
             })
        .def("createStraw",
             [](const Transform3& transform,
                const std::shared_ptr<const LineBounds>& bounds) {
               return Surface::makeShared<StrawSurface>(transform, bounds);
             })
        .def("createStraw",
             [](const std::shared_ptr<const LineBounds>& bounds,
                const SurfacePlacementBase& detelement) {
               return Surface::makeShared<StrawSurface>(bounds, detelement);
             })
        .def("createPerigee",
             [](const Vector3& vertex) {
               return Surface::makeShared<PerigeeSurface>(vertex);
             })
        .def("createPlane",
             [](const Transform3& transform,
                const std::shared_ptr<const PlanarBounds>& bounds) {
               return Surface::makeShared<PlaneSurface>(transform, bounds);
             })
        .def("createPlane",
             [](const std::shared_ptr<const PlanarBounds>& pbounds,
                const SurfacePlacementBase& detelement) {
               return Surface::makeShared<PlaneSurface>(pbounds, detelement);
             });

    py::class_<CylinderSurface, Surface, std::shared_ptr<CylinderSurface>>(
        m, "CylinderSurface");

    py::class_<DiscSurface, Surface, std::shared_ptr<DiscSurface>>(
        m, "DiscSurface");

    py::class_<PlaneSurface, Surface, std::shared_ptr<PlaneSurface>>(
        m, "PlaneSurface");

    py::class_<PerigeeSurface, Surface, std::shared_ptr<PerigeeSurface>>(
        m, "PerigeeSurface");

    py::class_<StrawSurface, Surface, std::shared_ptr<StrawSurface>>(
        m, "StrawSurface");
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
      m, "SurfaceHierarchyMap")
      .def(py::init<>())
      .def(py::init<std::vector<SurfaceHierarchyMap::InputElement>>())
      .def("__len__",
           [](const SurfaceHierarchyMap& self) { return self.size(); })
      .def("__getitem__",
           [](const SurfaceHierarchyMap& self, std::size_t index) {
             return self.valueAt(index);
           });
}

}  // namespace ActsPython
