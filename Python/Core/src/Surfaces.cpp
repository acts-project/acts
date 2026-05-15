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
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/LineSurface.hpp"
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

#include <sstream>

#include <pybind11/eval.h>
#include <pybind11/operators.h>
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
    py::class_<BoundaryTolerance>(m, "BoundaryTolerance")
        .def_static("infinite", &BoundaryTolerance::Infinite)
        .def_static("none", &BoundaryTolerance::None)
        .def_static("absoluteEuclidean", &BoundaryTolerance::AbsoluteEuclidean)
        .def("isInfinite", &BoundaryTolerance::isInfinite)
        .def("isNone", &BoundaryTolerance::isNone)
        .def("hasAbsoluteEuclidean", &BoundaryTolerance::hasAbsoluteEuclidean)
        .def("hasChi2Bound", &BoundaryTolerance::hasChi2Bound)
        .def("hasChi2Cartesian", &BoundaryTolerance::hasChi2Cartesian);
  }

  {
    py::class_<SurfaceBounds, std::shared_ptr<SurfaceBounds>>(m,
                                                              "SurfaceBounds")
        .def_property_readonly("type", &SurfaceBounds::type)
        .def("isCartesian", &SurfaceBounds::isCartesian)
        .def("boundToCartesianJacobian",
             &SurfaceBounds::boundToCartesianJacobian)
        .def("boundToCartesianMetric", &SurfaceBounds::boundToCartesianMetric)
        .def("values", &SurfaceBounds::values)
        .def("inside", py::overload_cast<const Vector2&>(&SurfaceBounds::inside,
                                                         py::const_))
        .def("inside",
             py::overload_cast<const Vector2&, const BoundaryTolerance&>(
                 &SurfaceBounds::inside, py::const_))
        .def("closestPoint", &SurfaceBounds::closestPoint)
        .def("distance", &SurfaceBounds::distance)
        .def("center", &SurfaceBounds::center)
        .def(py::self == py::self)
        .def("__str__", [](const SurfaceBounds& self) {
          std::ostringstream oss;
          oss << self;
          return oss.str();
        });

    py::enum_<SurfaceBounds::BoundsType>(m, "SurfaceBoundsType")
        .value("Cone", SurfaceBounds::BoundsType::eCone)
        .value("Cylinder", SurfaceBounds::BoundsType::eCylinder)
        .value("Diamond", SurfaceBounds::BoundsType::eDiamond)
        .value("Disc", SurfaceBounds::BoundsType::eDisc)
        .value("Ellipse", SurfaceBounds::BoundsType::eEllipse)
        .value("Line", SurfaceBounds::BoundsType::eLine)
        .value("Rectangle", SurfaceBounds::BoundsType::eRectangle)
        .value("Trapezoid", SurfaceBounds::BoundsType::eTrapezoid)
        .value("Triangle", SurfaceBounds::BoundsType::eTriangle)
        .value("DiscTrapezoid", SurfaceBounds::BoundsType::eDiscTrapezoid)
        .value("ConvexPolygon", SurfaceBounds::BoundsType::eConvexPolygon)
        .value("Annulus", SurfaceBounds::BoundsType::eAnnulus)
        .value("Boundless", SurfaceBounds::BoundsType::eBoundless)
        .value("Other", SurfaceBounds::BoundsType::eOther);
  }

  {
    py::class_<CylinderBounds, SurfaceBounds, std::shared_ptr<CylinderBounds>>(
        m, "CylinderBounds")
        .def(py::init<double, double>())
        .def(py::init(
            [](const std::array<double, CylinderBounds::eSize>& values) {
              return CylinderBounds(values);
            }))
        .def("__len__",
             [](const CylinderBounds&) { return CylinderBounds::eSize; })
        .def("__getitem__",
             [](const CylinderBounds& self, int index) {
               if (index < 0 || index >= CylinderBounds::eSize) {
                 throw py::index_error("CylinderBounds index out of range");
               }
               return self.get(CylinderBounds::BoundValues(index));
             })
        .def(
            "get",
            [](const CylinderBounds& self, CylinderBounds::BoundValues bValue) {
              return self.get(bValue);
            })
        .def("__str__", [](const CylinderBounds& self) {
          std::ostringstream oss;
          oss << self;
          return oss.str();
        });

    py::enum_<CylinderBounds::BoundValues>(m, "CylinderBoundsValue")
        .value("R", CylinderBounds::BoundValues::eR)
        .value("HalfLengthZ", CylinderBounds::BoundValues::eHalfLengthZ)
        .value("HalfPhiSector", CylinderBounds::BoundValues::eHalfPhiSector)
        .value("AveragePhi", CylinderBounds::BoundValues::eAveragePhi)
        .value("BevelMinZ", CylinderBounds::BoundValues::eBevelMinZ)
        .value("BevelMaxZ", CylinderBounds::BoundValues::eBevelMaxZ);
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
        .def("__len__",
             [](const AnnulusBounds&) { return AnnulusBounds::eSize; })
        .def("__getitem__",
             [](const AnnulusBounds& self, int index) {
               if (index < 0 || index >= AnnulusBounds::eSize) {
                 throw py::index_error("AnnulusBounds index out of range");
               }
               return self.get(AnnulusBounds::BoundValues(index));
             })
        .def("get",
             [](const AnnulusBounds& self, AnnulusBounds::BoundValues bValue) {
               return self.get(bValue);
             })
        .def("__str__", [](const AnnulusBounds& self) {
          std::ostringstream oss;
          oss << self;
          return oss.str();
        });

    py::enum_<AnnulusBounds::BoundValues>(m, "AnnulusBoundsValue")
        .value("MinR", AnnulusBounds::BoundValues::eMinR)
        .value("MaxR", AnnulusBounds::BoundValues::eMaxR)
        .value("MinPhiRel", AnnulusBounds::BoundValues::eMinPhiRel)
        .value("MaxPhiRel", AnnulusBounds::BoundValues::eMaxPhiRel)
        .value("OriginX", AnnulusBounds::BoundValues::eOriginX)
        .value("OriginY", AnnulusBounds::BoundValues::eOriginY);

    py::class_<RadialBounds, DiscBounds, std::shared_ptr<RadialBounds>>(
        m, "RadialBounds")
        .def(py::init<double, double>())
        .def(
            py::init([](const std::array<double, RadialBounds::eSize>& values) {
              return RadialBounds(values);
            }))
        .def("__len__", [](const RadialBounds&) { return RadialBounds::eSize; })
        .def("__getitem__",
             [](const RadialBounds& self, int index) {
               if (index < 0 || index >= RadialBounds::eSize) {
                 throw py::index_error("RadialBounds index out of range");
               }
               return self.get(RadialBounds::BoundValues(index));
             })
        .def("get",
             [](const RadialBounds& self, RadialBounds::BoundValues bValue) {
               return self.get(bValue);
             })
        .def("__str__", [](const RadialBounds& self) {
          std::ostringstream oss;
          oss << self;
          return oss.str();
        });

    py::enum_<RadialBounds::BoundValues>(m, "RadialBoundsValue")
        .value("MinR", RadialBounds::BoundValues::eMinR)
        .value("MaxR", RadialBounds::BoundValues::eMaxR)
        .value("AveragePhi", RadialBounds::BoundValues::eAveragePhi)
        .value("HalfPhiSector", RadialBounds::BoundValues::eHalfPhiSector);
  }

  {
    py::class_<LineBounds, SurfaceBounds, std::shared_ptr<LineBounds>>(
        m, "LineBounds")
        .def(py::init<double, double>())
        .def(py::init([](const std::array<double, LineBounds::eSize>& values) {
          return LineBounds(values);
        }))
        .def("__len__", [](const LineBounds&) { return LineBounds::eSize; })
        .def("__getitem__",
             [](const LineBounds& self, int index) {
               if (index < 0 || index >= LineBounds::eSize) {
                 throw py::index_error("LineBounds index out of range");
               }
               return self.get(LineBounds::BoundValues(index));
             })
        .def("get",
             [](const LineBounds& self, LineBounds::BoundValues bValue) {
               return self.get(bValue);
             })
        .def("__str__", [](const LineBounds& self) {
          std::ostringstream oss;
          oss << self;
          return oss.str();
        });

    py::enum_<LineBounds::BoundValues>(m, "LineBoundsValue")
        .value("R", LineBounds::BoundValues::eR)
        .value("HalfLengthZ", LineBounds::BoundValues::eHalfLengthZ);
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
        .def("__len__",
             [](const RectangleBounds&) { return RectangleBounds::eSize; })
        .def("__getitem__",
             [](const RectangleBounds& self, int index) {
               if (index < 0 || index >= RectangleBounds::eSize) {
                 throw py::index_error("RectangleBounds index out of range");
               }
               return self.get(RectangleBounds::BoundValues(index));
             })
        .def("get",
             [](const RectangleBounds& self,
                RectangleBounds::BoundValues bValue) {
               return self.get(bValue);
             })
        .def("__str__", [](const RectangleBounds& self) {
          std::ostringstream oss;
          oss << self;
          return oss.str();
        });

    py::enum_<RectangleBounds::BoundValues>(m, "RectangleBoundsValue")
        .value("MinX", RectangleBounds::BoundValues::eMinX)
        .value("MinY", RectangleBounds::BoundValues::eMinY)
        .value("MaxX", RectangleBounds::BoundValues::eMaxX)
        .value("MaxY", RectangleBounds::BoundValues::eMaxY);

    py::class_<TrapezoidBounds, PlanarBounds, std::shared_ptr<TrapezoidBounds>>(
        m, "TrapezoidBounds")
        .def(py::init<double, double, double, double>())
        .def(py::init(
            [](const std::array<double, TrapezoidBounds::eSize>& values) {
              return TrapezoidBounds(values);
            }))
        .def("__len__",
             [](const TrapezoidBounds&) { return TrapezoidBounds::eSize; })
        .def("__getitem__",
             [](const TrapezoidBounds& self, int index) {
               if (index < 0 || index >= TrapezoidBounds::eSize) {
                 throw py::index_error("TrapezoidBounds index out of range");
               }
               return self.get(TrapezoidBounds::BoundValues(index));
             })
        .def("get",
             [](const TrapezoidBounds& self,
                TrapezoidBounds::BoundValues bValue) {
               return self.get(bValue);
             })
        .def("__str__", [](const TrapezoidBounds& self) {
          std::ostringstream oss;
          oss << self;
          return oss.str();
        });

    py::enum_<TrapezoidBounds::BoundValues>(m, "TrapezoidBoundsValue")
        .value("HalfLengthXnegY",
               TrapezoidBounds::BoundValues::eHalfLengthXnegY)
        .value("HalfLengthXposY",
               TrapezoidBounds::BoundValues::eHalfLengthXposY)
        .value("HalfLengthY", TrapezoidBounds::BoundValues::eHalfLengthY)
        .value("RotationAngle", TrapezoidBounds::BoundValues::eRotationAngle);
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
        .def("normal", &Surface::normal)
        .def("insideBounds", &Surface::insideBounds)
        .def("isOnSurface", &Surface::isOnSurface)
        .def("localToGlobal", &Surface::localToGlobal)
        .def("localToGlobalTransform", &Surface::localToGlobalTransform)
        .def("toString", &Surface::toString)
        .def_property_readonly("type", &Surface::type)
        .def_property_readonly("name", &Surface::name)
        .def_property_readonly("bounds", &Surface::bounds,
                               py::return_value_policy::reference_internal)
        .def_property_readonly("thickness", &Surface::thickness)
        .def_property_readonly("isSensitive", &Surface::isSensitive)
        .def_property_readonly("isAlignable", &Surface::isAlignable)
        .def("visualize", &Surface::visualize)
        .def_property_readonly("surfaceMaterial",
                               &Surface::surfaceMaterialSharedPtr)
        .def_static("createCylinder",
                    [](const Transform3& transform,
                       const std::shared_ptr<const CylinderBounds>& bounds) {
                      return Surface::makeShared<CylinderSurface>(transform,
                                                                  bounds);
                    })
        .def_static("createDisc",
                    [](const Transform3& transform,
                       const std::shared_ptr<const DiscBounds>& bounds) {
                      return Surface::makeShared<DiscSurface>(transform,
                                                              bounds);
                    })
        .def_static("createDisc",
                    [](const std::shared_ptr<const DiscBounds>& bounds,
                       const SurfacePlacementBase& detelement) {
                      return Surface::makeShared<DiscSurface>(bounds,
                                                              detelement);
                    })
        .def_static("createStraw",
                    [](const Transform3& transform,
                       const std::shared_ptr<const LineBounds>& bounds) {
                      return Surface::makeShared<StrawSurface>(transform,
                                                               bounds);
                    })
        .def_static("createStraw",
                    [](const std::shared_ptr<const LineBounds>& bounds,
                       const SurfacePlacementBase& detelement) {
                      return Surface::makeShared<StrawSurface>(bounds,
                                                               detelement);
                    })
        .def_static("createPerigee",
                    [](const Vector3& vertex) {
                      return Surface::makeShared<PerigeeSurface>(vertex);
                    })
        .def_static("createPlane",
                    [](const Transform3& transform,
                       const std::shared_ptr<const PlanarBounds>& bounds) {
                      return Surface::makeShared<PlaneSurface>(transform,
                                                               bounds);
                    })
        .def_static("createPlane",
                    [](const std::shared_ptr<const PlanarBounds>& pbounds,
                       const SurfacePlacementBase& detelement) {
                      return Surface::makeShared<PlaneSurface>(pbounds,
                                                               detelement);
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

    py::class_<LineSurface, Surface, std::shared_ptr<LineSurface>>(
        m, "LineSurface");
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
