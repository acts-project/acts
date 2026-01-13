// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeStack.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryVisitor.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Geometry/VolumeResizeStrategy.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/RangeXD.hpp"
#include "Acts/Visualization/ViewConfig.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <array>
#include <memory>
#include <numbers>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace {
struct GeometryIdentifierHookBinding : public GeometryIdentifierHook {
  py::object callable;

  GeometryIdentifier decorateIdentifier(GeometryIdentifier identifier,
                                        const Surface& surface) const override {
    return callable(identifier, surface.getSharedPtr())
        .cast<GeometryIdentifier>();
  }
};

struct MaterialSurfaceSelector {
  std::vector<const Surface*> surfaces = {};

  /// @param surface is the test surface
  void operator()(const Surface* surface) {
    if (surface->surfaceMaterial() != nullptr &&
        !rangeContainsValue(surfaces, surface)) {
      surfaces.push_back(surface);
    }
  }
};

#define _INVOKE(method, name, arg)                                  \
  pybind11::gil_scoped_acquire gil;                                 \
  pybind11::function override = pybind11::get_override(this, name); \
  if (!override) {                                                  \
    method(arg);                                                    \
    return;                                                         \
  }                                                                 \
  override(&arg);

// We only implement the mutable visitor here, because pybind11 always casts
// away const ness in any case
class PyTrackingGeometryVisitor : public TrackingGeometryMutableVisitor {
 public:
  void visitVolume(TrackingVolume& volume) override {
    _INVOKE(TrackingGeometryMutableVisitor::visitVolume, "visitVolume", volume);
  }

  void visitPortal(Portal& portal) override {
    _INVOKE(TrackingGeometryMutableVisitor::visitPortal, "visitPortal", portal);
  }

  void visitLayer(Layer& layer) override {
    _INVOKE(TrackingGeometryMutableVisitor::visitLayer, "visitLayer", layer);
  }

  void visitSurface(Surface& surface) override {
    _INVOKE(TrackingGeometryMutableVisitor::visitSurface, "visitSurface",
            surface);
  }

  void visitBoundarySurface(
      Acts::BoundarySurfaceT<Acts::TrackingVolume>& boundary) override {
    _INVOKE(Acts::TrackingGeometryMutableVisitor::visitBoundarySurface,
            "visitBoundarySurface", boundary);
  }
};

#undef _INVOKE

}  // namespace

namespace ActsPython {

/// This adds the geometry bindings to the python module
/// @param m the module to add the bindings to
void addGeometry(py::module_& m) {
  {
    py::class_<GeometryContext>(m, "GeometryContext")
        .def(py::init([]() {
          // Issue Python warning about deprecated default constructor
          auto warnings = py::module_::import("warnings");
          auto builtins = py::module_::import("builtins");
          warnings.attr("warn")(
              "GeometryContext::dangerouslyDefaultConstruct() is deprecated. Use "
              "GeometryContext.dangerouslyDefaultConstruct() instead to "
              "make empty context construction explicit.",
              builtins.attr("DeprecationWarning"));
          return GeometryContext::dangerouslyDefaultConstruct();
        }))  // Keep for backward compatibility but warn
        .def_static("dangerouslyDefaultConstruct",
                    &GeometryContext::dangerouslyDefaultConstruct,
                    "Create a default GeometryContext (empty, no alignment "
                    "data)");

    py::class_<GeometryIdentifier>(m, "GeometryIdentifier")
        .def(py::init<>())
        .def(py::init<GeometryIdentifier::Value>())
        .def(py::init([](int volume, int boundary, int layer, int approach,
                         int sensitive, int extra) {
               return GeometryIdentifier()
                   .withVolume(volume)
                   .withBoundary(boundary)
                   .withLayer(layer)
                   .withApproach(approach)
                   .withSensitive(sensitive)
                   .withExtra(extra);
             }),
             py::arg("volume") = 0, py::arg("boundary") = 0,
             py::arg("layer") = 0, py::arg("approach") = 0,
             py::arg("sensitive") = 0, py::arg("extra") = 0)
        .def_property(
            "layer", &GeometryIdentifier::layer,
            [](GeometryIdentifier& self, GeometryIdentifier::Value value) {
              self = self.withLayer(value);
            })
        .def_property(
            "volume", &GeometryIdentifier::volume,
            [](GeometryIdentifier& self, GeometryIdentifier::Value value) {
              self = self.withVolume(value);
            })
        .def_property(
            "boundary", &GeometryIdentifier::boundary,
            [](GeometryIdentifier& self, GeometryIdentifier::Value value) {
              self = self.withBoundary(value);
            })
        .def_property(
            "approach", &GeometryIdentifier::approach,
            [](GeometryIdentifier& self, GeometryIdentifier::Value value) {
              self = self.withApproach(value);
            })
        .def_property(
            "sensitive", &GeometryIdentifier::sensitive,
            [](GeometryIdentifier& self, GeometryIdentifier::Value value) {
              self = self.withSensitive(value);
            })
        .def_property(
            "extra", &GeometryIdentifier::extra,
            [](GeometryIdentifier& self, GeometryIdentifier::Value value) {
              self = self.withExtra(value);
            })
        .def_property_readonly("value", &GeometryIdentifier::value)
        .def("__str__", [](const GeometryIdentifier& self) {
          std::stringstream ss;
          ss << self;
          return ss.str();
        });
  }

  {
    py::class_<DetectorElementBase, std::shared_ptr<DetectorElementBase>>(
        m, "DetectorElementBase");
  }

  {
    py::enum_<VolumeBounds::BoundsType>(m, "VolumeBoundsType")
        .value("Cone", VolumeBounds::BoundsType::eCone)
        .value("Cuboid", VolumeBounds::BoundsType::eCuboid)
        .value("CutoutCylinder", VolumeBounds::BoundsType::eCutoutCylinder)
        .value("Cylinder", VolumeBounds::BoundsType::eCylinder)
        .value("GenericCuboid", VolumeBounds::BoundsType::eGenericCuboid)
        .value("Trapezoid", VolumeBounds::BoundsType::eTrapezoid)
        .value("Other", VolumeBounds::BoundsType::eOther);
  }

  {
    auto trkGeo =
        py::class_<TrackingGeometry, std::shared_ptr<TrackingGeometry>>(
            m, "TrackingGeometry")
            .def(py::init(
                [](const MutableTrackingVolumePtr& volPtr,
                   const std::shared_ptr<const IMaterialDecorator>& matDec,
                   const GeometryIdentifierHook& hook,
                   Acts::Logging::Level level) {
                  auto logger =
                      Acts::getDefaultLogger("TrackingGeometry", level);
                  auto obj = std::make_shared<Acts::TrackingGeometry>(
                      volPtr, matDec.get(), hook, *logger);
                  return obj;
                }))
            .def("visitSurfaces",
                 [](TrackingGeometry& self, py::function& func) {
                   self.visitSurfaces(func);
                 })
            .def("geoIdSurfaceMap", &TrackingGeometry::geoIdSurfaceMap)
            .def("extractMaterialSurfaces",
                 [](TrackingGeometry& self) {
                   MaterialSurfaceSelector selector;
                   self.visitSurfaces(selector, false);
                   return selector.surfaces;
                 })
            .def_property_readonly("highestTrackingVolume",
                                   &TrackingGeometry::highestTrackingVolumePtr)
            .def("visualize", &TrackingGeometry::visualize, py::arg("helper"),
                 py::arg("gctx"), py::arg("viewConfig") = s_viewVolume,
                 py::arg("portalViewConfig") = s_viewPortal,
                 py::arg("sensitiveViewConfig") = s_viewSensitive);

    using apply_ptr_t =
        void (TrackingGeometry::*)(TrackingGeometryMutableVisitor&);

    trkGeo.def("apply", static_cast<apply_ptr_t>(&TrackingGeometry::apply));
  }

  {
    py::class_<VolumeBounds, std::shared_ptr<VolumeBounds>>(m, "VolumeBounds")
        .def("type", &VolumeBounds::type)
        .def("__str__", [](const VolumeBounds& self) {
          std::stringstream ss;
          ss << self;
          return ss.str();
        });

    auto cvb =
        py::class_<CylinderVolumeBounds, std::shared_ptr<CylinderVolumeBounds>,
                   VolumeBounds>(m, "CylinderVolumeBounds")
            .def(py::init<double, double, double, double, double, double,
                          double>(),
                 "rmin"_a, "rmax"_a, "halfz"_a, "halfphi"_a = std::numbers::pi,
                 "avgphi"_a = 0., "bevelMinZ"_a = 0., "bevelMaxZ"_a = 0.);

    py::enum_<CylinderVolumeBounds::Face>(cvb, "Face")
        .value("PositiveDisc", CylinderVolumeBounds::Face::PositiveDisc)
        .value("NegativeDisc", CylinderVolumeBounds::Face::NegativeDisc)
        .value("OuterCylinder", CylinderVolumeBounds::Face::OuterCylinder)
        .value("InnerCylinder", CylinderVolumeBounds::Face::InnerCylinder)
        .value("NegativePhiPlane", CylinderVolumeBounds::Face::NegativePhiPlane)
        .value("PositivePhiPlane",
               CylinderVolumeBounds::Face::PositivePhiPlane);
  }

  {
    py::class_<Volume, std::shared_ptr<Volume>>(m, "Volume");

    py::class_<TrackingVolume, Volume, std::shared_ptr<TrackingVolume>>(
        m, "TrackingVolume")
        .def(py::init<const Transform3&, std::shared_ptr<VolumeBounds>,
                      std::string>());
  }

  {
    py::class_<GeometryIdentifierHook, std::shared_ptr<GeometryIdentifierHook>>(
        m, "GeometryIdentifierHook")
        .def(py::init([](py::object callable) {
          auto hook = std::make_shared<GeometryIdentifierHookBinding>();
          hook->callable = std::move(callable);
          return hook;
        }));
  }

  {
    py::class_<TrackingGeometryMutableVisitor, PyTrackingGeometryVisitor,
               std::shared_ptr<TrackingGeometryMutableVisitor>>(
        m, "TrackingGeometryMutableVisitor")
        .def(py::init<>());
  }

  py::class_<ExtentEnvelope>(m, "ExtentEnvelope")
      .def(py::init<>())
      .def(py::init<const Envelope&>())
      .def(py::init([](Envelope x, Envelope y, Envelope z, Envelope r,
                       Envelope phi, Envelope rPhi, Envelope theta,
                       Envelope eta, Envelope mag) {
             return ExtentEnvelope({.x = x,
                                    .y = y,
                                    .z = z,
                                    .r = r,
                                    .phi = phi,
                                    .rPhi = rPhi,
                                    .theta = theta,
                                    .eta = eta,
                                    .mag = mag});
           }),
           py::arg("x") = zeroEnvelope, py::arg("y") = zeroEnvelope,
           py::arg("z") = zeroEnvelope, py::arg("r") = zeroEnvelope,
           py::arg("phi") = zeroEnvelope, py::arg("rPhi") = zeroEnvelope,
           py::arg("theta") = zeroEnvelope, py::arg("eta") = zeroEnvelope,
           py::arg("mag") = zeroEnvelope)
      .def_static("Zero", &ExtentEnvelope::Zero)
      .def("__getitem__", [](ExtentEnvelope& self,
                             AxisDirection bValue) { return self[bValue]; })
      .def("__setitem__", [](ExtentEnvelope& self, AxisDirection bValue,
                             const Envelope& value) { self[bValue] = value; })
      .def("__str__", [](const ExtentEnvelope& self) {
        std::array<std::string, numAxisDirections()> values;

        std::stringstream ss;
        for (AxisDirection val : allAxisDirections()) {
          ss << val << "=(" << self[val][0] << ", " << self[val][1] << ")";
          values.at(toUnderlying(val)) = ss.str();
          ss.str("");
        }

        ss.str("");
        ss << "ExtentEnvelope(";
        ss << boost::algorithm::join(values, ", ");
        ss << ")";
        return ss.str();
      });

  py::class_<Extent>(m, "Extent")
      .def(py::init<const ExtentEnvelope&>(),
           py::arg("envelope") = ExtentEnvelope::Zero())
      .def("range",
           [](const Extent& self, AxisDirection bval) -> std::array<double, 2> {
             return {self.min(bval), self.max(bval)};
           })
      .def("setRange",
           [](Extent& self, AxisDirection bval,
              const std::array<double, 2>& range) {
             self.set(bval, range[0], range[1]);
           })
      .def("__str__", &Extent::toString);

  {
    py::enum_<VolumeAttachmentStrategy>(m, "VolumeAttachmentStrategy")
        .value("Gap", VolumeAttachmentStrategy::Gap)
        .value("First", VolumeAttachmentStrategy::First)
        .value("Second", VolumeAttachmentStrategy::Second)
        .value("Midpoint", VolumeAttachmentStrategy::Midpoint);

    py::enum_<VolumeResizeStrategy>(m, "VolumeResizeStrategy")
        .value("Gap", VolumeResizeStrategy::Gap)
        .value("Expand", VolumeResizeStrategy::Expand);
  }

  py::class_<PortalShellBase>(m, "PortalShellBase");

  py::class_<ProtoLayer>(m, "ProtoLayer")
      .def(py::init<const GeometryContext&,
                    const std::vector<std::shared_ptr<Surface>>&,
                    const Transform3&>(),
           "gctx"_a, "surfaces"_a, "transform"_a = Transform3::Identity())
      .def("min", &ProtoLayer::min, "bval"_a, "addenv"_a = true)
      .def("max", &ProtoLayer::max, "bval"_a, "addenv"_a = true)
      .def_property_readonly("surfaces", &ProtoLayer::surfaces);
}

}  // namespace ActsPython
