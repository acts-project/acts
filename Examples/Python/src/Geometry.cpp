// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/CuboidalContainerBuilder.hpp"
#include "Acts/Detector/CylindricalContainerBuilder.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/DetectorVolumeBuilder.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/IndexedRootVolumeFinderBuilder.hpp"
#include "Acts/Detector/KdtSurfacesProvider.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/VolumeStructureBuilder.hpp"
#include "Acts/Detector/interface/IDetectorBuilder.hpp"
#include "Acts/Detector/interface/IDetectorComponentBuilder.hpp"
#include "Acts/Detector/interface/IExternalStructureBuilder.hpp"
#include "Acts/Detector/interface/IGeometryIdGenerator.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"
#include "Acts/Detector/interface/IRootVolumeFinderBuilder.hpp"
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
#include "ActsExamples/Geometry/VolumeAssociationTest.hpp"
#include "ActsPython/Utilities/Context.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <array>
#include <memory>
#include <numbers>
#include <unordered_map>
#include <vector>

#include <boost/algorithm/string/join.hpp>
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

struct MaterialSurfaceSelector {
  std::vector<const Acts::Surface*> surfaces = {};

  /// @param surface is the test surface
  void operator()(const Acts::Surface* surface) {
    if (surface->surfaceMaterial() != nullptr &&
        !rangeContainsValue(surfaces, surface)) {
      surfaces.push_back(surface);
    }
  }
};

struct IdentifierSurfacesCollector {
  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> surfaces;
  /// @param surface is the test surface
  void operator()(const Acts::Surface* surface) {
    surfaces[surface->geometryId()] = surface;
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
class PyTrackingGeometryVisitor : public Acts::TrackingGeometryMutableVisitor {
 public:
  void visitVolume(Acts::TrackingVolume& volume) override {
    _INVOKE(Acts::TrackingGeometryMutableVisitor::visitVolume, "visitVolume",
            volume);
  }

  void visitPortal(Acts::Portal& portal) override {
    _INVOKE(Acts::TrackingGeometryMutableVisitor::visitPortal, "visitPortal",
            portal);
  }

  void visitLayer(Acts::Layer& layer) override {
    _INVOKE(Acts::TrackingGeometryMutableVisitor::visitLayer, "visitLayer",
            layer);
  }

  void visitSurface(Acts::Surface& surface) override {
    _INVOKE(Acts::TrackingGeometryMutableVisitor::visitSurface, "visitSurface",
            surface);
  }
};

#undef _INVOKE

}  // namespace

namespace ActsPython {

void addBlueprint(Context& ctx);

void addGeometry(Context& ctx) {
  auto m = ctx.get("main");
  {
    py::class_<Acts::GeometryIdentifier>(m, "GeometryIdentifier")
        .def(py::init<>())
        .def(py::init<Acts::GeometryIdentifier::Value>())
        .def_property("layer", &Acts::GeometryIdentifier::layer,
                      [](Acts::GeometryIdentifier& self,
                         Acts::GeometryIdentifier::Value value) {
                        self = self.withLayer(value);
                      })
        .def_property("volume", &Acts::GeometryIdentifier::volume,
                      [](Acts::GeometryIdentifier& self,
                         Acts::GeometryIdentifier::Value value) {
                        self = self.withVolume(value);
                      })
        .def_property("boundary", &Acts::GeometryIdentifier::boundary,
                      [](Acts::GeometryIdentifier& self,
                         Acts::GeometryIdentifier::Value value) {
                        self = self.withBoundary(value);
                      })
        .def_property("approach", &Acts::GeometryIdentifier::approach,
                      [](Acts::GeometryIdentifier& self,
                         Acts::GeometryIdentifier::Value value) {
                        self = self.withApproach(value);
                      })
        .def_property("sensitive", &Acts::GeometryIdentifier::sensitive,
                      [](Acts::GeometryIdentifier& self,
                         Acts::GeometryIdentifier::Value value) {
                        self = self.withSensitive(value);
                      })
        .def_property("extra", &Acts::GeometryIdentifier::extra,
                      [](Acts::GeometryIdentifier& self,
                         Acts::GeometryIdentifier::Value value) {
                        self = self.withExtra(value);
                      })
        .def_property_readonly("value", &Acts::GeometryIdentifier::value)
        .def("__str__", [](const Acts::GeometryIdentifier& self) {
          std::stringstream ss;
          ss << self;
          return ss.str();
        });
  }

  {
    py::class_<Acts::Surface, std::shared_ptr<Acts::Surface>>(m, "Surface")
        // Can't bind directly because GeometryObject is virtual base of Surface
        .def_property_readonly(
            "geometryId",
            [](const Acts::Surface& self) { return self.geometryId(); })
        .def("center", &Acts::Surface::center)
        .def_property_readonly("type", &Acts::Surface::type)
        .def("visualize", &Acts::Surface::visualize)
        .def_property_readonly("surfaceMaterial",
                               &Acts::Surface::surfaceMaterialSharedPtr);
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
    py::enum_<Acts::VolumeBounds::BoundsType>(m, "VolumeBoundsType")
        .value("Cone", Acts::VolumeBounds::BoundsType::eCone)
        .value("Cuboid", Acts::VolumeBounds::BoundsType::eCuboid)
        .value("CutoutCylinder",
               Acts::VolumeBounds::BoundsType::eCutoutCylinder)
        .value("Cylinder", Acts::VolumeBounds::BoundsType::eCylinder)
        .value("GenericCuboid", Acts::VolumeBounds::BoundsType::eGenericCuboid)
        .value("Trapezoid", Acts::VolumeBounds::BoundsType::eTrapezoid)
        .value("Other", Acts::VolumeBounds::BoundsType::eOther);
  }

  {
    auto trkGeo =
        py::class_<Acts::TrackingGeometry,
                   std::shared_ptr<Acts::TrackingGeometry>>(m,
                                                            "TrackingGeometry")
            .def(py::init(
                [](const Acts::MutableTrackingVolumePtr& volPtr,
                   std::shared_ptr<const Acts::IMaterialDecorator> matDec,
                   const Acts::GeometryIdentifierHook& hook,
                   Acts::Logging::Level level) {
                  auto logger =
                      Acts::getDefaultLogger("TrackingGeometry", level);
                  auto obj = std::make_shared<Acts::TrackingGeometry>(
                      volPtr, matDec.get(), hook, *logger);
                  return obj;
                }))
            .def("visitSurfaces",
                 [](Acts::TrackingGeometry& self, py::function& func) {
                   self.visitSurfaces(func);
                 })
            .def("geoIdSurfaceMap", &Acts::TrackingGeometry::geoIdSurfaceMap)
            .def("extractMaterialSurfaces",
                 [](Acts::TrackingGeometry& self) {
                   MaterialSurfaceSelector selector;
                   self.visitSurfaces(selector, false);
                   return selector.surfaces;
                 })
            .def_property_readonly(
                "highestTrackingVolume",
                &Acts::TrackingGeometry::highestTrackingVolumePtr)
            .def("visualize", &Acts::TrackingGeometry::visualize,
                 py::arg("helper"), py::arg("gctx"),
                 py::arg("viewConfig") = Acts::s_viewVolume,
                 py::arg("portalViewConfig") = Acts::s_viewPortal,
                 py::arg("sensitiveViewConfig") = Acts::s_viewSensitive);

    using apply_ptr_t =
        void (Acts::TrackingGeometry::*)(Acts::TrackingGeometryMutableVisitor&);

    trkGeo.def("apply",
               static_cast<apply_ptr_t>(&Acts::TrackingGeometry::apply));
  }

  {
    py::class_<Acts::VolumeBounds, std::shared_ptr<Acts::VolumeBounds>>(
        m, "VolumeBounds")
        .def("type", &Acts::VolumeBounds::type)
        .def("__str__", [](const Acts::VolumeBounds& self) {
          std::stringstream ss;
          ss << self;
          return ss.str();
        });

    auto cvb =
        py::class_<Acts::CylinderVolumeBounds,
                   std::shared_ptr<Acts::CylinderVolumeBounds>,
                   Acts::VolumeBounds>(m, "CylinderVolumeBounds")
            .def(py::init<double, double, double, double, double, double,
                          double>(),
                 "rmin"_a, "rmax"_a, "halfz"_a, "halfphi"_a = std::numbers::pi,
                 "avgphi"_a = 0., "bevelMinZ"_a = 0., "bevelMaxZ"_a = 0.);

    py::enum_<Acts::CylinderVolumeBounds::Face>(cvb, "Face")
        .value("PositiveDisc", Acts::CylinderVolumeBounds::Face::PositiveDisc)
        .value("NegativeDisc", Acts::CylinderVolumeBounds::Face::NegativeDisc)
        .value("OuterCylinder", Acts::CylinderVolumeBounds::Face::OuterCylinder)
        .value("InnerCylinder", Acts::CylinderVolumeBounds::Face::InnerCylinder)
        .value("NegativePhiPlane",
               Acts::CylinderVolumeBounds::Face::NegativePhiPlane)
        .value("PositivePhiPlane",
               Acts::CylinderVolumeBounds::Face::PositivePhiPlane);
  }

  {
    py::class_<Acts::Volume, std::shared_ptr<Acts::Volume>>(m, "Volume");

    py::class_<Acts::TrackingVolume, Acts::Volume,
               std::shared_ptr<Acts::TrackingVolume>>(m, "TrackingVolume")
        .def(py::init<const Acts::Transform3&,
                      std::shared_ptr<Acts::VolumeBounds>, std::string>());
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
    py::class_<Acts::TrackingGeometryMutableVisitor, PyTrackingGeometryVisitor,
               std::shared_ptr<Acts::TrackingGeometryMutableVisitor>>(
        m, "TrackingGeometryMutableVisitor")
        .def(py::init<>());
  }

  py::class_<Acts::ExtentEnvelope>(m, "ExtentEnvelope")
      .def(py::init<>())
      .def(py::init<const Acts::Envelope&>())
      .def(py::init([](Acts::Envelope x, Acts::Envelope y, Acts::Envelope z,
                       Acts::Envelope r, Acts::Envelope phi,
                       Acts::Envelope rPhi, Acts::Envelope theta,
                       Acts::Envelope eta, Acts::Envelope mag) {
             return Acts::ExtentEnvelope({.x = x,
                                          .y = y,
                                          .z = z,
                                          .r = r,
                                          .phi = phi,
                                          .rPhi = rPhi,
                                          .theta = theta,
                                          .eta = eta,
                                          .mag = mag});
           }),
           py::arg("x") = Acts::zeroEnvelope, py::arg("y") = Acts::zeroEnvelope,
           py::arg("z") = Acts::zeroEnvelope, py::arg("r") = Acts::zeroEnvelope,
           py::arg("phi") = Acts::zeroEnvelope,
           py::arg("rPhi") = Acts::zeroEnvelope,
           py::arg("theta") = Acts::zeroEnvelope,
           py::arg("eta") = Acts::zeroEnvelope,
           py::arg("mag") = Acts::zeroEnvelope)
      .def_static("Zero", &Acts::ExtentEnvelope::Zero)
      .def("__getitem__",
           [](Acts::ExtentEnvelope& self, Acts::AxisDirection bValue) {
             return self[bValue];
           })
      .def("__setitem__",
           [](Acts::ExtentEnvelope& self, Acts::AxisDirection bValue,
              const Acts::Envelope& value) { self[bValue] = value; })
      .def("__str__", [](const Acts::ExtentEnvelope& self) {
        std::array<std::string, Acts::numAxisDirections()> values;

        std::stringstream ss;
        for (Acts::AxisDirection val : Acts::allAxisDirections()) {
          ss << val << "=(" << self[val][0] << ", " << self[val][1] << ")";
          values.at(Acts::toUnderlying(val)) = ss.str();
          ss.str("");
        }

        ss.str("");
        ss << "ExtentEnvelope(";
        ss << boost::algorithm::join(values, ", ");
        ss << ")";
        return ss.str();
      });

  py::class_<Acts::Extent>(m, "Extent")
      .def(py::init<const Acts::ExtentEnvelope&>(),
           py::arg("envelope") = Acts::ExtentEnvelope::Zero())
      .def("range",
           [](const Acts::Extent& self,
              Acts::AxisDirection bval) -> std::array<double, 2> {
             return {self.min(bval), self.max(bval)};
           })
      .def("setRange",
           [](Acts::Extent& self, Acts::AxisDirection bval,
              const std::array<double, 2>& range) {
             self.set(bval, range[0], range[1]);
           })
      .def("__str__", &Acts::Extent::toString);

  {
    py::enum_<Acts::VolumeAttachmentStrategy>(m, "VolumeAttachmentStrategy")
        .value("Gap", Acts::VolumeAttachmentStrategy::Gap)
        .value("First", Acts::VolumeAttachmentStrategy::First)
        .value("Second", Acts::VolumeAttachmentStrategy::Second)
        .value("Midpoint", Acts::VolumeAttachmentStrategy::Midpoint);

    py::enum_<Acts::VolumeResizeStrategy>(m, "VolumeResizeStrategy")
        .value("Gap", Acts::VolumeResizeStrategy::Gap)
        .value("Expand", Acts::VolumeResizeStrategy::Expand);
  }

  py::class_<Acts::PortalShellBase>(m, "PortalShellBase");

  addBlueprint(ctx);
}

void addExperimentalGeometry(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  using namespace Acts::Experimental;

  // Detector volume definition
  py::class_<DetectorVolume, std::shared_ptr<DetectorVolume>>(m,
                                                              "DetectorVolume")
      .def("surfaces", &DetectorVolume::surfaces)
      .def("surfacePtrs", &DetectorVolume::surfacePtrs);

  // Detector definition
  py::class_<Detector, std::shared_ptr<Detector>>(m, "Detector")
      .def("volumes", &Detector::volumes)
      .def("volumePtrs", &Detector::volumePtrs)
      .def("numberVolumes",
           [](const Detector& self) { return self.volumes().size(); })
      .def("extractMaterialSurfaces",
           [](const Detector& self) {
             MaterialSurfaceSelector selector;
             self.visitSurfaces(selector);
             return selector.surfaces;
           })
      .def("geoIdSurfaceMap",
           [](const Detector& self) {
             IdentifierSurfacesCollector collector;
             self.visitSurfaces(collector);
             return collector.surfaces;
           })
      .def("cylindricalVolumeRepresentation",
           [](const Detector& self, const Acts::GeometryContext& gctx) {
             // Loop over the volumes and gather the extent
             Acts::Extent extent;
             for (const auto& volume : self.volumes()) {
               extent.extend(volume->extent(gctx));
             }
             auto bounds = std::make_shared<Acts::CylinderVolumeBounds>(
                 0., extent.max(Acts::AxisDirection::AxisR),
                 extent.max(Acts::AxisDirection::AxisZ));

             return std::make_shared<Acts::Volume>(Acts::Transform3::Identity(),
                                                   std::move(bounds));
           });

  // Portal definition
  py::class_<Acts::Experimental::Portal,
             std::shared_ptr<Acts::Experimental::Portal>>(m, "Portal");

  {
    // The surface hierarchy map
    using SurfaceHierarchyMap =
        Acts::GeometryHierarchyMap<std::shared_ptr<Acts::Surface>>;

    py::class_<SurfaceHierarchyMap, std::shared_ptr<SurfaceHierarchyMap>>(
        m, "SurfaceHierarchyMap");

    // Extract volume / layer surfaces
    mex.def("extractVolumeLayerSurfaces", [](const SurfaceHierarchyMap& smap,
                                             bool sensitiveOnly) {
      std::map<
          unsigned int,
          std::map<unsigned int, std::vector<std::shared_ptr<Acts::Surface>>>>
          surfaceVolumeLayerMap;
      for (const auto& surface : smap) {
        auto gid = surface->geometryId();
        // Exclusion criteria
        if (sensitiveOnly && gid.sensitive() == 0) {
          continue;
        };
        surfaceVolumeLayerMap[gid.volume()][gid.layer()].push_back(surface);
      }
      // Return the surface volume map
      return surfaceVolumeLayerMap;
    });
  }

  {
    // Be able to construct a proto binning
    py::class_<Acts::ProtoAxis>(m, "ProtoAxis")
        .def(py::init<Acts::AxisBoundaryType, const std::vector<double>&>(),
             "bType"_a, "e"_a)
        .def(py::init<Acts::AxisBoundaryType, double, double, std::size_t>(),
             "bType"_a, "minE"_a, "maxE"_a, "nbins"_a)
        .def(py::init<Acts::AxisBoundaryType, std::size_t>(), "bType"_a,
             "nbins"_a);

    py::class_<Acts::DirectedProtoAxis>(m, "DirectedProtoAxis")
        .def(py::init<Acts::AxisDirection, Acts::AxisBoundaryType,
                      const std::vector<double>&>(),
             "bValue"_a, "bType"_a, "e"_a)
        .def(py::init<Acts::AxisDirection, Acts::AxisBoundaryType, double,
                      double, std::size_t>(),
             "bValue"_a, "bType"_a, "minE"_a, "maxE"_a, "nbins"_a)
        .def(py::init<Acts::AxisDirection, Acts::AxisBoundaryType,
                      std::size_t>(),
             "bValue"_a, "bType"_a, "nbins"_a);
  }

  {
    // The internal layer structure builder
    py::class_<Acts::Experimental::IInternalStructureBuilder,
               std::shared_ptr<Acts::Experimental::IInternalStructureBuilder>>(
        m, "IInternalStructureBuilder");

    auto lsBuilder =
        py::class_<Acts::Experimental::LayerStructureBuilder,
                   Acts::Experimental::IInternalStructureBuilder,
                   std::shared_ptr<Acts::Experimental::LayerStructureBuilder>>(
            m, "LayerStructureBuilder")
            .def(py::init(
                [](const Acts::Experimental::LayerStructureBuilder::Config&
                       config,
                   const std::string& name, Acts::Logging::Level level) {
                  return std::make_shared<
                      Acts::Experimental::LayerStructureBuilder>(
                      config, Acts::getDefaultLogger(name, level));
                }));

    auto lsConfig =
        py::class_<Acts::Experimental::LayerStructureBuilder::Config>(lsBuilder,
                                                                      "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(lsConfig, surfacesProvider, supports, binnings,
                       quarterSegments, auxiliary);

    // The internal layer structure builder
    py::class_<Acts::Experimental::ISurfacesProvider,
               std::shared_ptr<Acts::Experimental::ISurfacesProvider>>(
        m, "ISurfacesProvider");

    py::class_<Acts::Experimental::LayerStructureBuilder::SurfacesHolder,
               Acts::Experimental::ISurfacesProvider,
               std::shared_ptr<
                   Acts::Experimental::LayerStructureBuilder::SurfacesHolder>>(
        lsBuilder, "SurfacesHolder")
        .def(py::init<std::vector<std::shared_ptr<Acts::Surface>>>());
  }

  {
    using RangeXDDim1 = Acts::RangeXD<1u, double>;
    using KdtSurfacesDim1Bin100 = Acts::Experimental::KdtSurfaces<1u, 100u>;
    using KdtSurfacesProviderDim1Bin100 =
        Acts::Experimental::KdtSurfacesProvider<1u, 100u>;

    py::class_<RangeXDDim1>(m, "RangeXDDim1")
        .def(py::init([](const std::array<double, 2u>& irange) {
          RangeXDDim1 range;
          range[0].shrink(irange[0], irange[1]);
          return range;
        }));

    py::class_<KdtSurfacesDim1Bin100, std::shared_ptr<KdtSurfacesDim1Bin100>>(
        m, "KdtSurfacesDim1Bin100")
        .def(py::init<const Acts::GeometryContext&,
                      const std::vector<std::shared_ptr<Acts::Surface>>&,
                      const std::array<Acts::AxisDirection, 1u>&>())
        .def("surfaces", py::overload_cast<const RangeXDDim1&>(
                             &KdtSurfacesDim1Bin100::surfaces, py::const_));

    py::class_<KdtSurfacesProviderDim1Bin100,
               Acts::Experimental::ISurfacesProvider,
               std::shared_ptr<KdtSurfacesProviderDim1Bin100>>(
        m, "KdtSurfacesProviderDim1Bin100")
        .def(py::init<std::shared_ptr<KdtSurfacesDim1Bin100>,
                      const Acts::Extent&>());
  }

  {
    using RangeXDDim2 = Acts::RangeXD<2u, double>;
    using KdtSurfacesDim2Bin100 = Acts::Experimental::KdtSurfaces<2u, 100u>;
    using KdtSurfacesProviderDim2Bin100 =
        Acts::Experimental::KdtSurfacesProvider<2u, 100u>;

    py::class_<RangeXDDim2>(m, "RangeXDDim2")
        .def(py::init([](const std::array<double, 2u>& range0,
                         const std::array<double, 2u>& range1) {
          RangeXDDim2 range;
          range[0].shrink(range0[0], range0[1]);
          range[1].shrink(range1[0], range1[1]);
          return range;
        }));

    py::class_<KdtSurfacesDim2Bin100, std::shared_ptr<KdtSurfacesDim2Bin100>>(
        m, "KdtSurfacesDim2Bin100")
        .def(py::init<const Acts::GeometryContext&,
                      const std::vector<std::shared_ptr<Acts::Surface>>&,
                      const std::array<Acts::AxisDirection, 2u>&>())
        .def("surfaces", py::overload_cast<const RangeXDDim2&>(
                             &KdtSurfacesDim2Bin100::surfaces, py::const_));

    py::class_<KdtSurfacesProviderDim2Bin100,
               Acts::Experimental::ISurfacesProvider,
               std::shared_ptr<KdtSurfacesProviderDim2Bin100>>(
        m, "KdtSurfacesProviderDim2Bin100")
        .def(py::init<std::shared_ptr<KdtSurfacesDim2Bin100>,
                      const Acts::Extent&>());
  }

  {
    using RangeXDDim3 = Acts::RangeXD<3u, double>;

    py::class_<RangeXDDim3>(m, "RangeXDDim3")
        .def(py::init([](const std::array<double, 2u>& range0,
                         const std::array<double, 2u>& range1,
                         const std::array<double, 2u>& range2) {
          RangeXDDim3 range;
          range[0].shrink(range0[0], range0[1]);
          range[1].shrink(range1[0], range1[1]);
          range[2].shrink(range2[0], range2[1]);
          return range;
        }));
  }

  {
    // The external volume structure builder
    py::class_<Acts::Experimental::IExternalStructureBuilder,
               std::shared_ptr<Acts::Experimental::IExternalStructureBuilder>>(
        m, "IExternalStructureBuilder");

    auto vsBuilder =
        py::class_<Acts::Experimental::VolumeStructureBuilder,
                   Acts::Experimental::IExternalStructureBuilder,
                   std::shared_ptr<Acts::Experimental::VolumeStructureBuilder>>(
            m, "VolumeStructureBuilder")
            .def(py::init(
                [](const Acts::Experimental::VolumeStructureBuilder::Config&
                       config,
                   const std::string& name, Acts::Logging::Level level) {
                  return std::make_shared<VolumeStructureBuilder>(
                      config, Acts::getDefaultLogger(name, level));
                }));

    auto vsConfig =
        py::class_<Acts::Experimental::VolumeStructureBuilder::Config>(
            vsBuilder, "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(vsConfig, boundsType, boundValues, transform, auxiliary);
  }

  {
    py::class_<Acts::Experimental::IGeometryIdGenerator,
               std::shared_ptr<Acts::Experimental::IGeometryIdGenerator>>(
        m, "IGeometryIdGenerator");

    auto geoIdGen =
        py::class_<Acts::Experimental::GeometryIdGenerator,
                   Acts::Experimental::IGeometryIdGenerator,
                   std::shared_ptr<Acts::Experimental::GeometryIdGenerator>>(
            m, "GeometryIdGenerator")
            .def(py::init([](Acts::Experimental::GeometryIdGenerator::Config&
                                 config,
                             const std::string& name,
                             Acts::Logging::Level level) {
              return std::make_shared<Acts::Experimental::GeometryIdGenerator>(
                  config, Acts::getDefaultLogger(name, level));
            }));

    auto geoIdGenConfig =
        py::class_<Acts::Experimental::GeometryIdGenerator::Config>(geoIdGen,
                                                                    "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(geoIdGenConfig, containerMode, containerId,
                       resetSubCounters, overrideExistingIds);
  }

  {
    // Put them together to a detector volume
    py::class_<Acts::Experimental::IDetectorComponentBuilder,
               std::shared_ptr<Acts::Experimental::IDetectorComponentBuilder>>(
        m, "IDetectorComponentBuilder");

    auto dvBuilder =
        py::class_<Acts::Experimental::DetectorVolumeBuilder,
                   Acts::Experimental::IDetectorComponentBuilder,
                   std::shared_ptr<Acts::Experimental::DetectorVolumeBuilder>>(
            m, "DetectorVolumeBuilder")
            .def(py::init(
                [](const Acts::Experimental::DetectorVolumeBuilder::Config&
                       config,
                   const std::string& name, Acts::Logging::Level level) {
                  return std::make_shared<
                      Acts::Experimental::DetectorVolumeBuilder>(
                      config, Acts::getDefaultLogger(name, level));
                }))
            .def("construct",
                 &Acts::Experimental::DetectorVolumeBuilder::construct);

    auto dvConfig =
        py::class_<Acts::Experimental::DetectorVolumeBuilder::Config>(dvBuilder,
                                                                      "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(dvConfig, name, internalsBuilder, externalsBuilder,
                       geoIdGenerator, auxiliary);
  }

  {
    // The external volume structure builder
    py::class_<Acts::Experimental::IRootVolumeFinderBuilder,
               std::shared_ptr<Acts::Experimental::IRootVolumeFinderBuilder>>(
        m, "IRootVolumeFinderBuilder");

    auto irvBuilder =
        py::class_<Acts::Experimental::IndexedRootVolumeFinderBuilder,
                   Acts::Experimental::IRootVolumeFinderBuilder,
                   std::shared_ptr<
                       Acts::Experimental::IndexedRootVolumeFinderBuilder>>(
            m, "IndexedRootVolumeFinderBuilder")
            .def(py::init<std::vector<Acts::AxisDirection>>());
  }

  {
    // Cylindrical container builder
    auto ccBuilder =
        py::class_<
            Acts::Experimental::CylindricalContainerBuilder,
            Acts::Experimental::IDetectorComponentBuilder,
            std::shared_ptr<Acts::Experimental::CylindricalContainerBuilder>>(
            m, "CylindricalContainerBuilder")
            .def(py::init([](const Acts::Experimental::
                                 CylindricalContainerBuilder::Config& config,
                             const std::string& name,
                             Acts::Logging::Level level) {
              return std::make_shared<CylindricalContainerBuilder>(
                  config, Acts::getDefaultLogger(name, level));
            }))
            .def("construct", &CylindricalContainerBuilder::construct);

    auto ccConfig =
        py::class_<CylindricalContainerBuilder::Config>(ccBuilder, "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(ccConfig, builders, binning, rootVolumeFinderBuilder,
                       geoIdGenerator, geoIdReverseGen, auxiliary);
  }

  {
    // Cuboidal container builder
    auto ccBuilder =
        py::class_<Acts::Experimental::CuboidalContainerBuilder,
                   Acts::Experimental::IDetectorComponentBuilder,
                   std::shared_ptr<CuboidalContainerBuilder>>(
            m, "CuboidalContainerBuilder")
            .def(py::init(
                [](const Acts::Experimental::CuboidalContainerBuilder::Config&
                       config,
                   const std::string& name, Acts::Logging::Level level) {
                  return std::make_shared<
                      Acts::Experimental::CuboidalContainerBuilder>(
                      config, Acts::getDefaultLogger(name, level));
                }))
            .def("construct",
                 &Acts::Experimental::CuboidalContainerBuilder::construct);

    auto ccConfig =
        py::class_<Acts::Experimental::CuboidalContainerBuilder::Config>(
            ccBuilder, "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(ccConfig, builders, binning, rootVolumeFinderBuilder,
                       geoIdGenerator, geoIdReverseGen, auxiliary);
  }

  {
    // Detector builder
    auto dBuilder =
        py::class_<Acts::Experimental::DetectorBuilder,
                   std::shared_ptr<DetectorBuilder>>(m, "DetectorBuilder")
            .def(py::init(
                [](const Acts::Experimental::DetectorBuilder::Config& config,
                   const std::string& name, Acts::Logging::Level level) {
                  return std::make_shared<Acts::Experimental::DetectorBuilder>(
                      config, Acts::getDefaultLogger(name, level));
                }))
            .def("construct", &Acts::Experimental::DetectorBuilder::construct);

    auto dConfig = py::class_<Acts::Experimental::DetectorBuilder::Config>(
                       dBuilder, "Config")
                       .def(py::init<>());
    ACTS_PYTHON_STRUCT(dConfig, name, builder, geoIdGenerator,
                       materialDecorator, auxiliary);
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::VolumeAssociationTest, mex,
                                "VolumeAssociationTest", name, ntests,
                                randomNumbers, randomRange, detector);

  py::class_<Acts::ProtoLayer>(m, "ProtoLayer")
      .def(py::init<const Acts::GeometryContext&,
                    const std::vector<std::shared_ptr<Acts::Surface>>&,
                    const Acts::Transform3&>(),
           "gctx"_a, "surfaces"_a, "transform"_a = Acts::Transform3::Identity())
      .def("min", &Acts::ProtoLayer::min, "bval"_a, "addenv"_a = true)
      .def("max", &Acts::ProtoLayer::max, "bval"_a, "addenv"_a = true)
      .def_property_readonly("surfaces", &Acts::ProtoLayer::surfaces);
}

}  // namespace ActsPython
