// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/CylindricalContainerBuilder.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/DetectorVolumeBuilder.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/VolumeStructureBuilder.hpp"
#include "Acts/Detector/interface/IDetectorBuilder.hpp"
#include "Acts/Detector/interface/IDetectorComponentBuilder.hpp"
#include "Acts/Detector/interface/IExternalStructureBuilder.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

#include <array>
#include <memory>
#include <vector>

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
}

void addExperimentalGeometry(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  using namespace Acts::Experimental;

  // Detector definition
  py::class_<Detector, std::shared_ptr<Detector>>(m, "Detector");

  // Detector volume definition
  py::class_<DetectorVolume, std::shared_ptr<DetectorVolume>>(m,
                                                              "DetectorVolume");

  {
    // The surface hiearchy map
    using SurfaceHierarchyMap =
        Acts::GeometryHierarchyMap<std::shared_ptr<Surface>>;

    py::class_<SurfaceHierarchyMap, std::shared_ptr<SurfaceHierarchyMap>>(
        m, "SurfaceHierarchyMap");

    // Extract volume / layer surfaces
    mex.def("extractVolumeLayerSurfaces", [](const SurfaceHierarchyMap& smap,
                                             bool sensitiveOnly) {
      std::map<unsigned int,
               std::map<unsigned int, std::vector<std::shared_ptr<Surface>>>>
          surfaceVolumeLayerMap;
      for (const auto& surface : smap) {
        auto gid = surface->geometryId();
        // Exclusion criteria
        if (sensitiveOnly and gid.sensitive() == 0) {
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
    py::class_<ProtoBinning>(m, "ProtoBinning")
        .def(py::init<Acts::BinningValue, Acts::detail::AxisBoundaryType,
                      const std::vector<Acts::ActsScalar>&, std::size_t>())
        .def(py::init<Acts::BinningValue, Acts::detail::AxisBoundaryType,
                      Acts::ActsScalar, Acts::ActsScalar, std::size_t,
                      std::size_t>());
  }

  {
    // The internal layer structure builder
    py::class_<Acts::Experimental::IInternalStructureBuilder,
               std::shared_ptr<Acts::Experimental::IInternalStructureBuilder>>(
        m, "IInternalStructureBuilder");

    auto lsBuilder =
        py::class_<LayerStructureBuilder,
                   Acts::Experimental::IInternalStructureBuilder,
                   std::shared_ptr<LayerStructureBuilder>>(
            m, "LayerStructureBuilder")
            .def(py::init([](const LayerStructureBuilder::Config& config,
                             const std::string& name,
                             Acts::Logging::Level level) {
              return std::make_shared<LayerStructureBuilder>(
                  config, getDefaultLogger(name, level));
            }));

    auto lsConfig =
        py::class_<LayerStructureBuilder::Config>(lsBuilder, "Config")
            .def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(lsConfig, LayerStructureBuilder::Config);
    ACTS_PYTHON_MEMBER(surfacesProvider);
    ACTS_PYTHON_MEMBER(supports);
    ACTS_PYTHON_MEMBER(binnings);
    ACTS_PYTHON_MEMBER(nSegments);
    ACTS_PYTHON_MEMBER(auxilliary);
    ACTS_PYTHON_STRUCT_END();

    // The internal layer structure builder
    py::class_<Acts::Experimental::ISurfacesProvider,
               std::shared_ptr<Acts::Experimental::ISurfacesProvider>>(
        m, "ISurfacesProvider");

    py::class_<LayerStructureBuilder::SurfacesHolder,
               Acts::Experimental::ISurfacesProvider,
               std::shared_ptr<LayerStructureBuilder::SurfacesHolder>>(
        lsBuilder, "SurfacesHolder")
        .def(py::init<std::vector<std::shared_ptr<Surface>>>());
  }

  {
    // The external volume structure builder
    py::class_<Acts::Experimental::IExternalStructureBuilder,
               std::shared_ptr<Acts::Experimental::IExternalStructureBuilder>>(
        m, "IExternalStructureBuilder");

    auto vsBuilder =
        py::class_<VolumeStructureBuilder,
                   Acts::Experimental::IExternalStructureBuilder,
                   std::shared_ptr<VolumeStructureBuilder>>(
            m, "VolumeStructureBuilder")
            .def(py::init([](const VolumeStructureBuilder::Config& config,
                             const std::string& name,
                             Acts::Logging::Level level) {
              return std::make_shared<VolumeStructureBuilder>(
                  config, getDefaultLogger(name, level));
            }));

    auto vsConfig =
        py::class_<VolumeStructureBuilder::Config>(vsBuilder, "Config")
            .def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(vsConfig, VolumeStructureBuilder::Config);
    ACTS_PYTHON_MEMBER(boundsType);
    ACTS_PYTHON_MEMBER(boundValues);
    ACTS_PYTHON_MEMBER(auxilliary);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    // Put them together to a detector volume
    py::class_<Acts::Experimental::IDetectorComponentBuilder,
               std::shared_ptr<Acts::Experimental::IDetectorComponentBuilder>>(
        m, "IDetectorComponentBuilder");

    auto dvBuilder =
        py::class_<DetectorVolumeBuilder,
                   Acts::Experimental::IDetectorComponentBuilder,
                   std::shared_ptr<DetectorVolumeBuilder>>(
            m, "DetectorVolumeBuilder")
            .def(py::init([](const DetectorVolumeBuilder::Config& config,
                             const std::string& name,
                             Acts::Logging::Level level) {
              return std::make_shared<DetectorVolumeBuilder>(
                  config, getDefaultLogger(name, level));
            }))
            .def("construct", &DetectorVolumeBuilder::construct);

    auto dvConfig =
        py::class_<DetectorVolumeBuilder::Config>(dvBuilder, "Config")
            .def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(dvConfig, DetectorVolumeBuilder::Config);
    ACTS_PYTHON_MEMBER(name);
    ACTS_PYTHON_MEMBER(internalsBuilder);
    ACTS_PYTHON_MEMBER(externalsBuilder);
    ACTS_PYTHON_MEMBER(auxilliary);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    // Cylindrical container builder
    auto ccBuilder =
        py::class_<CylindricalContainerBuilder,
                   Acts::Experimental::IDetectorComponentBuilder,
                   std::shared_ptr<CylindricalContainerBuilder>>(
            m, "CylindricalContainerBuilder")
            .def(py::init([](const CylindricalContainerBuilder::Config& config,
                             const std::string& name,
                             Acts::Logging::Level level) {
              return std::make_shared<CylindricalContainerBuilder>(
                  config, getDefaultLogger(name, level));
            }))
            .def("construct", &CylindricalContainerBuilder::construct);

    auto ccConfig =
        py::class_<CylindricalContainerBuilder::Config>(ccBuilder, "Config")
            .def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(ccConfig, CylindricalContainerBuilder::Config);
    ACTS_PYTHON_MEMBER(builders);
    ACTS_PYTHON_MEMBER(binning);
    ACTS_PYTHON_MEMBER(auxilliary);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    // Detector builder
    auto dBuilder =
        py::class_<DetectorBuilder, std::shared_ptr<DetectorBuilder>>(
            m, "DetectorBuilder")
            .def(py::init([](const DetectorBuilder::Config& config,
                             const std::string& name,
                             Acts::Logging::Level level) {
              return std::make_shared<DetectorBuilder>(
                  config, getDefaultLogger(name, level));
            }))
            .def("construct", &DetectorBuilder::construct);

    auto dConfig = py::class_<DetectorBuilder::Config>(dBuilder, "Config")
                       .def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(dConfig, DetectorBuilder::Config);
    ACTS_PYTHON_MEMBER(name);
    ACTS_PYTHON_MEMBER(builder);
    ACTS_PYTHON_MEMBER(auxilliary);
    ACTS_PYTHON_STRUCT_END();
  }
}

}  // namespace Acts::Python
