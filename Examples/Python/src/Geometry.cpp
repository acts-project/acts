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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

struct IdentifierSurfacesCollector {
  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> surfaces;
  /// @param surface is the test surface
  void operator()(const Acts::Surface* surface) {
    surfaces[surface->geometryId()] = surface;
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

}  // namespace

namespace ActsPython {

void addBlueprint(Context& ctx);

void addGeometryLegacy(Context& ctx) {

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

  // Add the surface hierarchy map
  using SurfaceHierarchyMap =
      Acts::GeometryHierarchyMap<std::shared_ptr<Acts::Surface>>;

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


}

}  // namespace ActsPython
