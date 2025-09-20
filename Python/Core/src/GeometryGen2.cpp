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
#include "Acts/Material/IMaterialDecorator.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
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

using namespace Acts;
using namespace Acts::Experimental;
using namespace ActsExamples;

namespace {

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

struct IdentifierSurfacesCollector {
  std::unordered_map<GeometryIdentifier, const Surface*> surfaces;
  /// @param surface is the test surface
  void operator()(const Surface* surface) {
    surfaces[surface->geometryId()] = surface;
  }
};

}  // namespace

namespace ActsPython {
/// This adds the geometry building bindings for the Gen2 geometry
/// @param m the module to add the bindings to
void addGeometryGen2(py::module_& m) {
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
           [](const Detector& self, const GeometryContext& gctx) {
             // Loop over the volumes and gather the extent
             Extent extent;
             for (const auto& volume : self.volumes()) {
               extent.extend(volume->extent(gctx));
             }
             auto bounds = std::make_shared<CylinderVolumeBounds>(
                 0., extent.max(AxisDirection::AxisR),
                 extent.max(AxisDirection::AxisZ));

             return std::make_shared<Volume>(Transform3::Identity(),
                                             std::move(bounds));
           });

  // Portal definition
  py::class_<Experimental::Portal, std::shared_ptr<Experimental::Portal>>(
      m, "Portal");

  {
    // The surface hierarchy map
    using SurfaceHierarchyMap = GeometryHierarchyMap<std::shared_ptr<Surface>>;

    // Extract volume / layer surfaces
    m.def("extractVolumeLayerSurfaces", [](const SurfaceHierarchyMap& smap,
                                           bool sensitiveOnly) {
      std::map<unsigned int,
               std::map<unsigned int, std::vector<std::shared_ptr<Surface>>>>
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
    py::class_<Experimental::IInternalStructureBuilder,
               std::shared_ptr<Experimental::IInternalStructureBuilder>>(
        m, "IInternalStructureBuilder");

    auto lsBuilder =
        py::class_<LayerStructureBuilder,
                   Experimental::IInternalStructureBuilder,
                   std::shared_ptr<LayerStructureBuilder>>(
            m, "LayerStructureBuilder")
            .def(py::init([](const LayerStructureBuilder::Config& config,
                             const std::string& name, Logging::Level level) {
              return std::make_shared<LayerStructureBuilder>(
                  config, getDefaultLogger(name, level));
            }));

    auto lsConfig =
        py::class_<LayerStructureBuilder::Config>(lsBuilder, "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(lsConfig, surfacesProvider, supports, binnings,
                       quarterSegments, auxiliary);

    // The internal layer structure builder
    py::class_<Experimental::ISurfacesProvider,
               std::shared_ptr<Experimental::ISurfacesProvider>>(
        m, "ISurfacesProvider");

    py::class_<LayerStructureBuilder::SurfacesHolder,
               Experimental::ISurfacesProvider,
               std::shared_ptr<LayerStructureBuilder::SurfacesHolder>>(
        lsBuilder, "SurfacesHolder")
        .def(py::init<std::vector<std::shared_ptr<Surface>>>());
  }

  {
    using RangeXDDim1 = RangeXD<1u, double>;
    using KdtSurfacesDim1Bin100 = Experimental::KdtSurfaces<1u, 100u>;
    using KdtSurfacesProviderDim1Bin100 =
        Experimental::KdtSurfacesProvider<1u, 100u>;

    py::class_<RangeXDDim1>(m, "RangeXDDim1")
        .def(py::init([](const std::array<double, 2u>& irange) {
          RangeXDDim1 range;
          range[0].shrink(irange[0], irange[1]);
          return range;
        }));

    py::class_<KdtSurfacesDim1Bin100, std::shared_ptr<KdtSurfacesDim1Bin100>>(
        m, "KdtSurfacesDim1Bin100")
        .def(py::init<const GeometryContext&,
                      const std::vector<std::shared_ptr<Surface>>&,
                      const std::array<AxisDirection, 1u>&>())
        .def("surfaces", py::overload_cast<const RangeXDDim1&>(
                             &KdtSurfacesDim1Bin100::surfaces, py::const_));

    py::class_<KdtSurfacesProviderDim1Bin100, Experimental::ISurfacesProvider,
               std::shared_ptr<KdtSurfacesProviderDim1Bin100>>(
        m, "KdtSurfacesProviderDim1Bin100")
        .def(py::init<std::shared_ptr<KdtSurfacesDim1Bin100>, const Extent&>());
  }

  {
    using RangeXDDim2 = RangeXD<2u, double>;
    using KdtSurfacesDim2Bin100 = Experimental::KdtSurfaces<2u, 100u>;
    using KdtSurfacesProviderDim2Bin100 =
        Experimental::KdtSurfacesProvider<2u, 100u>;

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
        .def(py::init<const GeometryContext&,
                      const std::vector<std::shared_ptr<Surface>>&,
                      const std::array<AxisDirection, 2u>&>())
        .def("surfaces", py::overload_cast<const RangeXDDim2&>(
                             &KdtSurfacesDim2Bin100::surfaces, py::const_));

    py::class_<KdtSurfacesProviderDim2Bin100, Experimental::ISurfacesProvider,
               std::shared_ptr<KdtSurfacesProviderDim2Bin100>>(
        m, "KdtSurfacesProviderDim2Bin100")
        .def(py::init<std::shared_ptr<KdtSurfacesDim2Bin100>, const Extent&>());
  }

  {
    using RangeXDDim3 = RangeXD<3u, double>;

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
    py::class_<Experimental::IExternalStructureBuilder,
               std::shared_ptr<Experimental::IExternalStructureBuilder>>(
        m, "IExternalStructureBuilder");

    auto vsBuilder =
        py::class_<VolumeStructureBuilder,
                   Experimental::IExternalStructureBuilder,
                   std::shared_ptr<VolumeStructureBuilder>>(
            m, "VolumeStructureBuilder")
            .def(py::init([](const VolumeStructureBuilder::Config& config,
                             const std::string& name, Logging::Level level) {
              return std::make_shared<VolumeStructureBuilder>(
                  config, getDefaultLogger(name, level));
            }));

    auto vsConfig =
        py::class_<VolumeStructureBuilder::Config>(vsBuilder, "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(vsConfig, boundsType, boundValues, transform, auxiliary);
  }

  {
    py::class_<Experimental::IGeometryIdGenerator,
               std::shared_ptr<Experimental::IGeometryIdGenerator>>(
        m, "IGeometryIdGenerator");

    auto geoIdGen =
        py::class_<Experimental::GeometryIdGenerator,
                   Experimental::IGeometryIdGenerator,
                   std::shared_ptr<Experimental::GeometryIdGenerator>>(
            m, "GeometryIdGenerator")
            .def(py::init([](Experimental::GeometryIdGenerator::Config& config,
                             const std::string& name, Logging::Level level) {
              return std::make_shared<Experimental::GeometryIdGenerator>(
                  config, getDefaultLogger(name, level));
            }));

    auto geoIdGenConfig = py::class_<Experimental::GeometryIdGenerator::Config>(
                              geoIdGen, "Config")
                              .def(py::init<>());
    ACTS_PYTHON_STRUCT(geoIdGenConfig, containerMode, containerId,
                       resetSubCounters, overrideExistingIds);
  }

  {
    // Put them together to a detector volume
    py::class_<Experimental::IDetectorComponentBuilder,
               std::shared_ptr<Experimental::IDetectorComponentBuilder>>(
        m, "IDetectorComponentBuilder");

    auto dvBuilder =
        py::class_<DetectorVolumeBuilder,
                   Experimental::IDetectorComponentBuilder,
                   std::shared_ptr<DetectorVolumeBuilder>>(
            m, "DetectorVolumeBuilder")
            .def(py::init([](const DetectorVolumeBuilder::Config& config,
                             const std::string& name, Logging::Level level) {
              return std::make_shared<DetectorVolumeBuilder>(
                  config, getDefaultLogger(name, level));
            }))
            .def("construct", &DetectorVolumeBuilder::construct);

    auto dvConfig =
        py::class_<DetectorVolumeBuilder::Config>(dvBuilder, "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(dvConfig, name, internalsBuilder, externalsBuilder,
                       geoIdGenerator, auxiliary);
  }

  {
    // The external volume structure builder
    py::class_<Experimental::IRootVolumeFinderBuilder,
               std::shared_ptr<Experimental::IRootVolumeFinderBuilder>>(
        m, "IRootVolumeFinderBuilder");

    auto irvBuilder =
        py::class_<
            Experimental::IndexedRootVolumeFinderBuilder,
            Experimental::IRootVolumeFinderBuilder,
            std::shared_ptr<Experimental::IndexedRootVolumeFinderBuilder>>(
            m, "IndexedRootVolumeFinderBuilder")
            .def(py::init<std::vector<AxisDirection>>());
  }

  {
    // Cylindrical container builder
    auto ccBuilder =
        py::class_<CylindricalContainerBuilder,
                   Experimental::IDetectorComponentBuilder,
                   std::shared_ptr<CylindricalContainerBuilder>>(
            m, "CylindricalContainerBuilder")
            .def(py::init([](const CylindricalContainerBuilder::Config& config,
                             const std::string& name, Logging::Level level) {
              return std::make_shared<CylindricalContainerBuilder>(
                  config, getDefaultLogger(name, level));
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
        py::class_<CuboidalContainerBuilder,
                   Experimental::IDetectorComponentBuilder,
                   std::shared_ptr<CuboidalContainerBuilder>>(
            m, "CuboidalContainerBuilder")
            .def(py::init([](const CuboidalContainerBuilder::Config& config,
                             const std::string& name, Logging::Level level) {
              return std::make_shared<CuboidalContainerBuilder>(
                  config, getDefaultLogger(name, level));
            }))
            .def("construct", &CuboidalContainerBuilder::construct);

    auto ccConfig =
        py::class_<CuboidalContainerBuilder::Config>(ccBuilder, "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(ccConfig, builders, binning, rootVolumeFinderBuilder,
                       geoIdGenerator, geoIdReverseGen, auxiliary);
  }

  {
    // Detector builder
    auto dBuilder =
        py::class_<DetectorBuilder, std::shared_ptr<DetectorBuilder>>(
            m, "DetectorBuilder")
            .def(py::init([](const DetectorBuilder::Config& config,
                             const std::string& name, Logging::Level level) {
              return std::make_shared<DetectorBuilder>(
                  config, getDefaultLogger(name, level));
            }))
            .def("construct", &DetectorBuilder::construct);

    auto dConfig = py::class_<DetectorBuilder::Config>(dBuilder, "Config")
                       .def(py::init<>());
    ACTS_PYTHON_STRUCT(dConfig, name, builder, geoIdGenerator,
                       materialDecorator, auxiliary);
  }
}

}  // namespace ActsPython
