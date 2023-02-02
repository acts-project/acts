// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/KDTreeTrackingGeometryBuilder.hpp"

#include "Acts/Geometry/CylinderLayer.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <functional>
#include <iostream>

namespace {

/// @brief Anonymous helper function to measure a surface
///
/// @tparam kDIM the dimension of the tree
/// @param gtcx the Geometry context
/// @param surface the Surface to be measured
/// @param bValues the binning values for the KDTree dimensions
/// @param polyhedron use polyhedron otherwise binning position
///
/// @return a newly created KDTree
template <size_t kDIM>
std::array<Acts::ActsScalar, kDIM> measure(
    const Acts::GeometryContext& gctx, const Acts::Surface& surface,
    const std::array<Acts::BinningValue, kDIM>& bValues, bool polyhedron) {
  std::array<Acts::ActsScalar, kDIM> ms;
  if (polyhedron) {
    auto ext = surface.polyhedronRepresentation(gctx, 1u).extent();
    for (std::size_t idim = 0; idim < kDIM; ++idim) {
      ms[idim] = ext.medium(bValues[idim]);
    }
  } else {
    for (std::size_t idim = 0; idim < kDIM; ++idim) {
      auto binningPosition = surface.binningPosition(gctx, bValues[idim]);
      ms[idim] = Acts::VectorHelpers::cast(binningPosition, bValues[idim]);
    }
  }
  return ms;
}

/// @brief Anonymous helper function to generate a KDTree
///
/// @tparam kDIM the dimension of the tree
/// @param gtcx the Geometry context
/// @param surfaces the Surfaces
/// @param bValues the binning values for the KDTree dimensions
/// @param polyhedron use polyhedron otherwise binning position
///
/// @return a newly created KDTree
template <size_t kDIM>
Acts::KDTree<kDIM, std::shared_ptr<Acts::Surface>, Acts::ActsScalar, std::array,
             100>
generateKDTree(const Acts::GeometryContext& gctx,
               const std::vector<std::shared_ptr<Acts::Surface>>& surfaces,
               const std::array<Acts::BinningValue, kDIM>& bValues,
               bool polyhedron) {
  // The measured surface defintion
  using MeasuredSurface = std::pair<std::array<Acts::ActsScalar, kDIM>,
                                    std::shared_ptr<Acts::Surface>>;
  // Prepare all the surfaces
  std::vector<MeasuredSurface> surfacesMeasured;
  surfacesMeasured.reserve(surfaces.size());
  for (const auto& s : surfaces) {
    auto ext = s->polyhedronRepresentation(gctx, 1u).extent();
    surfacesMeasured.push_back(
        MeasuredSurface{measure<kDIM>(gctx, *s, bValues, polyhedron), s});
  }

  // Create the KDTree
  return Acts::KDTree<kDIM, std::shared_ptr<Acts::Surface>, Acts::ActsScalar,
                      std::array, 100>(std::move(surfacesMeasured));
}

}  // namespace

Acts::KDTreeTrackingGeometryBuilder::KDTreeTrackingGeometryBuilder(
    const Acts::KDTreeTrackingGeometryBuilder::Config& cfg,
    std::unique_ptr<const Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {
  // Inconsistency - bail out
  if (m_cfg.bValues.size() == 0u or m_cfg.bValues.size() > 3u) {
    throw std::invalid_argument(
        "KDTreeTrackingGeometryBuilder works only in binnings 1D to 3D.");
  }
  // Legacy harmonization
  m_cfg.protoDetector.harmonize();
}

std::unique_ptr<const Acts::TrackingGeometry>
Acts::KDTreeTrackingGeometryBuilder::trackingGeometry(
    const GeometryContext& gctx) const {
  // Fill the KDTree accorindly to the binning setup
  SurfaceKDT surfaceKDT;
  if (m_cfg.bValues.size() == 1) {
    surfaceKDT.kdt1 = generateKDTree<1u>(
        gctx, m_cfg.surfaces, {m_cfg.bValues[0u]}, m_cfg.polyhedronMeasurement);
  } else if (m_cfg.bValues.size() == 2u) {
    surfaceKDT.kdt2 = generateKDTree<2u>(gctx, m_cfg.surfaces,
                                         {m_cfg.bValues[0u], m_cfg.bValues[1u]},
                                         m_cfg.polyhedronMeasurement);
  } else {
    surfaceKDT.kdt3 = generateKDTree<3u>(
        gctx, m_cfg.surfaces,
        {m_cfg.bValues[0u], m_cfg.bValues[1u], m_cfg.bValues[2u]},
        m_cfg.polyhedronMeasurement);
  }

  // Walk through the proto builder
  auto protoWorld = m_cfg.protoDetector.worldVolume;

  KDTreeTrackingGeometryBuilder::Cache cCache;

  auto worldTrackingVolume =
      translateVolume(cCache, gctx, surfaceKDT, protoWorld);

  ACTS_INFO("Retrieved " << cCache.surfaceCounter
                         << " surfaces from the KDTree.");

  // return the geometry to the service
  return std::make_unique<const Acts::TrackingGeometry>(worldTrackingVolume);
}

std::shared_ptr<Acts::TrackingVolume>
Acts::KDTreeTrackingGeometryBuilder::translateVolume(
    Cache& cCache, const GeometryContext& gctx, const SurfaceKDT& kdt,
    const ProtoVolume& ptVolume, const std::string& indent) const {
  ACTS_DEBUG(indent << "Processing ProtoVolume: " << ptVolume.name);
  std::vector<std::shared_ptr<const TrackingVolume>> translatedVolumes = {};

  // Simple gap volume
  if (not ptVolume.container.has_value()) {
    ACTS_VERBOSE(indent << "> empty volume to be built");
    MutableTrackingVolumeVector mtv = {};
    auto tVolume = m_cfg.trackingVolumeHelper->createGapTrackingVolume(
        gctx, mtv, nullptr, ptVolume.extent, 0, Surface::SurfaceType::Other,
        ptVolume.name);
    ACTS_DEBUG(indent << "> translated into gap volume bounds: "
                      << tVolume->volumeBounds());
    return tVolume;
  } else {
    // Get the container information
    auto& cts = ptVolume.container.value();

    // Container information must be present
    if (cts.constituentVolumes.empty()) {
      throw std::invalid_argument(
          "KDTreeTrackingGeometryBuilder: no consituents given.");
    }

    // This volume is a volume container
    if (not cts.layerContainer) {
      ACTS_VERBOSE(indent << "> volume container with "
                          << cts.constituentVolumes.size() << " constituents.");
      for (auto& cVolume : cts.constituentVolumes) {
        auto dtVolume = translateVolume(cCache, gctx, kdt, cVolume,
                                        indent + m_cfg.hierarchyIndent);
        translatedVolumes.push_back(dtVolume);
      }
      auto tVolume = m_cfg.trackingVolumeHelper->createContainerTrackingVolume(
          gctx, translatedVolumes);
      ACTS_DEBUG(indent << "> translated into container volume bounds: "
                        << tVolume->volumeBounds());
      return tVolume;
    } else {
      // This volume is a layer container
      std::vector<std::shared_ptr<const Layer>> layers = {};
      ACTS_VERBOSE(indent << "> layer container with "
                          << cts.constituentVolumes.size() << " layers.");
      for (auto& plVolume : cts.constituentVolumes) {
        if (plVolume.internal.has_value()) {
          // Retrieve the layer surfaces from the KD tree
          std::vector<std::shared_ptr<const Surface>> layerSurfaces;
          if (m_cfg.bValues.size() == 1u) {
            RangeXD<1u, ActsScalar> range;
            range[0u] = plVolume.extent.range(m_cfg.bValues[0u]);
            ACTS_VERBOSE(indent + ">> looking in range = " << range.toString());
            auto kdtSurfaces = kdt.kdt1.value().rangeSearchWithKey(range);
            layerSurfaces.reserve(kdtSurfaces.size());
            for (const auto& kdts : kdtSurfaces) {
              layerSurfaces.push_back(kdts.second);
            }
          } else if (m_cfg.bValues.size() == 2u) {
            RangeXD<2u, ActsScalar> range;
            range[0u] = plVolume.extent.range(m_cfg.bValues[0u]);
            range[1u] = plVolume.extent.range(m_cfg.bValues[1u]);
            ACTS_VERBOSE(indent + ">> looking in range = " << range.toString());
            auto kdtSurfaces = kdt.kdt2.value().rangeSearchWithKey(range);
            layerSurfaces.reserve(kdtSurfaces.size());
            for (const auto& kdts : kdtSurfaces) {
              layerSurfaces.push_back(kdts.second);
            }

          } else if (m_cfg.bValues.size() == 3u) {
            RangeXD<3u, ActsScalar> range;
            range[0u] = plVolume.extent.range(m_cfg.bValues[0u]);
            range[1u] = plVolume.extent.range(m_cfg.bValues[1u]);
            range[2u] = plVolume.extent.range(m_cfg.bValues[2u]);
            ACTS_VERBOSE(indent + ">> looking in range = " << range.toString());
            auto kdtSurfaces = kdt.kdt3.value().rangeSearchWithKey(range);
            layerSurfaces.reserve(kdtSurfaces.size());
            for (const auto& kdts : kdtSurfaces) {
              layerSurfaces.push_back(kdts.second);
            }
          }
          ACTS_VERBOSE(indent + ">> found " << layerSurfaces.size()
                                            << " surfaces in the KDTree.");

          layers.push_back(translateLayer(cCache, gctx, layerSurfaces, plVolume,
                                          indent + m_cfg.hierarchyIndent));
        } else {
          ACTS_WARNING(indent << "> layer type volume has no internal "
                                 "description, layer not built.");
        }
      }
      // Create a new tracking volume with those layers
      auto tVolume = m_cfg.trackingVolumeHelper->createTrackingVolume(
          gctx, layers, {}, nullptr, ptVolume.extent, ptVolume.name);
      ACTS_DEBUG(indent << "> translated into bounds: "
                        << tVolume->volumeBounds());
      return tVolume;
    }
  }
  return nullptr;
}

/// @return a new tracking volume
std::shared_ptr<const Acts::Layer>
Acts::KDTreeTrackingGeometryBuilder::translateLayer(
    Cache& cCache, const GeometryContext& gctx,
    const std::vector<std::shared_ptr<const Surface>>& layerSurfaces,
    const ProtoVolume& plVolume, const std::string& indent) const {
  ACTS_DEBUG(indent + "Processing ProtoVolume: " << plVolume.name);

  // This is only called if the volume has internal structure
  auto& its = plVolume.internal.value();
  cCache.surfaceCounter += layerSurfaces.size();
  std::shared_ptr<const Acts::Layer> tLayer = nullptr;

  if (layerSurfaces.size() == 1u) {
    auto surface = layerSurfaces[0u];
    const auto& transform = surface->transform(gctx);
    if (its.layerType == Acts::Surface::SurfaceType::Cylinder) {
      ACTS_VERBOSE(indent +
                   ">> creating cylinder layer from a single surface.");
      // Get the bounds
      auto cylinderBounds =
          dynamic_cast<const CylinderBounds*>(&(surface->bounds()));
      auto cylinderBoundsClone =
          std::make_shared<const CylinderBounds>(*cylinderBounds);
      auto cylinderLayer =
          CylinderLayer::create(transform, cylinderBoundsClone, nullptr, 1.);
      cylinderLayer->assignSurfaceMaterial(surface->surfaceMaterialSharedPtr());
      tLayer = cylinderLayer;
    } else if (its.layerType == Acts::Surface::SurfaceType::Disc) {
      ACTS_VERBOSE(indent +
                   ">> creating cylinder layer from a single surface.");
      // Get the bounds
      auto radialBounds =
          dynamic_cast<const RadialBounds*>(&(surface->bounds()));
      auto radialBoundsClone =
          std::make_shared<const RadialBounds>(*radialBounds);
      auto discLayer =
          DiscLayer::create(transform, radialBoundsClone, nullptr, 1.);
      discLayer->assignSurfaceMaterial(surface->surfaceMaterialSharedPtr());
      tLayer = discLayer;
    } else if (its.layerType == Acts::Surface::SurfaceType::Plane) {
      ACTS_VERBOSE(indent + ">> creating plane layer from a single surface.");
      // auto radialBounds =
      //     dynamic_cast<const RadialBounds*>(&(surface->bounds()));
      // auto radialBoundsClone =
      //     std::make_shared<const RadialBounds>(*radialBounds);
      // auto discLayer =
      //     DiscLayer::create(transform, radialBoundsClone, nullptr, 1.);
      // discLayer->assignSurfaceMaterial(surface->surfaceMaterialSharedPtr());
      // tLayer = discLayer;
    } else {
      throw std::invalid_argument(
          "KDTreeTrackingGeometryBuilder: layer type is neither cylinder nor "
          "disk.");
    }

  } else if (layerSurfaces.size() > 1u) {
    Acts::BinningType bType0 = Acts::equidistant;
    Acts::BinningType bType1 = Acts::equidistant;
    std::size_t bins0 = 0;
    std::size_t bins1 = 0;
    // In case explicit binning is given
    if (its.surfaceBinning.size() == 2u) {
      bType0 = its.surfaceBinning[0u].type;
      bType1 = its.surfaceBinning[1u].type;
      // In case explicit bin numbers are given in addition
      if (bType0 == Acts::equidistant and bType1 == Acts::equidistant and
          its.surfaceBinning[0u].bins() > 1u and
          its.surfaceBinning[1u].bins() > 1u) {
        bins0 = its.surfaceBinning[0u].bins();
        bins1 = its.surfaceBinning[1u].bins();
        ACTS_VERBOSE(indent + ">> binning provided externally to be "
                     << bins0 << " x " << bins1 << ".");
      }
    }

    Acts::ProtoLayer pLayer(gctx, layerSurfaces);
    pLayer.envelope = plVolume.extent.envelope();
    if (its.layerType == Acts::Surface::SurfaceType::Cylinder) {
      ACTS_VERBOSE(indent + ">> creating cylinder layer with "
                   << layerSurfaces.size() << " surfaces.");
      // Forced equidistant or auto-binned
      tLayer = (bins0 * bins1 > 0)
                   ? m_cfg.layerCreator->cylinderLayer(gctx, layerSurfaces,
                                                       bins0, bins1, pLayer)
                   : m_cfg.layerCreator->cylinderLayer(gctx, layerSurfaces,
                                                       bType0, bType1, pLayer);

    } else if (its.layerType == Acts::Surface::SurfaceType::Disc) {
      ACTS_VERBOSE(indent + ">> creating disc layer with "
                   << layerSurfaces.size() << " surfaces.");
      // Forced equidistant or auto-binned
      tLayer = (bins0 * bins1 > 0)
                   ? m_cfg.layerCreator->discLayer(gctx, layerSurfaces, bins0,
                                                   bins1, pLayer)

                   : m_cfg.layerCreator->discLayer(gctx, layerSurfaces, bType0,
                                                   bType1, pLayer);
    } else if (its.layerType == Acts::Surface::SurfaceType::Plane) {
      ACTS_VERBOSE(indent + ">> creating plane layer with "
                   << layerSurfaces.size() << " surfaces.");
      // Forced equidistant or auto-binned
      tLayer =
          m_cfg.layerCreator->planeLayer(gctx, layerSurfaces, 1u, 1u, binZ);
    } else {
      throw std::invalid_argument(
          "KDTreeTrackingGeometryBuilder: layer type is neither cylinder nor "
          "disk.");
    }
  }
  /** @TODO PlaneLayer needs representing volume
  if (tLayer != nullptr and tLayer->representingVolume() != nullptr) {
    ACTS_DEBUG(indent << "> translated into layer bounds: "
                      << tLayer->representingVolume()->volumeBounds());
  } else {
    throw std::runtime_error(
        "KDTreeTrackingGeometryBuilder: layer was not built.");
  }
  */
  return tLayer;
}
