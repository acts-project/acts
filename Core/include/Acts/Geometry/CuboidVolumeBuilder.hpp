// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include <memory>
#include <vector>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ITrackingVolumeBuilder.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class TrackingVolume;
class VolumeBounds;
class RectangleBounds;
class ISurfaceMaterial;
class IVolumeMaterial;
class DetectorElementBase;
class PlaneSurface;

/// @brief This class builds a box detector with a configurable amount of
/// surfaces in it. The idea is to allow a quick configuration of a detector for
/// mostly unit test applications. Therefore this class does not demand to be a
/// universal construction factory but a raw first draft of the idea of factory
/// that may be extended in the future.
class CuboidVolumeBuilder : public ITrackingVolumeBuilder {
 public:
  /// @brief This struct stores the data for the construction of a single
  /// PlaneSurface
  struct SurfaceConfig {
    // Center position
    Vector3D position;
    // Rotation
    RotationMatrix3D rotation = RotationMatrix3D::Identity();
    // Bounds
    std::shared_ptr<const RectangleBounds> rBounds = nullptr;
    // Attached material
    std::shared_ptr<const ISurfaceMaterial> surMat = nullptr;
    // Thickness
    double thickness = 0.;
    // Constructor function for optional detector elements
    // Arguments are transform, rectangle bounds and thickness.
    std::function<DetectorElementBase*(std::shared_ptr<const Transform3D>,
                                       std::shared_ptr<const RectangleBounds>,
                                       double)>
        detElementConstructor;
  };

  /// @brief This struct stores the data for the construction of a PlaneLayer
  /// that has a single PlaneSurface encapsulated
  struct LayerConfig {
    // Configuration of the surface
    SurfaceConfig surfaceCfg;
    // Encapsulated surface
    std::shared_ptr<const PlaneSurface> surface = nullptr;
    // Boolean flag if layer is active
    bool active = false;
    // Bins in Y direction
    size_t binsY = 1;
    // Bins in Z direction
    size_t binsZ = 1;
  };

  /// @brief This struct stores the data for the construction of a cuboid
  /// TrackingVolume with a given number of PlaneLayers
  struct VolumeConfig {
    // Center position
    Vector3D position;
    // Lengths in x,y,z
    Vector3D length;
    // Configurations of its layers
    std::vector<LayerConfig> layerCfg;
    // Stored layers
    std::vector<std::shared_ptr<const Layer>> layers;
    // Configurations of confined volumes
    std::vector<VolumeConfig> volumeCfg;
    // Stored confined volumes
    std::vector<std::shared_ptr<TrackingVolume>> trackingVolumes;
    // Name of the volume
    std::string name = "Volume";
    // Material
    std::shared_ptr<const IVolumeMaterial> volumeMaterial = nullptr;
  };

  /// @brief This struct stores the configuration of the tracking geometry
  struct Config {
    // Center position
    Vector3D position = Vector3D(0., 0., 0.);
    // Length in x,y,z
    Vector3D length = Vector3D(0., 0., 0.);
    // Configuration of its volumes
    std::vector<VolumeConfig> volumeCfg = {};
  };

  /// @brief Default constructor without a configuration
  CuboidVolumeBuilder() = default;

  /// @brief Constructor that sets the config
  ///
  /// @param [in] cfg Configuration of the detector
  CuboidVolumeBuilder(Config& cfg) : m_cfg(cfg) {}

  /// @brief Setter of the config
  ///
  /// @param [in] cfg Configuration that is set
  void setConfig(Config& cfg) { m_cfg = cfg; }

  /// @brief This function creates a surface with a given configuration. A
  /// detector element is attached if the template parameter is non-void.
  ///
  /// @param [in] gctx the geometry context for this building
  /// @param [in] cfg Configuration of the surface
  ///
  /// @return Pointer to the created surface
  std::shared_ptr<const PlaneSurface> buildSurface(
      const GeometryContext& gctx, const SurfaceConfig& cfg) const;

  /// @brief This function creates a layer with a surface encaspulated with a
  /// given configuration. The surface gets a detector element attached if the
  /// template parameter is non-void.
  ///
  /// @param [in] gctx the geometry context for this building
  /// @param [in, out] cfg Configuration of the layer and the surface
  ///
  /// @return Pointer to the created layer
  std::shared_ptr<const Layer> buildLayer(const GeometryContext& gctx,
                                          LayerConfig& cfg) const;

  /// @brief This function creates a TrackingVolume with a configurable number
  /// of layers and surfaces. Each surface gets a detector element attached if
  /// the template parameter is non-void.
  ///
  /// @param [in] gctx the geometry context for this building
  /// @param [in, out] cfg Configuration of the TrackingVolume
  ///
  /// @return Pointer to the created TrackingVolume
  std::shared_ptr<TrackingVolume> buildVolume(const GeometryContext& gctx,
                                              VolumeConfig& cfg) const;

  /// @brief This function evaluates the minimum and maximum of the binning as
  /// given by the configurations of the surfaces and layers. The ordering
  /// depends on the binning value specified in the configuration of the volume.
  ///
  /// @param [in] gctx the geometry context for this building
  /// @param [in] cfg Container with the given surfaces and layers
  ///
  /// @return Pair containing the minimum and maximum along the binning
  /// direction
  std::pair<double, double> binningRange(const GeometryContext& gctx,
                                         const VolumeConfig& cfg) const;

  void
  sortVolumes(std::vector<std::pair<TrackingVolumePtr, Vector3D>>& tapVec,
              BinningValue bValue) const;
};

  /// @brief This function builds a world TrackingVolume based on a given
  /// configuration
  ///
  /// @param [in] gctx the geometry context for this building
  ///
  /// @return Pointer to the created TrackingGeometry
  std::shared_ptr<TrackingVolume> trackingVolume(
      const GeometryContext& gctx,
      std::shared_ptr<const TrackingVolume> /*unused*/,
      std::shared_ptr<const VolumeBounds> /*unused*/) const override;

 private:
  /// Configuration of the world volume
  Config m_cfg;
};
<<<<<<< HEAD
=======

std::shared_ptr<const PlaneSurface>
CuboidVolumeBuilder::buildSurface(const SurfaceConfig& cfg) const
{
  PlaneSurface* surface;

  // Build transformation
  Transform3D trafo(Transform3D::Identity() * cfg.rotation);
  trafo.translation() = cfg.position;

  // Create and store surface
  if (cfg.detElementConstructor) {
    surface = new PlaneSurface(cfg.rBounds,
                               *(cfg.detElementConstructor(std::make_tuple(
                                   std::make_shared<const Transform3D>(trafo),
                                   cfg.rBounds,
                                   cfg.thickness))));
  } else {
    surface = new PlaneSurface(std::make_shared<const Transform3D>(trafo),
                               cfg.rBounds);
  }
  surface->setAssociatedMaterial(cfg.surMat);
  return std::shared_ptr<const PlaneSurface>(surface);
}

std::shared_ptr<const Layer>
CuboidVolumeBuilder::buildLayer(LayerConfig& cfg) const
{
  // Build the surface
  if (cfg.surface == nullptr) {
    cfg.surface = buildSurface(cfg.surfaceCfg);
  }
  // Build transformation centered at the surface position
  Transform3D trafo(Transform3D::Identity() * cfg.surfaceCfg.rotation);
  trafo.translation() = cfg.surfaceCfg.position;

  LayerCreator::Config lCfg;
  lCfg.surfaceArrayCreator
      = std::shared_ptr<const SurfaceArrayCreator>(new SurfaceArrayCreator());
  LayerCreator layerCreator(lCfg);

  return layerCreator.planeLayer({cfg.surface},
                                 cfg.binsY,
                                 cfg.binsZ,
                                 BinningValue::binX,
                                 boost::none,
                                 std::make_shared<const Transform3D>(trafo));
}

std::pair<double, double>
CuboidVolumeBuilder::binningRange(const VolumeConfig& cfg) const
{
  // Construct return value
  std::pair<double, double> minMax = std::make_pair(
      std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
  for (const auto& layercfg : cfg.layerCfg) {
    // Test if new extreme is found and set it
    if (layercfg.surfaceCfg.position.x() - 1. * units::_um < minMax.first) {
      minMax.first = layercfg.surfaceCfg.position.x() - 1. * units::_um;
    }
    if (layercfg.surfaceCfg.position.x() + 1. * units::_um > minMax.second) {
      minMax.second = layercfg.surfaceCfg.position.x() + 1. * units::_um;
    }
  }
  return minMax;
}

std::shared_ptr<TrackingVolume>
CuboidVolumeBuilder::buildVolume(VolumeConfig& cfg) const
{
  // Build transformation
  Transform3D trafo(Transform3D::Identity());
  trafo.translation() = cfg.position;
  // Set bounds
  auto bounds = std::make_shared<const CuboidVolumeBounds>(
      cfg.length.x() * 0.5, cfg.length.y() * 0.5, cfg.length.z() * 0.5);

  //~ if (cfg.layerCfg.empty()) {
  //~ // Build dummy layer if no layer is given (tmp solution)
  //~ SurfaceConfig sCfg;
  //~ sCfg.position = cfg.position;
  //~ // Rotation of the surfaces
  //~ double   rotationAngle = M_PI * 0.5;
  //~ Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
  //~ Vector3D yPos(0., 1., 0.);
  //~ Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
  //~ sCfg.rotation.col(0) = xPos;
  //~ sCfg.rotation.col(1) = yPos;
  //~ sCfg.rotation.col(2) = zPos;
  //~ // Bounds
  //~ sCfg.rBounds = std::make_shared<const RectangleBounds>(
  //~ RectangleBounds(cfg.length.y() * 0.5, cfg.length.z() * 0.5));

  //~ LayerConfig lCfg;
  //~ lCfg.surfaceCfg     = sCfg;
  //~ lCfg.layerThickness = 1. * units::_mm;

  //~ cfg.layerCfg.push_back(lCfg);
  //~ }

  // Gather the layers
  LayerVector layVec;
  if (cfg.layers.empty()) {
    cfg.layers.reserve(cfg.layerCfg.size());

    for (auto& layerCfg : cfg.layerCfg) {
      cfg.layers.push_back(buildLayer(layerCfg));
      layVec.push_back(cfg.layers.back());
    }
  } else {
    for (auto& lay : cfg.layers) {
      layVec.push_back(lay);
    }
  }

  // Build confined volumes
  if (cfg.trackingVolumes.empty())
    for (VolumeConfig vc : cfg.volumeCfg)
      cfg.trackingVolumes.push_back(buildVolume(vc));

  std::shared_ptr<TrackingVolume> trackVolume;
  if (layVec.empty()) {
    // Build TrackingVolume
    trackVolume
        = TrackingVolume::create(std::make_shared<const Transform3D>(trafo),
                                 bounds,
                                 cfg.material,
                                 nullptr,
                                 {},
                                 cfg.trackingVolumes,
                                 {},
                                 cfg.name);
  } else {
    // Build layer array
    std::pair<double, double> minMax = binningRange(cfg);
    LayerArrayCreator layArrCreator(
        getDefaultLogger("LayerArrayCreator", Logging::INFO));
    std::unique_ptr<const LayerArray> layArr(
        layArrCreator.layerArray(layVec,
                                 minMax.first,
                                 minMax.second,
                                 BinningType::arbitrary,
                                 BinningValue::binX));

    // Build TrackingVolume
    trackVolume
        = TrackingVolume::create(std::make_shared<const Transform3D>(trafo),
                                 bounds,
                                 cfg.material,
                                 std::move(layArr),
                                 layVec,
                                 cfg.trackingVolumes,
                                 {},
                                 cfg.name);
  }

  return trackVolume;
}

MutableTrackingVolumePtr
    CuboidVolumeBuilder::trackingVolume(TrackingVolumePtr /*unused*/,
                                        VolumeBoundsPtr /*unused*/) const
{
  // Build volumes
  if (cfg.volumes.empty()) {
    cfg.volumes.reserve(cfg.volumeCfg.size());
    for (VolumeConfig volCfg : cfg.volumeCfg) {
      cfg.volumes.push_back(buildVolume(volCfg));
    }
  }

  // Glue volumes
  for (unsigned int i = 0; i < cfg.volumes.size() - 1; i++) {
    cfg.volumes[i + 1]->glueTrackingVolume(BoundarySurfaceFace::negativeFaceYZ,
                                           cfg.volumes[i],
                                           BoundarySurfaceFace::positiveFaceYZ);
    cfg.volumes[i]->glueTrackingVolume(BoundarySurfaceFace::positiveFaceYZ,
                                       cfg.volumes[i + 1],
                                       BoundarySurfaceFace::negativeFaceYZ);
  }

  // Translation
  Transform3D trafo(Transform3D::Identity());
  trafo.translation() = m_cfg.position;

  // Size of the volume
  auto volume = std::make_shared<const CuboidVolumeBounds>(
      m_cfg.length.x() * 0.5, m_cfg.length.y() * 0.5, m_cfg.length.z() * 0.5);

  // Build vector of confined volumes
  std::vector<std::pair<TrackingVolumePtr, Vector3D>> tapVec;
  tapVec.reserve(cfg.volumeCfg.size());
  for (auto& tVol : cfg.volumes) {
    tapVec.push_back(std::make_pair(tVol, tVol->center()));
  }

  // Set bin boundaries along binning
  std::vector<float> binBoundaries;
  binBoundaries.push_back(cfg.volumes[0]->center().x()
                          - cfg.volumeCfg[0].length.x() * 0.5);
  for (size_t i = 0; i < cfg.volumes.size(); i++) {
    binBoundaries.push_back(cfg.volumes[i]->center().x()
                            + cfg.volumeCfg[i].length.x() * 0.5);
  }

  // Build binning
  BinningData binData(BinningOption::open, BinningValue::binX, binBoundaries);
  std::unique_ptr<const BinUtility> bu(new BinUtility(binData));

  // Build TrackingVolume array
  std::shared_ptr<const TrackingVolumeArray> trVolArr(
      new BinnedArrayXD<TrackingVolumePtr>(tapVec, std::move(bu)));

  // Create world volume
  MutableTrackingVolumePtr mtvp(TrackingVolume::create(
      std::make_shared<const Transform3D>(trafo), volume, trVolArr, "World"));

  mtvp->sign(GeometrySignature::Global);
  return mtvp;
}
>>>>>>> Boundary connections added
}  // namespace Acts
