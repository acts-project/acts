// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/BoundarySurfaceFace.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ITrackingVolumeHelper.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace Acts {

class Layer;
class TrackingVolume;
class VolumeBounds;
class CuboidVolumeBounds;
class IVolumeMaterial;
class ILayerArrayCreator;
class ITrackingVolumeArrayCreator;

/// @class CuboidVolumeHelper
///
/// The concrete implementation for cubic TrackingVolume
/// objects of the ITrackingVolumeCreator interface
///
/// @note this builder is restricted to a telescope along x,y or z
///
class CuboidVolumeHelper : public ITrackingVolumeHelper {
 public:
  /// @struct Config
  /// Nested configuration struct for this CuboidVolumeHelper
  struct Config {
    /// a tool for coherent LayerArray creation
    std::shared_ptr<const ILayerArrayCreator> layerArrayCreator = nullptr;
    /// Helper Tool to create TrackingVolume
    std::shared_ptr<const ITrackingVolumeArrayCreator>
        trackingVolumeArrayCreator = nullptr;
    /// Thickness of passive layers
    double passiveLayerThickness = 1;
    /// Envelope
    ActsScalar envelope = 1;
    /// The Binning value
    BinningValue bValue = binZ;
  };

  /// Constructor
  /// @param cvhConfig is the configuration struct for this builder
  /// @param logger logging instance
  CuboidVolumeHelper(const Config& cvhConfig,
                     std::unique_ptr<const Logger> logger =
                         getDefaultLogger("CuboidVolumeHelper", Logging::INFO));

  ~CuboidVolumeHelper() override = default;

  /// Create a TrackingVolume* from a set of layers and (optional) parameters
  ///
  /// @param gctx is the geometry context for witch the volume is built
  /// @param layers vector of static layers confined by the TrackingVolume
  /// if no bounds or HepTransform is given, they define the size
  /// together with the volume enevlope parameters
  /// @param volumeMaterial material properties for this TrackingVolume
  /// @param volumeBounds: confinement of this TrackingVolume
  /// @param mtvVector (optiona) Vector of confined TrackingVolumes
  /// @param transform (optional) placement of this TrackingVolume
  /// @param volumeName  volume name to be given
  /// @param bType (optional) BinningType - arbitrary(default) or equidistant
  ///
  /// @return shared pointer to a new TrackingVolume
  MutableTrackingVolumePtr createTrackingVolume(
      const GeometryContext& gctx, const LayerVector& layers,
      std::shared_ptr<const IVolumeMaterial> volumeMaterial,
      VolumeBoundsPtr volumeBounds, MutableTrackingVolumeVector mtvVector = {},
      const Transform3& transform = Transform3::Identity(),
      const std::string& volumeName = "UndefinedVolume",
      BinningType bType = arbitrary) const override;

  /// Create a TrackingVolume* from a set of layers and (optional) parameters
  ///
  /// @param gctx is the geometry context for witch the volume is built
  /// @param layers vector of static layers confined by the TrackingVolume
  /// if no bounds or HepTransform is given, they define the size
  /// together with the volume enevlope parameters
  /// @param volumeMaterial material properties for this TrackingVolume
  /// @param mtvVector Vector of confined TrackingVolumes
  /// @param dimensions the dimensions parameters
  /// @param volumeName  volume name to be given
  /// @param bType (optional) BinningType - arbitrary(default) or equidistant
  ///
  /// @return shared pointer to a new TrackingVolume
  MutableTrackingVolumePtr createTrackingVolume(
      const GeometryContext& gctx, const LayerVector& layers,
      MutableTrackingVolumeVector mtvVector,
      std::shared_ptr<const IVolumeMaterial> volumeMaterial,
      const Extent& dimension,
      const std::string& volumeName = "UndefinedVolume",
      BinningType bType = arbitrary) const override;

  /// Create a gap volume from dimensions and
  /// @note this TrackingVolume is restricted to Translation only
  ///
  /// @param [in] gctx the geometry context for this building
  /// @param mtvVector Vector of confined TrackingVolumes
  /// @param volumeMaterial dense material properties for this TrackingVolume
  /// @param dimensions the dimensions parameters
  /// @param materialLayers number of material layers (aequidistant binning)
  /// @param cylinder type of layers
  /// @param volumeName  volume name to be given
  ///
  /// @return shared pointer to a new TrackingVolume
  MutableTrackingVolumePtr createGapTrackingVolume(
      const GeometryContext& gctx, MutableTrackingVolumeVector& mtvVector,
      std::shared_ptr<const IVolumeMaterial> volumeMaterial,
      const Extent& dimension, unsigned int materialLayers,
      Surface::SurfaceType layerType = Surface::SurfaceType::Plane,
      const std::string& volumeName = "UndefinedVolume") const override;

  /// Create a gap volume from dimensions and
  ///
  /// @param [in] gctx the geometry context for this building
  /// @param mtvVector Vector of confined TrackingVolumes
  /// @param volumeMaterial dense material properties for this TrackingVolume
  /// @param dimensions the dimensions parameters
  /// @param layerPositions custom layer positions
  /// @param layerType is the type for the layer to be built
  /// @param volumeName  : volume name to be given
  /// @param bType (optional) BinningType - arbitrary(default) or equidistant
  ///
  /// @return shared pointer to a new TrackingVolume
  MutableTrackingVolumePtr createGapTrackingVolume(
      const GeometryContext& gctx, MutableTrackingVolumeVector& mtvVector,
      std::shared_ptr<const IVolumeMaterial> volumeMaterial,
      const Extent& dimension, const std::vector<double>& layerPositions,
      Surface::SurfaceType layerType = Surface::SurfaceType::Plane,
      const std::string& volumeName = "UndefinedVolume",
      BinningType bType = arbitrary) const override;

  /// Create a container volumes from sub volumes
  ///
  /// @param [in] gctx the geometry context for this building
  /// @param volumes the volumes to be contained
  ///
  /// @note this tool only offers limited functionality of
  /// simply creating a container volume in one direction
  ///
  /// @return shared pointer to a new TrackingVolume
  MutableTrackingVolumePtr createContainerTrackingVolume(
      const GeometryContext& gctx,
      const TrackingVolumeVector& volumes) const override;

  /// Set configuration method
  ///
  /// @param cvhConfig is the configurtion struct assigned
  void setConfiguration(const Config& cvhConfig);

  /// Get configuration method
  Config getConfiguration() const;

  /// Set logging instance
  ///
  /// @param newLogger is the logger isntance to be set
  void setLogger(std::unique_ptr<const Logger> newLogger);

 protected:
  /// Configuration object
  Config m_cfg;

 private:
  /// Private access method to the logging instance
  const Logger& logger() const { return *m_logger; }

  /// the looging instance
  std::unique_ptr<const Logger> m_logger;
};

inline CuboidVolumeHelper::Config CuboidVolumeHelper::getConfiguration() const {
  return m_cfg;
}
}  // namespace Acts
