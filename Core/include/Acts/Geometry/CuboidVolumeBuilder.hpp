// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ITrackingVolumeBuilder.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace Acts {

class TrackingVolume;
class VolumeBounds;
class RectangleBounds;
class ISurfaceMaterial;
class IVolumeMaterial;
class Surface;
class Layer;

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
    /// Center position of the surface
    Vector3 position{};
    // Rotation
    /// Rotation matrix defining surface orientation
    RotationMatrix3 rotation = RotationMatrix3::Identity();
    // Bounds
    /// Rectangle bounds defining surface dimensions
    std::shared_ptr<const RectangleBounds> rBounds = nullptr;
    // Attached material
    /// Surface material description
    std::shared_ptr<const ISurfaceMaterial> surMat = nullptr;
    // Thickness
    /// Material thickness of the surface
    double thickness = 0.;
    /// Constructor function for optional detector elements
    /// Arguments are transform, rectangle bounds and thickness.
    std::function<SurfacePlacementBase*(
        const Transform3&, std::shared_ptr<const RectangleBounds>, double)>
        detElementConstructor;
  };

  /// @brief This struct stores the data for the construction of a PlaneLayer
  struct LayerConfig {
    // Configuration of the surface
    /// Configuration objects for surfaces in this layer
    std::vector<SurfaceConfig> surfaceCfg;
    // Encapsulated surface
    /// Pre-built surfaces to be included in this layer
    std::vector<std::shared_ptr<const Surface>> surfaces;
    // Boolean flag if layer is active
    /// Flag indicating whether layer participates in track finding
    bool active = false;
    // Bins in Y direction
    /// Number of bins for surface arrangement in Y direction
    std::size_t binsY = 1;
    // Bins in Z direction
    /// Number of bins for surface arrangement in Z direction
    std::size_t binsZ = 1;
    // Envelope in X
    /// Envelope extensions in X direction [negative, positive]
    std::array<double, 2u> envelopeX{0, 0};
    // Envelope in Y
    /// Envelope extensions in Y direction [negative, positive]
    std::array<double, 2u> envelopeY{0, 0};
    // Envelope in Z
    /// Envelope extensions in Z direction [negative, positive]
    std::array<double, 2u> envelopeZ{0, 0};
    // An optional rotation for this
    /// Optional rotation transformation for this layer
    std::optional<RotationMatrix3> rotation{std::nullopt};
    // Dimension for the binning
    /// Axis direction for surface binning within the layer
    AxisDirection binningDimension = AxisDirection::AxisX;
  };

  /// @brief This struct stores the data for the construction of a cuboid
  /// TrackingVolume with a given number of PlaneLayers
  struct VolumeConfig {
    /// Center position of the tracking volume
    Vector3 position{};
    /// Dimensions of the volume in x, y, z directions
    Vector3 length{};
    /// Configuration objects for layers within this volume
    std::vector<LayerConfig> layerCfg;
    /// Pre-built layers to be included in this volume
    std::vector<std::shared_ptr<const Layer>> layers;
    /// Configuration objects for sub-volumes within this volume
    std::vector<VolumeConfig> volumeCfg;
    /// Pre-built sub-volumes contained within this volume
    std::vector<std::shared_ptr<TrackingVolume>> trackingVolumes;
    /// Identifier name for this tracking volume
    std::string name = "Volume";
    /// Volume material description for this tracking volume
    std::shared_ptr<const IVolumeMaterial> volumeMaterial = nullptr;
    /// Axis direction for layer binning within the volume
    AxisDirection binningDimension = AxisDirection::AxisX;
  };

  /// @brief This struct stores the configuration of the tracking geometry
  struct Config {
    /// Center position of the world volume
    Vector3 position = Vector3(0., 0., 0.);
    /// Dimensions of the world volume in x, y, z directions
    Vector3 length = Vector3(0., 0., 0.);
    /// Configuration objects for all volumes in the detector
    std::vector<VolumeConfig> volumeCfg = {};
  };

  /// @brief Default constructor without a configuration
  CuboidVolumeBuilder() = default;

  /// @brief Constructor that sets the config
  ///
  /// @param [in] cfg Configuration of the detector
  explicit CuboidVolumeBuilder(Config& cfg) : m_cfg(cfg) {}

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
  std::shared_ptr<const Surface> buildSurface(const GeometryContext& gctx,
                                              const SurfaceConfig& cfg) const;

  /// @brief This function creates a layer with a surface encapsulated with a
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

  /// Sort volumes along a given axis direction
  /// @param tapVec Vector of tracking volume and position pairs to sort
  /// @param bValue Axis direction for sorting
  void sortVolumes(std::vector<std::pair<TrackingVolumePtr, Vector3>>& tapVec,
                   AxisDirection bValue) const;

  /// @brief This function builds a world TrackingVolume based on a given
  /// configuration
  ///
  /// @param [in] gctx the geometry context for this building
  ///
  /// @return Pointer to the created TrackingGeometry
  std::shared_ptr<TrackingVolume> trackingVolume(
      const GeometryContext& gctx,
      std::shared_ptr<const TrackingVolume> /*oppositeVolume*/,
      std::shared_ptr<const VolumeBounds> /*outsideBounds*/) const override;

 private:
  /// Configuration of the world volume
  Config m_cfg;
};
}  // namespace Acts
