// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include <memory>
#include <vector>
#include "Acts/Tools/ITrackingVolumeBuilder.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class TrackingVolume;
class VolumeBounds;
class Material;
class RectangleBounds;
class SurfaceMaterial;
class DetectorElementBase;
class PlaneSurface;

/// @brief This class builds a box detector with a configurable amount of
/// surfaces in it. The idea is to allow a quick configuration of a detector for
/// mostly unit test applications. Therefore this class does not demand to be a
/// universal construction factory but a raw first draft of the idea of factory
/// that may be extended in the future.
class CuboidVolumeBuilder : public ITrackingVolumeBuilder
{
public:
  /// @brief This struct stores the data for the construction of a single
  /// PlaneSurface
  struct SurfaceConfig
  {
    // Center position
    Vector3D position;
    // Rotation
    RotationMatrix3D rotation = RotationMatrix3D::Identity();
    // Bounds
    std::shared_ptr<const RectangleBounds> rBounds = nullptr;
    // Attached material
    std::shared_ptr<const SurfaceMaterial> surMat = nullptr;
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
  struct LayerConfig
  {
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
  struct VolumeConfig
  {
    // Center position
    Vector3D position;
    // Lengths in x,y,z
    Vector3D length;
    // Configurations of its layers
    std::vector<LayerConfig> layerCfg;
    // Stored layers
    std::vector<std::shared_ptr<const Layer>> layers;
    // Name of the volume
    std::string name = "Volume";
    // Material
    std::shared_ptr<const Material> material = nullptr;
  };

  /// @brief This struct stores the configuration of the tracking geometry
  struct Config
  {
    // Center position
    Vector3D position;
    // Length in x,y,z
    Vector3D length;
    // Configuration of its volumes
    std::vector<VolumeConfig> volumeCfg;
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
  void
  setConfig(Config& cfg)
  {
    m_cfg = cfg;
  }

  /// @brief This function creates a surface with a given configuration. A
  /// detector element is attached if the template parameter is non-void.
  ///
  /// @param [in] cfg Configuration of the surface
  /// @return Pointer to the created surface
  std::shared_ptr<const PlaneSurface>
  buildSurface(const SurfaceConfig& cfg) const;

  /// @brief This function creates a layer with a surface encaspulated with a
  /// given configuration. The surface gets a detector element attached if the
  /// template parameter is non-void.
  ///
  /// @param [in, out] cfg Configuration of the layer and the surface
  /// @return Pointer to the created layer
  std::shared_ptr<const Layer>
  buildLayer(LayerConfig& cfg) const;

  /// @brief This function creates a TrackingVolume with a configurable number
  /// of layers and surfaces. Each surface gets a detector element attached if
  /// the template parameter is non-void.
  ///
  /// @param [in, out] cfg Configuration of the TrackingVolume
  /// @return Pointer to the created TrackingVolume
  std::shared_ptr<TrackingVolume>
  buildVolume(VolumeConfig& cfg) const;

  /// @brief This function evaluates the minimum and maximum of the binning as
  /// given by the configurations of the surfaces and layers. The ordering
  /// depends on the binning value specified in the configuration of the volume.
  ///
  /// @param [in] cfg Container with the given surfaces and layers
  /// @return Pair containing the minimum and maximum along the binning
  /// direction
  std::pair<double, double>
  binningRange(const VolumeConfig& cfg) const;

  /// @brief This function builds a world TrackingVolume based on a given
  /// configuration
  ///
  /// @return Pointer to the created TrackingGeometry
  std::shared_ptr<TrackingVolume>
  trackingVolume(std::shared_ptr<const TrackingVolume> /*unused*/,
                 std::shared_ptr<const VolumeBounds> /*unused*/) const override;

private:
  /// Configuration of the world volume
  Config m_cfg;
};

}  // namespace Acts
