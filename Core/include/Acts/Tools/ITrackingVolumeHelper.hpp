// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ITrackingVolumeHelper.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
// Geometry module
#include "Acts/Utilities/BinningType.hpp"
// Core module
#include <memory>
#include <string>
#include <vector>
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class Layer;
class TrackingVolume;
class VolumeBounds;
class Material;

using LayerPtr                 = std::shared_ptr<const Layer>;
using TrackingVolumePtr        = std::shared_ptr<const TrackingVolume>;
using MutableTrackingVolumePtr = std::shared_ptr<TrackingVolume>;
using VolumeBoundsPtr          = std::shared_ptr<const VolumeBounds>;

using LayerVector          = std::vector<LayerPtr>;
using TrackingVolumeVector = std::vector<TrackingVolumePtr>;

///  @class ITrackingVolumeHelper
///
/// Interface class ITrackingVolumeHelper tools, it inherits from IAlgTool.
/// The ITrackingVolumeHelper is a tool to pack a set of layers into a volume,
/// or - to wrap several volumes into a container volume.
///
/// @todo add documentation how this is done
///
/// TrackingVolumes only exist as std::shared_ptr
///
class ITrackingVolumeHelper
{
public:
  /// Virtual destructor
  virtual ~ITrackingVolumeHelper() = default;

  /// Create a TrackingVolume* from a set of layers and (optional) parameters
  ///
  /// @param layers vector of static layers confined by the TrackingVolume
  /// if no bounds or HepTransform is given, they define the size
  /// together with the volume enevlope parameters
  /// @param matprop dense material properties for this TrackingVolume
  /// @param volBounds (optional) bounds of this TrackingVolume - ownership
  /// given
  /// @param transform (optional) placement of this TrackingVolume - ownership
  /// given
  /// @param volumeName  volume name to be given
  /// @param btype (optional) BinningType - arbitrary(default) or equidistant
  ///
  /// @return shared pointer to a new TrackingVolume
  virtual MutableTrackingVolumePtr
  createTrackingVolume(const LayerVector&                 layers,
                       std::shared_ptr<const Material>    matprop,
                       VolumeBoundsPtr                    volBounds,
                       std::shared_ptr<const Transform3D> transform = nullptr,
                       const std::string& volumeName = "UndefinedVolume",
                       BinningType        btype = arbitrary) const = 0;

  /// Create a TrackingVolume* from a set of layers and (optional) parameters
  ///
  /// @param layers vector of static layers confined by the TrackingVolume
  /// if no bounds or HepTransform is given, they define the size
  /// together with the volume enevlope parameters
  /// @param matprop dense material properties for this TrackingVolume
  /// @param loc0Min, loc0Max, loc1Min, loc1Max : local position in space,
  /// this TrackingVolume is restricted to Translation only
  /// @param volumeName  volume name to be given
  /// @param btype (optional) BinningType - arbitrary(default) or equidistant
  ///
  /// @return shared pointer to a new TrackingVolume
  virtual MutableTrackingVolumePtr
  createTrackingVolume(const LayerVector&              layers,
                       std::shared_ptr<const Material> matprop,
                       double                          loc0Min,
                       double                          loc0Max,
                       double                          loc1Min,
                       double                          loc1Max,
                       const std::string& volumeName = "UndefinedVolume",
                       BinningType        btype = arbitrary) const = 0;

  /// Create a gap volume from dimensions and
  ///
  /// @param matprop dense material properties for this TrackingVolume
  /// @param loc0Min, loc0Max, loc1Min, loc1Max : local position in space,
  /// this TrackingVolume is restricted to Translation only
  /// @param materialLayers number of material layers (aequidistant binning)
  /// @param cylinder type of layers
  /// @param volumeName  volume name to be given
  ///
  /// @return shared pointer to a new TrackingVolume
  virtual MutableTrackingVolumePtr
  createGapTrackingVolume(std::shared_ptr<const Material> matprop,
                          double                          loc0Min,
                          double                          loc0Max,
                          double                          loc1Min,
                          double                          loc1Max,
                          unsigned int                    materialLayers,
                          bool                            cylinder = true,
                          const std::string&              volumeName
                          = "UndefinedVolume") const = 0;

  /// Create a gap volume from dimensions and
  ///
  /// @param matprop dense material properties for this TrackingVolume
  /// @param loc0Min, loc0Max, loc1Min, loc1Max local position in space,
  /// @param layerPositions custom layer positions
  /// @param cylinder type of layers
  /// @param volumeName  : volume name to be given
  /// @param btype (optional) BinningType - arbitrary(default) or equidistant
  ///
  /// @return shared pointer to a new TrackingVolume
  virtual MutableTrackingVolumePtr
  createGapTrackingVolume(std::shared_ptr<const Material> matprop,
                          double                          loc0Min,
                          double                          loc0Max,
                          double                          loc1Min,
                          double                          loc1Max,
                          const std::vector<double>&      layerPositions,
                          bool                            cylinder = true,
                          const std::string& volumeName = "UndefinedVolume",
                          BinningType        btype = arbitrary) const = 0;

  /// Create a one level higher TrackingVolue
  ///
  /// @param volumes the volumes to be contained
  ///
  /// @return shared pointer to a new TrackingVolume
  virtual MutableTrackingVolumePtr
  createContainerTrackingVolume(const TrackingVolumeVector& volumes) const = 0;
};

}  // namespace