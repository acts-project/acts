// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"

#include <string>

namespace Acts {

/// A proto volume description being used to define an overall
/// structure of either a TrackingVolume or Experimental::DetectorVolume
struct ProtoVolume {
  /// Name of the proto volume
  std::string name = "";
  /// The extent of this layer
  Extent extent;
  /// Boolean to indicate this is a container, needed for legacy
  bool layerContainer = false;

  /// Set the layer type to indicate the layer volume
  Acts::Surface::SurfaceType layerType = Acts::Surface::SurfaceType::Other;
  /// The surface binninng fo the layer
  std::vector<BinningData> layerSurfaceBinning = {};
  /// Internal structure container
  std::vector<ProtoVolume> constituentVolumes = {};
  /// The constituent binning if this a container
  std::vector<BinningData> constituentBinning = {};

  /// Define an operator==
  ///
  /// @param ptVolume the proto volume to be checked
  bool operator==(const ProtoVolume& ptVolume) const;

  /// Harmonize the detector information, this can run in two
  /// modes, steered by the @param legacy boolean
  ///
  /// The legacy mode prepares everything for `Acts::TrackingVolume`,
  /// if off it creates a description for `Acts::Detector`.
  void harmonize(bool legacy = true);

  /// Extend the tracking volume with the its own constituents,
  /// upwards here means that extens are promoted to the mother
  ///
  /// @param ptVolume the protoVolume
  void extendUp(ProtoVolume& ptVolume);

  /// Extend the tracking volume with the its own constituents
  /// @param bValue the binning value that is propagated
  void propagateMinDown(BinningValue bValue);

  /// Extend the tracking volume with the its own constituents
  /// @param bValue the binning value that is propagated
  void propagateMaxDown(BinningValue bValue);

  /// Constrain the daughter volumes with this volume
  ///
  /// @param ptVolume is the proto volume from which the constrain
  /// is taken
  void constrainDown(const ProtoVolume& ptVolume);

  /// Write the tracking volume to screen
  /// @param indent the current indentation
  std::string toString(const std::string& indent = "") const;
};

/// A proto detector description being used to define an overall
/// structure of either a TrackingGeometry or Experimental::Detector
struct ProtoDetector {
  std::string name = "";
  ProtoVolume worldVolume;

  /// Harmonize the detector information, this can run in two
  /// modes, steered by the @param legacy boolean
  ///
  /// The legacy mode prepares everything for `Acts::TrackingVolume`,
  /// if off it creates a description for `Acts::Detector`.
  ///
  void harmonize(bool legacy = true) {
    worldVolume.extendUp(worldVolume);
    worldVolume.constrainDown(worldVolume);
    worldVolume.harmonize(legacy);
  }

  /// Define an operator==
  ///
  /// @param pd the proto detector to be checked
  bool operator==(const ProtoDetector& pd) const {
    return (name == pd.name and worldVolume == pd.worldVolume);
  }

  /// Write the tracking volume to screen
  /// @param indent the current indentation
  std::string toString(const std::string& indent = "") const;
};

}  // namespace Acts
