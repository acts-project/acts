// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <map>

#include <detray/core/detector.hpp>
#include <detray/definitions/grid_axis.hpp>
#include <detray/detectors/default_metadata.hpp>
#include <detray/io/frontend/payloads.hpp>

namespace Acts::Experimental {
class DetectorVolume;
}  // namespace Acts::Experimental

namespace ActsPlugins {

using DetrayMetaData = detray::default_metadata<detray::array<double>>;

using DetrayHostDetector = detray::detector<DetrayMetaData>;

namespace DetrayConversionUtils {

/// Detray conversion cache object
///
/// This object is used to synchronize link information between the
/// different converters (geometry, material, surface grids)
struct Cache {
  /// Explicit constructor with detector volumes
  ///
  /// @param detectorVolumes the number of detector volumes
  explicit Cache(
      const std::vector<const Acts::Experimental::DetectorVolume*>& dVolumes)
      : detectorVolumes(dVolumes) {}

  /// The volumes of the detector for index lookup
  std::vector<const Acts::Experimental::DetectorVolume*> detectorVolumes;
  /// This is a map to pass on volume link information
  std::map<Acts::GeometryIdentifier, unsigned long> volumeLinks;
  /// This is a multimap to pass volume local surface link information
  /// The portal splitting requires a multimap implementation here
  ///
  /// These are volume local, hence indexed per volumes
  std::map<std::size_t, std::multimap<Acts::GeometryIdentifier, unsigned long>>
      localSurfaceLinks;

  /// Find the position of the volume to point to
  ///
  /// @param volume the volume to find
  ///
  /// @note throws exception if volume is not found
  std::size_t volumeIndex(
      const Acts::Experimental::DetectorVolume* volume) const {
    if (auto candidate = std::ranges::find(detectorVolumes, volume);
        candidate != detectorVolumes.end()) {
      return std::distance(detectorVolumes.begin(), candidate);
    }
    throw std::invalid_argument("Volume not found in the cache");
  }
};

/// Convert the binning option
///
/// @param bOption the binning option
///
/// @return a detray binning option
detray::axis::bounds convertBinningOption(Acts::BinningOption bOption);

/// Convert the binning value
///
/// @param bValue the binning value
///
/// @return a detray binning value
detray::axis::label convertAxisDirection(Acts::AxisDirection bValue);

/// Convert the binning type
///
/// @param bType the binning type
///
/// @return a detray binning type
detray::axis::binning convertBinningType(Acts::BinningType bType);

/// Convert the binning data to an axis
///
/// @param bData the binning data to be converted
///
/// @return a detray axis payload
detray::io::axis_payload convertBinningData(const Acts::BinningData& bData);

}  // namespace DetrayConversionUtils
}  // namespace ActsPlugins
