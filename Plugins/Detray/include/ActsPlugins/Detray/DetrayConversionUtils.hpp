// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <map>
#include <tuple>

#include <detray/core/detector.hpp>
#include <detray/definitions/grid_axis.hpp>
#include <detray/detectors/default_metadata.hpp>
#include <detray/io/frontend/payloads.hpp>

namespace ActsPlugins {

/// Detray metadata type
using DetrayMetaData = detray::default_metadata<detray::array<double>>;

/// Detray host detector type
using DetrayHostDetector = detray::detector<DetrayMetaData>;

/// @ingroup detray_plugin
namespace DetrayConversionUtils {

/// @addtogroup detray_plugin
/// @{

/// Detray conversion cache object
///
/// This object is used to synchronize link information between the
/// different converters (geometry, material, surface grids)
struct Cache {
  /// This is a map to pass on volume link information
  std::map<Acts::GeometryIdentifier, unsigned long> volumeLinks;
  /// This is a multimap to pass volume local surface link information
  /// The portal splitting requires a multimap implementation here
  ///
  /// These are volume local, hence indexed per volumes
  std::map<std::size_t, std::multimap<Acts::GeometryIdentifier, unsigned long>>
      localSurfaceLinks;
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

/// Convert an IAxis to a detray axis payload
///
/// @param axis the axis to be converted
///
/// @return a detray axis payload
detray::io::axis_payload convertAxis(const Acts::IAxis& axis);

/// Convert a MaterialSlab to a detray material slab payload
///
/// @param slab the material slab to be converted
///
/// @return a detray material slab payload
detray::io::material_slab_payload convertMaterialSlab(
    const Acts::MaterialSlab& slab);

/// Convert a Transform3 to a detray transform payload
///
/// @param transform the transform to be converted
///
/// @return a detray transform payload
detray::io::transform_payload convertTransform(
    const Acts::Transform3& transform);

/// Convert a 1D BinUtility to a 2D BinUtility for Detray
///
/// Detray expects 2-dimensional grids. This function converts 1D grids
/// to 2D by adding a dummy second dimension. Currently supported 2D grids
/// are: x-y, r-phi, phi-z
///
/// @param bUtility the bin utility to be converted (may be 1D or 2D)
///
/// @return a tuple containing:
///   - the converted 2D BinUtility
///   - a boolean indicating if axes were swapped
std::tuple<Acts::BinUtility, bool> convertBinUtilityTo2D(
    const Acts::BinUtility& bUtility);

/// @}

}  // namespace DetrayConversionUtils
}  // namespace ActsPlugins
