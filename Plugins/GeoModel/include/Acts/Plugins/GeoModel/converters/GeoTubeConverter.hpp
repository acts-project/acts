// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Plugins/GeoModel/detail/GenericGeoShapeConverter.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <tuple>

#include <GeoModelKernel/GeoTube.h>

class GeoFullPhysVol;
class GeoTube;

namespace Acts {

namespace detail {
struct GeoTubeConverter {
  Surface::SurfaceType targetShape = Surface::SurfaceType::Straw;

  /// @brief Convert a GeoTube to a detector element and surface
  ///
  /// @param geoFPV The full physical volume to convert (contains shape)
  /// @param geoTube The GeoTube shape to convert
  /// @param absTransform from the GeoPhysVol
  /// @param bool sensitive
  ///
  /// @return The detector element and surface
  Result<GeoModelSensitiveSurface> operator()(const GeoFullPhysVol& geoFPV,
                                              const GeoTube& geoTube,
                                              const Transform3& absTransform,
                                              bool sensitive) const;
};
}  // namespace detail

/// @brief The GeoTube converter
///
/// This is a dedicated converter for GeoTube shapes
using GeoTubeConverter =
    detail::GenericGeoShapeConverter<GeoTube, detail::GeoTubeConverter>;

}  // namespace Acts
