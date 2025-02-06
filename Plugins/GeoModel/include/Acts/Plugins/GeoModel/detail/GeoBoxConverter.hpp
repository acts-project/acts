// SPDX-PackageName: "ACTS"
// SPDX-FileCopyrightText: 2016 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <tuple>

#include <GeoModelKernel/GeoBox.h>

namespace Acts::detail {

struct GeoBoxConverter {
  /// @brief Convert a GeoBox to a detector element and surface
  ///
  /// @param geoFPV The full physical volume to convert (contains shape)
  /// @param geoBox The GeoBox to convert
  /// @param absTransform from the GeoPhysVol
  /// @param bool sensitive
  ///
  /// @return The detector element and surface
  Result<GeoModelSensitiveSurface> operator()(const PVConstLink& geoPV,
                                              const GeoBox& geoBox,
                                              const Transform3& absTransform,
                                              bool sensitive) const;
};

}  // namespace Acts::detail
