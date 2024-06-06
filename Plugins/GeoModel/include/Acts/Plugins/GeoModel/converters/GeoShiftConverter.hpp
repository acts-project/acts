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
#include "Acts/Plugins/GeoModel/converters/GeoBoxConverter.hpp"
#include "Acts/Plugins/GeoModel/converters/GeoTrdConverter.hpp"
#include "Acts/Plugins/GeoModel/converters/GeoTubeConverter.hpp"
#include "Acts/Plugins/GeoModel/detail/GenericGeoShapeConverter.hpp"
#include "Acts/Utilities/Result.hpp"
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/DiscSurface.hpp>
#include <Acts/Surfaces/TrapezoidBounds.hpp>
#include <Acts/Surfaces/CylinderBounds.hpp>
#include <Acts/Surfaces/StrawSurface.hpp>
#include <Acts/Surfaces/LineBounds.hpp>

#include <memory>
#include <tuple>

#include <GeoModelKernel/GeoShapeShift.h>

namespace Acts {

namespace detail {

struct GeoShiftConverter {
  /// @brief Convert a GeoBox to a detector element and surface
  ///
  /// @param geoFPV The full physical volume to convert (contains shape)
  /// @param geoBox The GeoBox to convert
  /// @param absTransform from the GeoPhysVol
  /// @param bool sensitive
  ///
  /// @return The detector element and surface
  Result<GeoModelSensitiveSurface>
  operator()(const GeoFullPhysVol& geoFPV, const GeoShapeShift& geoShift,
             const Transform3& absTransform, bool sensitive) const;
};
}  // namespace detail

/// @brief The GeoShift + Trd/Box/Tube converter
///
/// This is a dedicated converter for GeoBox shapes
using GeoShiftConverter =
    detail::GenericGeoShapeConverter<GeoShapeShift, detail::GeoShiftConverter>;

}  // namespace Acts
