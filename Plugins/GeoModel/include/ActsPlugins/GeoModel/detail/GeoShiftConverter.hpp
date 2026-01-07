// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/BoundFactory.hpp"
#include "ActsPlugins/GeoModel/GeoModelDetectorElement.hpp"

#include <memory>
#include <tuple>

#include <GeoModelKernel/GeoShapeShift.h>

namespace ActsPlugins::detail {

struct GeoShiftConverter {
  /// @brief Convert a GeoBox to a detector element and surface
  ///
  /// @param geoFPV The full physical volume to convert (contains shape)
  /// @param geoBox The GeoBox to convert
  /// @param absTransform from the GeoPhysVol
  /// @param bool sensitive
  ///
  /// @return The detector element and surface
  Acts::Result<GeoModelSensitiveSurface> operator()(
      const PVConstLink& geoPV, const GeoShapeShift& geoShift,
      const Acts::Transform3& absTransform,
      Acts::SurfaceBoundFactory& boundFactory, bool sensitive) const;
};

}  // namespace ActsPlugins::detail
