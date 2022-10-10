// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Experimental/detail/DetectorVolumeLinks.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <memory>
#include <tuple>
#include <vector>

namespace Acts {
namespace Experimental {
namespace detail {

struct CylindricalDetectorHelperOptions {
  /// @brief  Handling strategy
  enum class Handling { eStrict, eResize, eFill };
  /// @brief  Specify a handling strategy
  Handling handling = Handling::eStrict;
  /// @brief  Tolerance at which glueing is still attempted
  ActsScalar glueTolerance = std::numeric_limits<ActsScalar>::epsilon();
  /// ACTS log level
  Logging::Level logLevel = Logging::INFO;
};

/// @brief Method to attach two cylindrical detector volumes, i.e. to
/// create a shared, pointing boundary between the volumes
///
/// @param gctx the geometry context of this call
/// @param volumes the are the detector volumes
/// @param options is a boolean to steer whether the volumes should be resized
///
/// @note throws exceptions if misconfigured
/// @note currently relies on sorted volumes
///
void connectCylindricalVolumes(
    const GeometryContext& gctx,
    const std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    const CylindricalDetectorHelperOptions& options) noexcept(false);

}
}  // namespace Experimental
}  // namespace Acts