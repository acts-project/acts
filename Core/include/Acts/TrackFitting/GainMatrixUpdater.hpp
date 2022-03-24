// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

/// Kalman update step using the gain matrix formalism.
class GainMatrixUpdater {
 public:
  /// Run the Kalman update step for a single trajectory state.
  ///
  /// @tparam kMeasurementSizeMax
  /// @param[in] gctx The current geometry context object, e.g. alignment
  /// @param[in,out] trackState The track state
  /// @param[in] direction The navigation direction
  /// @param[in] logger Where to write logging information to
  Result<void> operator()(
      const GeometryContext& gctx, MultiTrajectory::TrackStateProxy trackState,
      NavigationDirection direction = NavigationDirection::forward,
      LoggerWrapper logger = getDummyLogger()) const;
};

}  // namespace Acts
