// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/CompositeSpacePoint.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
namespace Acts::Experimental {

/// @brief Interface concept to calibrate CompositeSpacePoints

template <typename Calibrator_t, typename SpacePoint_t>
concept CompositeSpacePointFastCalibrator =
    CompositeSpacePoint<SpacePoint_t> &&
    requires(const Calibrator_t calibrator, const Acts::CalibrationContext& ctx,
             const SpacePoint_t& spacePoint, const double t0) {
      /// @brief Returns the drift velocity of the straw measurement's radius - time relation
      ///        which is defined as the first derivative of the relation.
      /// @param ctx: Calibration context to access the calibration constants (Experiment specific)
      /// @param spacePoint: Reference to the calibrated space point
      { calibrator.driftVelocity(ctx, spacePoint, t0) } -> std::same_as<double>;
      /// @brief Returns the drift acceleration of the straw measurement's radius - time relation
      ///        which is defined as the second derivative of the relation
      /// @param ctx: Calibration context to access the calibration constants (Experiment specific)
      /// @param spacePoint: Reference to the calibrated space point
      {
        calibrator.driftAcceleration(ctx, spacePoint, t0)
      } -> std::same_as<double>;

      { calibrator.driftRadius(ctx, spacePoint, t0) } -> std::same_as<double>;
    };
}  // namespace Acts::Experimental
