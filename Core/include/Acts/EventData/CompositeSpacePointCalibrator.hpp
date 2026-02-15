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

/// @brief Interface concept to define the straw measurement calibrator used by the FastStrawLineFitter.
///         The basic assumption is that the straw radii are reconstructed from
///         a time measurement which is converted into a drift-radius using a
///         differentiable r-t relation. Due to e.g. late arrival of the
///         particle, the drift time may be subject to corrections which is
///         summarized by a general time shift t0. This shift can be estimated
///         during a straight line fit. The calibrator returns the updated drift
///         radius & the first and second derivative of the r-t relation called
///         driftVelocity and driftAcceleration, respectively.
template <typename Calibrator_t, typename SpacePoint_t>
concept CompositeSpacePointFastCalibrator =
    CompositeSpacePoint<SpacePoint_t> &&
    requires(const Calibrator_t calibrator, const Acts::CalibrationContext& ctx,
             const SpacePoint_t& spacePoint, const double t0) {
      /// @brief Returns the drift velocity of the straw measurement's radius - time relation
      ///        which is defined as the first derivative of the relation.
      /// @param ctx: Calibration context to access the calibration constants (Experiment specific)
      /// @param spacePoint: Reference to the calibrated spacepoint
      { calibrator.driftVelocity(ctx, spacePoint, t0) } -> std::same_as<double>;
      /// @brief Returns the drift acceleration of the straw measurement's radius - time relation
      ///        which is defined as the second derivative of the relation
      /// @param ctx: Calibration context to access the calibration constants (Experiment specific)
      /// @param spacePoint: Reference to the calibrated spacepoint
      {
        calibrator.driftAcceleration(ctx, spacePoint, t0)
      } -> std::same_as<double>;

      { calibrator.driftRadius(ctx, spacePoint, t0) } -> std::same_as<double>;
    };

/// @brief Interface concept for a CompositeSpacePointCalibrator. The spacepoint
///        container is parsed to the interface together with a reference track
///        position, direction and the time offset, all expressed at the
///        reference plane z=0.
template <typename Calibrator_t, typename UnCalibCont_t, typename CalibCont_t>
concept CompositeSpacePointCalibrator =
    CompositeSpacePointContainer<UnCalibCont_t> &&
    CompositeSpacePointContainer<CalibCont_t> &&
    requires(const Calibrator_t calibrator, const UnCalibCont_t& uncalibCont,
             CalibCont_t& calibCont, const Vector3& trackPos,
             const Vector3& trackDir, const double trackT0,
             const CalibrationContext& ctx,
             const Acts::RemovePointer_t<typename CalibCont_t::value_type>&
                 measurement) {
      ///  @brief Calibrate the entire input spacepoint container using the external track parameters
      ///  @param ctx: Calibration context to access the calibration constants (Experiment specific)
      ///  @param trackPos: Position of the track / segment
      ///  @param trackDir: Direction of the track / segment
      ///  @param trackT0: Time offset w.r.t. the nominal time of flight calculated by (globPos) / c
      ///  @param uncalibCont: Const reference to the calibrated input container to calibrate
      {
        calibrator.calibrate(ctx, trackPos, trackDir, trackT0, uncalibCont)
      } -> std::same_as<CalibCont_t>;
      /// @brief Update the signs of each straw measurement according to whether the measurement is on the
      ///        left or on the right-hand side of the given line
      ///  @param trackPos: Position of the track / segment
      ///  @param trackDir: Direction of the track / segment
      ///  @param calibCont: Mutable reference to the calibrate composite spacepoint container
      {
        calibrator.updateSigns(trackPos, trackDir, calibCont)
      } -> std::same_as<void>;
      /// @brief Returns the drift velocity of the straw measurement's radius - time relation
      ///        which is defined as the first derivative of the relation.
      ///  @param ctx: Calibration context to access the calibration constants (Experiment specific)
      ///  @param measurement: Reference to the calibrated spacepoint
      { calibrator.driftVelocity(ctx, measurement) } -> std::same_as<double>;
      /// @brief Returns the drift acceleration of the straw measurement's radius - time relation
      ///        which is defined as the second derivative of the relation
      ///  @param ctx: Calibration context to access the calibration constants (Experiment specific)
      ///  @param measurement: Reference to the calibrated spacepoint
      {
        calibrator.driftAcceleration(ctx, measurement)
      } -> std::same_as<double>;
    };

}  // namespace Acts::Experimental
