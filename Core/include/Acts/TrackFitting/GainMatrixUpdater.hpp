// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cassert>
#include <system_error>
#include <tuple>

namespace Acts {

/// Kalman update step using the gain matrix formalism.
class GainMatrixUpdater {
  struct InternalTrackState {
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Parameters predicted;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Covariance predictedCovariance;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Parameters filtered;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Covariance filteredCovariance;
    // This is used to build a covariance matrix view in the .cpp file
    double* calibrated;
    double* calibratedCovariance;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Projector projector;
    unsigned int calibratedSize;
  };

 public:
  /// Run the Kalman update step for a single trajectory state.
  ///
  /// @tparam kMeasurementSizeMax
  /// @param[in,out] trackState The track state
  /// @param[in] direction The navigation direction
  /// @param[in] logger Where to write logging information to
  template <typename traj_t>
  Result<void> operator()(const GeometryContext& /*gctx*/,
                          typename traj_t::TrackStateProxy trackState,
                          Direction direction = Direction::Forward,
                          const Logger& logger = getDummyLogger()) const {
    ACTS_VERBOSE("Invoked GainMatrixUpdater");

    // there should be a calibrated measurement
    assert(trackState.hasCalibrated());
    // we should have predicted state set
    assert(trackState.hasPredicted());
    // filtering should not have happened yet, but is allocated, therefore set
    assert(trackState.hasFiltered());

    // read-only handles. Types are eigen maps to backing storage
    // const auto predicted = trackState.predicted();
    // const auto predictedCovariance = trackState.predictedCovariance();

    ACTS_VERBOSE(
        "Predicted parameters: " << trackState.predicted().transpose());
    ACTS_VERBOSE("Predicted covariance:\n" << trackState.predictedCovariance());

    // read-write handles. Types are eigen maps into backing storage.
    // This writes directly into the trajectory storage
    // auto filtered = trackState.filtered();
    // auto filteredCovariance = trackState.filteredCovariance();

    auto [chi2, error] = visitMeasurement(
        InternalTrackState{
            trackState.predicted(),
            trackState.predictedCovariance(),
            trackState.filtered(),
            trackState.filteredCovariance(),
            // This abuses an incorrectly sized vector / matrix to access the
            // data pointer! This works (don't use the matrix as is!), but be
            // careful!
            trackState
                .template calibrated<
                    MultiTrajectoryTraits::MeasurementSizeMax>()
                .data(),
            trackState
                .template calibratedCovariance<
                    MultiTrajectoryTraits::MeasurementSizeMax>()
                .data(),
            trackState.projector(),
            trackState.calibratedSize(),
        },
        direction, logger);

    trackState.chi2() = chi2;

    return error ? Result<void>::failure(error) : Result<void>::success();
  }

 private:
  std::tuple<double, std::error_code> visitMeasurement(
      InternalTrackState trackState, Direction direction,
      const Logger& logger) const;
};

}  // namespace Acts
