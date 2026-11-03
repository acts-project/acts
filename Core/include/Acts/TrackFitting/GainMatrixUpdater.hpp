// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cassert>

namespace Acts {

/// Kalman update step using the gain matrix formalism.
/// Uses the reference surface's bound parametrization for coordinate residuals
/// and canonicalizes parameters and covariance together. States without a
/// reference surface retain the canonical global-spherical interpretation.
/// A failed update leaves filtered values, covariance and chi2 unchanged.
/// @ingroup track_fitting
class GainMatrixUpdater {
 public:
  GainMatrixUpdater() = default;

  /// @param useJosephFormulation Whether to use the Joseph formulation for the
  /// covariance update, which is more numerically stable at the cost of
  /// additional computations.
  explicit GainMatrixUpdater(bool useJosephFormulation)
      : m_useJosephFormulation(useJosephFormulation) {}

  /// Run the Kalman update step for a single trajectory state.
  ///
  /// @tparam traj_t Trajectory type
  /// @param[in] gctx Geometry context for bound-coordinate interpretation
  /// @param[in,out] trackState The track state
  /// @param[in] logger Where to write logging information to
  /// @return Success or failure of the update procedure
  template <typename traj_t>
  Result<void> operator()(const GeometryContext& gctx,
                          typename traj_t::TrackStateProxy trackState,
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

    return visitMeasurement(gctx, AnyMutableTrackStateProxy{trackState},
                            logger);
  }

 private:
  bool m_useJosephFormulation = false;

  Result<void> visitMeasurement(const GeometryContext& gctx,
                                AnyMutableTrackStateProxy trackState,
                                const Logger& logger) const;

  template <std::size_t N>
  Result<void> visitMeasurementImpl(const GeometryContext& gctx,
                                    AnyMutableTrackStateProxy trackState,
                                    const Logger& logger) const;
};

}  // namespace Acts
