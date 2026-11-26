// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/BoundParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cassert>

namespace Acts {

/// Kalman update step using the gain matrix formalism.
/// Uses the supplied bound parametrization of the track state's reference
/// surface for coordinate residuals and canonicalizes parameters and covariance
/// together. States without a reference surface use the canonical
/// global-spherical interpretation.
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
  /// @param[in] gctx Geometry context, used to check @p parametrization in
  ///            debug builds
  /// @param[in] parametrization Bound parametrization of the track state's
  ///            reference surface in @p gctx, or the canonical global-spherical
  ///            one for a state without a reference surface
  /// @param[in,out] trackState The track state
  /// @param[in] logger Where to write logging information to
  /// @return Success or failure of the update procedure
  template <typename traj_t>
  Result<void> operator()([[maybe_unused]] const GeometryContext& gctx,
                          const BoundParametrization& parametrization,
                          typename traj_t::TrackStateProxy trackState,
                          const Logger& logger = getDummyLogger()) const {
    ACTS_VERBOSE("Invoked GainMatrixUpdater");

    // the caller resolves the parametrization of the state's surface
    assert(parametrization ==
           (trackState.hasReferenceSurface()
                ? trackState.referenceSurface().boundParametrization(gctx)
                : BoundParametrization::globalSpherical()));

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

    return visitMeasurement(parametrization,
                            AnyMutableTrackStateProxy{trackState}, logger);
  }

 private:
  bool m_useJosephFormulation = false;

  Result<void> visitMeasurement(const BoundParametrization& parametrization,
                                AnyMutableTrackStateProxy trackState,
                                const Logger& logger) const;

  template <std::size_t N>
  Result<void> visitMeasurementImpl(const BoundParametrization& parametrization,
                                    AnyMutableTrackStateProxy trackState,
                                    const Logger& logger) const;
};

}  // namespace Acts
