// This file is part of the Acts project.
//
// Copyright (C) 2016-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

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
    unsigned int calibratedSize;
    // This is used to build a covariance matrix view in the .cpp file
    const double* calibrated;
    const double* calibratedCovariance;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Projector projector;

    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Parameters predicted;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Covariance predictedCovariance;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Parameters filtered;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Covariance filteredCovariance;
  };

 public:
  /// Run the Kalman update step for a single trajectory state.
  ///
  /// @tparam kMeasurementSizeMax
  /// @param[in,out] trackState The track state
  /// @param[in] logger Where to write logging information to
  template <typename traj_t>
  Result<void> operator()(const GeometryContext& /*gctx*/,
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

    auto [chi2, error] = visitMeasurement(
        InternalTrackState{
            trackState.calibratedSize(),
            // Note that we pass raw pointers here which are used in the correct
            // shape later
            trackState.effectiveCalibrated().data(),
            trackState.effectiveCalibratedCovariance().data(),
            trackState.projector(),
            trackState.predicted(),
            trackState.predictedCovariance(),
            trackState.filtered(),
            trackState.filteredCovariance(),
        },
        logger);

    trackState.chi2() = chi2;

    return error ? Result<void>::failure(error) : Result<void>::success();
  }

 private:
  std::tuple<double, std::error_code> visitMeasurement(
      InternalTrackState trackState, const Logger& logger) const;

  template <std::size_t N>
  std::tuple<double, std::error_code> visitMeasurementImpl(
      InternalTrackState trackState, const Logger& logger) const;
};

template <std::size_t N>
std::tuple<double, std::error_code> GainMatrixUpdater::visitMeasurementImpl(
    InternalTrackState trackState, const Logger& logger) const {
  double chi2 = 0;

  constexpr std::size_t kMeasurementSize = N;
  using ParametersVector = ActsVector<kMeasurementSize>;
  using CovarianceMatrix = ActsSquareMatrix<kMeasurementSize>;

  typename TrackStateTraits<kMeasurementSize, true>::Calibrated calibrated{
      trackState.calibrated};
  typename TrackStateTraits<kMeasurementSize, true>::CalibratedCovariance
      calibratedCovariance{trackState.calibratedCovariance};

  ACTS_VERBOSE("Measurement dimension: " << kMeasurementSize);
  ACTS_VERBOSE("Calibrated measurement: " << calibrated.transpose());
  ACTS_VERBOSE("Calibrated measurement covariance:\n" << calibratedCovariance);

  const auto H = trackState.projector
                     .template topLeftCorner<kMeasurementSize, eBoundSize>()
                     .eval();

  ACTS_VERBOSE("Measurement projector H:\n" << H);

  const auto K = (trackState.predictedCovariance * H.transpose() *
                  (H * trackState.predictedCovariance * H.transpose() +
                   calibratedCovariance)
                      .inverse())
                     .eval();

  ACTS_VERBOSE("Gain Matrix K:\n" << K);

  if (K.hasNaN()) {
    // set to error abort execution
    return {0, KalmanFitterError::UpdateFailed};
  }

  trackState.filtered =
      trackState.predicted + K * (calibrated - H * trackState.predicted);
  trackState.filteredCovariance =
      (BoundSquareMatrix::Identity() - K * H) * trackState.predictedCovariance;
  ACTS_VERBOSE("Filtered parameters: " << trackState.filtered.transpose());
  ACTS_VERBOSE("Filtered covariance:\n" << trackState.filteredCovariance);

  ParametersVector residual;
  residual = calibrated - H * trackState.filtered;
  ACTS_VERBOSE("Residual: " << residual.transpose());

  CovarianceMatrix m =
      ((CovarianceMatrix::Identity() - H * K) * calibratedCovariance).eval();

  chi2 = (residual.transpose() * m.inverse() * residual).value();

  ACTS_VERBOSE("Chi2: " << chi2);

  return {chi2, {}};
}

#define _EXTERN(N)                                                          \
  extern template std::tuple<double, std::error_code>                       \
  GainMatrixUpdater::visitMeasurementImpl<N>(InternalTrackState trackState, \
                                             const Logger& logger) const

_EXTERN(1);
_EXTERN(2);
_EXTERN(3);
_EXTERN(4);
_EXTERN(5);
_EXTERN(6);

#undef _EXTERN

}  // namespace Acts
