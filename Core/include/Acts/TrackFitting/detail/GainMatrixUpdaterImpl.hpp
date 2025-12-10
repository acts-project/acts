// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/AnyTrackState.hpp"
#include "Acts/EventData/TrackParameterHelpers.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <tuple>

namespace Acts {

template <std::size_t N>
std::tuple<double, std::error_code> GainMatrixUpdater::visitMeasurementImpl(
    AnyMutableTrackState trackState, const Logger& logger) const {
  double chi2 = 0;

  constexpr std::size_t kMeasurementSize = N;
  using ParametersVector = ActsVector<kMeasurementSize>;
  using CovarianceMatrix = ActsSquareMatrix<kMeasurementSize>;

  auto calibrated = trackState.calibrated<kMeasurementSize>();
  auto calibratedCovariance =
      trackState.calibratedCovariance<kMeasurementSize>();

  ACTS_VERBOSE("Measurement dimension: " << kMeasurementSize);
  ACTS_VERBOSE("Calibrated measurement: " << calibrated.transpose());
  ACTS_VERBOSE("Calibrated measurement covariance:\n" << calibratedCovariance);

  const auto validSubspaceIndices =
      trackState.template projectorSubspaceIndices<kMeasurementSize>();

  FixedBoundSubspaceHelper<kMeasurementSize> subspaceHelper(
      validSubspaceIndices);

  // TODO use subspace helper for projection instead
  const auto H = subspaceHelper.projector();

  ACTS_VERBOSE("Measurement projector H:\n" << H);

  auto filtered = trackState.filtered();
  auto filteredCovariance = trackState.filteredCovariance();
  const auto predicted = trackState.predicted();
  const auto predictedCovariance = trackState.predictedCovariance();

  const auto K =
      (predictedCovariance * H.transpose() *
       (H * predictedCovariance * H.transpose() + calibratedCovariance)
           .inverse())
          .eval();

  ACTS_VERBOSE("Gain Matrix K:\n" << K);

  if (K.hasNaN()) {
    // set to error abort execution
    return {0, KalmanFitterError::UpdateFailed};
  }

  filtered = predicted + K * (calibrated - H * predicted);
  // Normalize phi and theta
  filtered = normalizeBoundParameters(filtered);
  filteredCovariance =
      (BoundSquareMatrix::Identity() - K * H) * predictedCovariance;
  ACTS_VERBOSE("Filtered parameters: " << filtered.transpose());
  ACTS_VERBOSE("Filtered covariance:\n" << filteredCovariance);

  ParametersVector residual;
  residual = calibrated - H * filtered;
  ACTS_VERBOSE("Residual: " << residual.transpose());

  CovarianceMatrix m =
      ((CovarianceMatrix::Identity() - H * K) * calibratedCovariance).eval();

  chi2 = (residual.transpose() * m.inverse() * residual).value();

  ACTS_VERBOSE("Chi2: " << chi2);

  return {chi2, {}};
}

// Ensure thet the compiler does not implicitly instantiate the template

#define _EXTERN(N)                                                            \
  extern template std::tuple<double, std::error_code>                         \
  GainMatrixUpdater::visitMeasurementImpl<N>(AnyMutableTrackState trackState, \
                                             const Logger& logger) const

_EXTERN(1);
_EXTERN(2);
_EXTERN(3);
_EXTERN(4);
_EXTERN(5);
_EXTERN(6);

#undef _EXTERN

}  // namespace Acts
