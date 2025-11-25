// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameterHelpers.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <tuple>

namespace Acts {

template <std::size_t N>
std::tuple<double, std::error_code> GainMatrixUpdater::visitMeasurementImpl(
    const InternalTrackState& trackState, const Logger& logger) const {
  constexpr std::size_t kMeasurementSize = N;
  using ProjectedVector = ActsVector<kMeasurementSize>;
  using ProjectedMatrix = ActsSquareMatrix<kMeasurementSize>;

  typename TrackStateTraits<kMeasurementSize, true>::Parameters
      predictedParameters{trackState.predicted};
  typename TrackStateTraits<kMeasurementSize, true>::Covariance
      predictedCovariance{trackState.predictedCovariance};
  typename TrackStateTraits<kMeasurementSize, false>::Parameters
      filteredParameters{trackState.filtered};
  typename TrackStateTraits<kMeasurementSize, false>::Covariance
      filteredCovariance{trackState.filteredCovariance};
  typename TrackStateTraits<kMeasurementSize, true>::Calibrated calibrated{
      trackState.calibrated};
  typename TrackStateTraits<kMeasurementSize, true>::CalibratedCovariance
      calibratedCovariance{trackState.calibratedCovariance};

  ACTS_VERBOSE("Measurement dimension: " << kMeasurementSize);
  ACTS_VERBOSE("Calibrated measurement: " << calibrated.transpose());
  ACTS_VERBOSE("Calibrated measurement covariance:\n" << calibratedCovariance);

  const std::span<const std::uint8_t, kMeasurementSize> validSubspaceIndices(
      trackState.projector.begin(),
      trackState.projector.begin() + kMeasurementSize);
  const FixedBoundSubspaceHelper<kMeasurementSize> subspaceHelper(
      validSubspaceIndices);

  // TODO use subspace helper for projection instead
  const auto H = subspaceHelper.projector();

  ACTS_VERBOSE("Measurement projector H:\n" << H);

  const auto projectedPredictedCovariance =
      (H * predictedCovariance * H.transpose()).eval();
  const auto K =
      (predictedCovariance * H.transpose() *
       (projectedPredictedCovariance + calibratedCovariance).inverse())
          .eval();

  ACTS_VERBOSE("Gain Matrix K:\n" << K);

  if (K.hasNaN()) {
    // set to error abort execution
    return {0, KalmanFitterError::UpdateFailed};
  }

  filteredParameters =
      predictedParameters + K * (calibrated - H * predictedParameters);
  // Normalize phi and theta
  filteredParameters = normalizeBoundParameters(filteredParameters);
  filteredCovariance =
      (BoundSquareMatrix::Identity() - K * H) * predictedCovariance;
  ACTS_VERBOSE("Filtered parameters: " << filteredParameters.transpose());
  ACTS_VERBOSE("Filtered covariance:\n" << filteredCovariance);

  const ProjectedVector residual = calibrated - H * filteredParameters;
  ACTS_VERBOSE("Residual: " << residual.transpose());

  const ProjectedMatrix m = calibratedCovariance - projectedPredictedCovariance;
  const double chi2 = (residual.transpose() * m.inverse() * residual).value();
  ACTS_VERBOSE("Chi2: " << chi2);

  return {chi2, {}};
}

// Ensure thet the compiler does not implicitly instantiate the template

#define _EXTERN(N)                                    \
  extern template std::tuple<double, std::error_code> \
  GainMatrixUpdater::visitMeasurementImpl<N>(         \
      const InternalTrackState& trackState, const Logger& logger) const

_EXTERN(1);
_EXTERN(2);
_EXTERN(3);
_EXTERN(4);
_EXTERN(5);
_EXTERN(6);

#undef _EXTERN

}  // namespace Acts
