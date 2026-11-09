// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/detail/AllFinite.hpp"

#include <cmath>
#include <cstddef>

namespace Acts {

template <std::size_t N>
Result<void> GainMatrixUpdater::visitMeasurementImpl(
    const GeometryContext& gctx, AnyMutableTrackStateProxy trackState,
    const Logger& logger) const {
  constexpr std::size_t kMeasurementSize = N;
  using ProjectedVector = Vector<kMeasurementSize>;
  using ProjectedMatrix = SquareMatrix<kMeasurementSize>;

  const auto calibrated = trackState.calibrated<kMeasurementSize>();
  const auto calibratedCovariance =
      trackState.calibratedCovariance<kMeasurementSize>();

  ACTS_VERBOSE("Measurement dimension: " << kMeasurementSize);
  ACTS_VERBOSE("Calibrated measurement: " << calibrated.transpose());
  ACTS_VERBOSE("Calibrated measurement covariance:\n" << calibratedCovariance);

  const auto validSubspaceIndices =
      trackState.template projectorSubspaceIndices<kMeasurementSize>();

  const FixedBoundSubspaceHelper<kMeasurementSize> subspaceHelper(
      validSubspaceIndices);

  // TODO use subspace helper for projection instead
  const auto H = subspaceHelper.projector();

  ACTS_VERBOSE("Measurement projector H:\n" << H);

  const auto predicted = trackState.predicted();
  const auto predictedCovariance = trackState.predictedCovariance();
  const auto parametrization =
      trackState.hasReferenceSurface()
          ? trackState.referenceSurface().boundParametrization(gctx)
          : BoundParametrization::globalSpherical();

  // The facade interprets periodic direction residuals in the measured
  // subspace, as if the measurement were embedded in the predicted state.
  const ProjectedVector innovation = parametrization.subspaceDifference(
      calibrated, predicted, validSubspaceIndices);
  const ProjectedMatrix innovationInverse =
      (H * predictedCovariance * H.transpose() + calibratedCovariance)
          .inverse();
  const auto K =
      (predictedCovariance * H.transpose() * innovationInverse).eval();

  ACTS_VERBOSE("Gain Matrix K:\n" << K);

  if (!detail::allFinite(K)) {
    // set to error abort execution
    return Result<void>::failure(KalmanFitterError::UpdateFailed);
  }

  BoundVector updated = predicted + K * innovation;
  const auto tmp = (BoundMatrix::Identity() - K * H).eval();
  BoundMatrix updatedCovariance;
  if (!m_useJosephFormulation) {
    updatedCovariance = tmp * predictedCovariance;
  } else {
    updatedCovariance = tmp * predictedCovariance * tmp.transpose() +
                        K * calibratedCovariance * K.transpose();
  }
  if (!detail::allFinite(updated) || !detail::allFinite(updatedCovariance)) {
    return Result<void>::failure(KalmanFitterError::UpdateFailed);
  }
  const auto canonicalization =
      parametrization.canonicalizeInPlace(updated, updatedCovariance);
  if (!canonicalization.ok()) {
    return canonicalization.error();
  }

  // Chi2 uses the innovation in the measurement's original coordinate basis.
  // Folding the output coordinates must not change this statistic.
  const double chi2 =
      (innovation.transpose() * innovationInverse * innovation).value();
  if (!std::isfinite(chi2)) {
    return Result<void>::failure(KalmanFitterError::UpdateFailed);
  }
  trackState.filtered() = updated;
  trackState.filteredCovariance() = updatedCovariance;
  trackState.chi2() = chi2;
  ACTS_VERBOSE("Filtered parameters: " << trackState.filtered().transpose());
  ACTS_VERBOSE("Filtered covariance:\n" << trackState.filteredCovariance());
  ACTS_VERBOSE("Chi2: " << chi2);

  return Result<void>::success();
}

// Ensure that the compiler does not implicitly instantiate the template

#define _EXTERN(N)                                                         \
  extern template Result<void> GainMatrixUpdater::visitMeasurementImpl<N>( \
      const GeometryContext& gctx, AnyMutableTrackStateProxy trackState,   \
      const Logger& logger) const

_EXTERN(1);
_EXTERN(2);
_EXTERN(3);
_EXTERN(4);
_EXTERN(5);
_EXTERN(6);

#undef _EXTERN

}  // namespace Acts
