// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {
namespace detail {
/// Calculate the global track parameters covariance for a smoothed trajectory
/// stored in MultiTrajecty based on formulas at Journal of Physics: Conference
/// Series 219 (2010) 032028.
///
/// @tparam source_link_t The source link type of the trajectory
/// @tparam parameters_t The track parameters type
///
/// @param multiTraj The MultiTrajectory containing the trajectory to be
/// investigated
/// @param entryIndex The trajectory entry index
/// @param covarianceForMeasurementStates Indicator for whether to consider only
/// measurement states
///
/// @return The global track parameters covariance matrix
template <typename source_link_t, typename parameters_t = BoundParameters>
ActsMatrixX<BoundParametersScalar> globalTrackParametersCovariance(
    const Acts::MultiTrajectory<source_link_t>& multiTraj,
    const size_t& entryIndex, bool covarianceForMeasurementStates = true) {
  using CovMatrix_t = typename parameters_t::CovMatrix_t;
  using gain_matrix_t = CovMatrix_t;

  // The last smoothed state index
  size_t lastSmoothedIndex = SIZE_MAX;
  // The total number of smoothed states
  size_t nSmoothedStates = 0;
  // The order of smoothing for those measurement states (necessary to retrieve
  // elements relevant to measurements in the full global track parameters
  // covariance)
  std::vector<size_t> smoothingOrders;
  smoothingOrders.reserve(15);
  // Visit all the states
  multiTraj.visitBackwards(entryIndex, [&](const auto& ts) {
    if (ts.hasSmoothed()) {
      if (lastSmoothedIndex == SIZE_MAX) {
        lastSmoothedIndex = ts.index();
      }
      nSmoothedStates++;
      if (ts.typeFlags().test(TrackStateFlag::MeasurementFlag)) {
        smoothingOrders.push_back(nSmoothedStates);
      }
    }
  });

  // Set the size of global track parameters covariance for all smoothed states
  ActsMatrixX<BoundParametersScalar> fullGlobalTrackParamsCov;
  fullGlobalTrackParamsCov.resize(nSmoothedStates * eBoundParametersSize,
                                  nSmoothedStates * eBoundParametersSize);
  // Visit the smoothed states to calculate the full global track parameters
  // covariance
  size_t nProcessed = 0;
  auto prev_ts = multiTraj.getTrackState(lastSmoothedIndex);
  multiTraj.visitBackwards(lastSmoothedIndex, [&](const auto& ts) {
    const size_t iRow = fullGlobalTrackParamsCov.rows() -
                        eBoundParametersSize * (nProcessed + 1);
    // Fill the covariance of this state
    fullGlobalTrackParamsCov.block<eBoundParametersSize, eBoundParametersSize>(
        iRow, iRow) = ts.smoothedCovariance();
    // Fill the correlation between this state (indexed by i-1) and
    // beforehand smoothed states (indexed by j): C^n_{i-1, j}= G_{i-1} *
    // C^n_{i, j} for i <= j
    if (nProcessed > 0) {
      // Calculate the gain matrix
      gain_matrix_t G = ts.filteredCovariance() *
                        prev_ts.jacobian().transpose() *
                        prev_ts.predictedCovariance().inverse();
      // Loop over the beforehand smoothed states
      for (size_t iProcessed = 1; iProcessed <= nProcessed; iProcessed++) {
        const size_t iCol = iRow + eBoundParametersSize * iProcessed;
        CovMatrix_t prev_correlation =
            fullGlobalTrackParamsCov
                .block<eBoundParametersSize, eBoundParametersSize>(
                    iRow + eBoundParametersSize, iCol);
        CovMatrix_t correlation = G * prev_correlation;
        fullGlobalTrackParamsCov
            .block<eBoundParametersSize, eBoundParametersSize>(iRow, iCol) =
            correlation;
        fullGlobalTrackParamsCov
            .block<eBoundParametersSize, eBoundParametersSize>(iCol, iRow) =
            correlation.transpose();
      }
    }
    nProcessed++;
    prev_ts = ts;
  });

  // If necessary, extract only those elements for measurement states
  if (covarianceForMeasurementStates) {
    ActsMatrixX<BoundParametersScalar> shrinkedGlobalTrackParamsCov;
    size_t nMeasurementStates = smoothingOrders.size();
    // Set the size of the matrix
    shrinkedGlobalTrackParamsCov.resize(
        nMeasurementStates * eBoundParametersSize,
        nMeasurementStates * eBoundParametersSize);
    for (size_t i = 0; i < nMeasurementStates; i++) {
      for (size_t j = 0; j < nMeasurementStates; j++) {
        // Get the covariance/correlation
        size_t iRowFull =
            (nSmoothedStates - smoothingOrders.at(i)) * eBoundParametersSize;
        size_t iColFull =
            (nSmoothedStates - smoothingOrders.at(j)) * eBoundParametersSize;
        CovMatrix_t correlation =
            fullGlobalTrackParamsCov
                .block<eBoundParametersSize, eBoundParametersSize>(iRowFull,
                                                                   iColFull);
        // Fill the covariance/correlation
        size_t iRowShrinked =
            (nMeasurementStates - i - 1) * eBoundParametersSize;
        size_t iColShrinked =
            (nMeasurementStates - j - 1) * eBoundParametersSize;
        shrinkedGlobalTrackParamsCov
            .block<eBoundParametersSize, eBoundParametersSize>(
                iRowShrinked, iColShrinked) = correlation;
      }
    }
    return shrinkedGlobalTrackParamsCov;
  }

  return fullGlobalTrackParamsCov;
}

}  // namespace detail
}  // namespace Acts
