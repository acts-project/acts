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

#include <unordered_map>

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
///
/// @return The global track parameters covariance matrix and the starting
/// row/column for smoothed states
template <typename source_link_t, typename parameters_t = BoundParameters>
std::pair<ActsMatrixX<BoundParametersScalar>,
          std::unordered_map<size_t, size_t>>
globalTrackParametersCovariance(
    const Acts::MultiTrajectory<source_link_t>& multiTraj,
    const size_t& entryIndex) {
  using CovMatrix = typename parameters_t::CovMatrix_t;
  using GainMatrix = CovMatrix;

  // The last smoothed state index
  size_t lastSmoothedIndex = SIZE_MAX;
  // The total number of smoothed states
  size_t nSmoothedStates = 0;
  // Visit all the states
  multiTraj.visitBackwards(entryIndex, [&](const auto& ts) {
    if (ts.hasSmoothed()) {
      if (lastSmoothedIndex == SIZE_MAX) {
        lastSmoothedIndex = ts.index();
      }
      nSmoothedStates++;
    }
  });

  // Set the size of global track parameters covariance for all smoothed states
  ActsMatrixX<BoundParametersScalar> fullGlobalTrackParamsCov =
      ActsMatrixX<BoundParametersScalar>::Zero(
          nSmoothedStates * eBoundParametersSize,
          nSmoothedStates * eBoundParametersSize);
  // The index of state within the trajectory and the starting row/column for
  // this state in the global covariance matrix
  std::unordered_map<size_t, size_t> stateRowIndices;
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
      GainMatrix G = ts.filteredCovariance() * prev_ts.jacobian().transpose() *
                     prev_ts.predictedCovariance().inverse();
      // Loop over the beforehand smoothed states
      for (size_t iProcessed = 1; iProcessed <= nProcessed; iProcessed++) {
        const size_t iCol = iRow + eBoundParametersSize * iProcessed;
        CovMatrix prev_correlation =
            fullGlobalTrackParamsCov
                .block<eBoundParametersSize, eBoundParametersSize>(
                    iRow + eBoundParametersSize, iCol);
        CovMatrix correlation = G * prev_correlation;
        fullGlobalTrackParamsCov
            .block<eBoundParametersSize, eBoundParametersSize>(iRow, iCol) =
            correlation;
        fullGlobalTrackParamsCov
            .block<eBoundParametersSize, eBoundParametersSize>(iCol, iRow) =
            correlation.transpose();
      }
    }
    stateRowIndices.emplace(ts.index(), iRow);
    nProcessed++;
    prev_ts = ts;
  });

  return std::make_pair(fullGlobalTrackParamsCov, stateRowIndices);
}

}  // namespace detail
}  // namespace Acts
