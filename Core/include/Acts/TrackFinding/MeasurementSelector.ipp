// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFinding/MeasurementSelector.hpp"

#include <cmath>
#include <optional>

namespace Acts {

template <typename traj_t>
Result<
    std::pair<typename std::vector<typename traj_t::TrackStateProxy>::iterator,
              typename std::vector<typename traj_t::TrackStateProxy>::iterator>>
MeasurementSelector::select(
    const BoundParametrization& parametrization,
    std::vector<typename traj_t::TrackStateProxy>& candidates, bool& isOutlier,
    const Logger& logger) const {
  using Result = Result<std::pair<
      typename std::vector<typename traj_t::TrackStateProxy>::iterator,
      typename std::vector<typename traj_t::TrackStateProxy>::iterator>>;

  ACTS_VERBOSE("Invoked MeasurementSelector");

  // Return if no measurement
  if (candidates.empty()) {
    return Result::success(std::pair(candidates.begin(), candidates.end()));
  }

  // Candidates normally share one prediction. Find a usable one before
  // invoking parametrization operations, which require finite coordinates.
  // Keep the found prediction to avoid repeating the component lookup.
  auto prediction = candidates.begin();
  std::optional<decltype(prediction->predicted())> predicted;
  for (; prediction != candidates.end(); ++prediction) {
    if (auto candidatePredicted = prediction->predicted();
        detail::allFinite(candidatePredicted)) {
      predicted.emplace(candidatePredicted);
      break;
    }
  }
  if (!predicted.has_value()) {
    for (auto& candidate : candidates) {
      candidate.chi2() = std::numeric_limits<double>::quiet_NaN();
    }
    isOutlier = false;
    return Result::success(std::pair(candidates.begin(), candidates.begin()));
  }
  // All candidates are on the same surface, which also selects the cuts.
  // The caller resolves the surface's parametrization.
  const GeometryIdentifier geoID = prediction->referenceSurface().geometryId();
  const double theta = parametrization.theta(*predicted);
  // Find the appropriate cuts
  const auto cutsResult = getCuts(geoID, theta);
  if (!cutsResult.ok()) {
    return cutsResult.error();
  }
  const Cuts& cuts = *cutsResult;

  if (cuts.numMeasurements == 0ul) {
    return CombinatorialKalmanFilterError::MeasurementSelectionFailed;
  }

  double minChi2 = std::numeric_limits<double>::max();
  std::size_t minIndex = std::numeric_limits<std::size_t>::max();

  isOutlier = false;

  // Loop over all measurements to select the compatible measurements
  // Sort track states which do not satisfy the chi2 cut to the end.
  // When done trackStateIterEnd will point to the first element that
  // does not satisfy the chi2 cut.
  std::size_t passedCandidates = 0ul;
  // Candidates normally share the predicted storage of the first one, so the
  // finiteness of a prediction only has to be checked once per distinct
  // buffer. The predicted component is a writable view of the backend storage,
  // so equal data pointers refer to the same, here unmodified, parameters.
  // Candidates with separate storage, e.g. from a custom track state creator,
  // are checked individually.
  const double* checkedPredicted = predicted->data();
  bool checkedPredictedFinite = true;
  for (std::size_t i(0ul); i < candidates.size(); ++i) {
    auto& trackState = candidates[i];

    const auto candidatePredicted = trackState.predicted();
    if (candidatePredicted.data() != checkedPredicted) {
      checkedPredicted = candidatePredicted.data();
      checkedPredictedFinite = detail::allFinite(candidatePredicted);
    }

    // We access the dynamic size of the matrix here but use them later
    // through a template function which accesses the data pointer
    // with compile time size. That way the Eigen math operations are
    // still done with compile time size and no dynamic memory allocation
    // is needed.
    // A non-finite prediction must not reach the parametrization operations,
    // which require finite coordinates.
    double chi2 =
        checkedPredictedFinite
            ? calculateChi2(
                  parametrization, trackState.effectiveCalibrated().data(),
                  trackState.effectiveCalibratedCovariance().data(),
                  candidatePredicted, trackState.predictedCovariance(),
                  trackState.projectorSubspaceIndices(),
                  trackState.calibratedSize())
            : std::numeric_limits<double>::quiet_NaN();
    trackState.chi2() = chi2;

    // A diverging fit can inflate the predicted covariance until inverting it
    // overflows, making the chi2 non-finite. Such a candidate has to be
    // rejected explicitly: every comparison against NaN is false, so it would
    // otherwise pass the chi2 cut below while leaving minIndex unset, and it
    // would make the comparator of the sort below inconsistent.
    if (!std::isfinite(chi2)) {
      ACTS_VERBOSE("Non-finite chi2 " << chi2 << " on surface " << geoID
                                      << ", rejecting measurement candidate");
      continue;
    }

    if (chi2 < minChi2) {
      minChi2 = chi2;
      minIndex = i;
    }

    // only consider track states which pass the chi2 cut
    if (chi2 >= cuts.chi2Measurement[trackState.calibratedSize()]) {
      continue;
    }

    if (passedCandidates != i) {
      std::swap(candidates[passedCandidates], candidates[i]);
      if (minIndex == i) {
        minIndex = passedCandidates;
      }
    }
    ++passedCandidates;
  }

  // Handle if there are no measurements below the chi2 cut off
  if (passedCandidates == 0ul) {
    // minIndex stays unset if every candidate was rejected above.
    if (minIndex < candidates.size() &&
        minChi2 < cuts.chi2Outlier[candidates[minIndex].calibratedSize()]) {
      ACTS_VERBOSE(
          "No measurement candidate. Return an outlier measurement chi2="
          << minChi2);
      isOutlier = true;
      return Result::success(std::pair(candidates.begin() + minIndex,
                                       candidates.begin() + minIndex + 1));
    } else {
      ACTS_VERBOSE("No measurement candidate. Return empty chi2=" << minChi2);
      return Result::success(std::pair(candidates.begin(), candidates.begin()));
    }
  }

  if (passedCandidates <= 1ul) {
    // Return the single item range, so no sorting is necessary. The loop moves
    // the passed candidate to the front. It is not always the candidate with
    // minChi2, because the chi2 cut can depend on the measurement dimension.
    ACTS_VERBOSE("Returning only 1 element chi2=" << candidates.front().chi2());
    return Result::success(
        std::pair(candidates.begin(), candidates.begin() + 1));
  }

  std::sort(
      candidates.begin(), candidates.begin() + passedCandidates,
      [](const auto& tsa, const auto& tsb) { return tsa.chi2() < tsb.chi2(); });

  ACTS_VERBOSE("Number of selected measurements: "
               << passedCandidates << ", max: " << cuts.numMeasurements);

  return Result::success(std::pair(
      candidates.begin(),
      candidates.begin() + std::min(cuts.numMeasurements, passedCandidates)));
}

}  // namespace Acts
