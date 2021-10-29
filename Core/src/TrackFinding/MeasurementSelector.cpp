// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/MeasurementSelector.hpp"

namespace Acts {

#define NEW

#ifdef NEW
Result<std::pair<std::vector<MultiTrajectory::TrackStateProxy>::iterator,
                 std::vector<MultiTrajectory::TrackStateProxy>::iterator>>
MeasurementSelector::select(
    std::vector<MultiTrajectory::TrackStateProxy>& candidates, bool& isOutlier,
    LoggerWrapper logger) const {
  using Result = Result<
      std::pair<std::vector<MultiTrajectory::TrackStateProxy>::iterator,
                std::vector<MultiTrajectory::TrackStateProxy>::iterator>>;

  ACTS_VERBOSE("Invoked MeasurementSelector");

  // Return error if no measurement
  if (candidates.empty()) {
    return CombinatorialKalmanFilterError::MeasurementSelectionFailed;
  }

  // Get geoID of this surface
  auto surface = &candidates.front().referenceSurface();
  auto geoID = surface->geometryId();

  // Find the appropriate cuts
  auto cuts = m_config.find(geoID);
  if (cuts == m_config.end()) {
    // for now we consider missing cuts an unrecoverable error
    // TODO consider other options e.g. do not add measurements at all (not
    // even as outliers)
    return CombinatorialKalmanFilterError::MeasurementSelectionFailed;
  }

  const double chi2CutOff = cuts->chi2CutOff;
  const size_t numMeasurementsCutOff = cuts->numMeasurementsCutOff;
  ACTS_VERBOSE("Allowed maximum chi2: " << chi2CutOff);
  ACTS_VERBOSE(
      "Allowed maximum number of measurements: " << numMeasurementsCutOff);
  ACTS_VERBOSE("Number of measurement candidates: " << candidates.size());

  // std::vector<std::pair<size_t, double>> measChi2;
  // measChi2.resize(candidates.size());
  double minChi2 = std::numeric_limits<double>::max();
  size_t minIndex = 0;
  size_t index = 0;
  size_t nInitialCandidates = 0;
  // Loop over all measurements to select the compatible measurements
  for (auto& trackState : candidates) {
    // Take the parameter covariance
    const auto predicted = trackState.predicted();
    const auto predictedCovariance = trackState.predictedCovariance();
    visit_measurement(
        trackState.calibrated(), trackState.calibratedCovariance(),
        trackState.calibratedSize(),
        [&](const auto calibrated, const auto calibratedCovariance) {
          constexpr size_t kMeasurementSize =
              decltype(calibrated)::RowsAtCompileTime;

          using ParametersVector = ActsVector<kMeasurementSize>;

          // Take the projector (measurement mapping function)
          const auto H =
              trackState.projector()
                  .template topLeftCorner<kMeasurementSize, eBoundSize>()
                  .eval();

          // Get the residuals
          ParametersVector res;
          res = calibrated - H * predicted;

          // Get the chi2
          double& chi2 = trackState.chi2();
          chi2 = (res.transpose() *
                  ((calibratedCovariance +
                    H * predictedCovariance * H.transpose()))
                      .inverse() *
                  res)
                     .eval()(0, 0);
          ACTS_VERBOSE("Chi2: " << chi2);

          // use track state chi2 storage for this
          if (chi2 < chi2CutOff) {
            // // // measChi2.at(nInitialCandidates) = {index, chi2};
            nInitialCandidates++;
          }
          // Search for the measurement with the min chi2
          if (chi2 < minChi2) {
            minChi2 = chi2;
            minIndex = index;
          }
        });

    index++;
  }

  // If there is no selected measurement, return the measurement with the best
  // chi2 and tag it as an outlier
  if (nInitialCandidates == 0) {
    ACTS_VERBOSE("No measurement candidate. Return an outlier measurement.");
    isOutlier = true;
    auto outlierIt = std::next(candidates.begin(), minIndex);
    // reset this state's chi2
    // outlierIt->chi2() = 0;
    // return single item range, no sorting necessary
    return Result::success(std::pair{outlierIt, std::next(outlierIt, 1)});
  }

  // Get the number of measurement candidates with provided constraint
  // considered
  // size_t nFinalCandidates =
  // std::min(nInitialCandidates, numMeasurementsCutOff);

  std::sort(
      candidates.begin(), candidates.end(),
      [](const auto& tsa, const auto& tsb) { return tsa.chi2() < tsb.chi2(); });

  auto endIterator = candidates.begin();
  auto maxIterator = candidates.end();
  if (candidates.size() > numMeasurementsCutOff) {
    maxIterator = std::next(candidates.begin(), numMeasurementsCutOff);
  }

  for (; endIterator != maxIterator; ++endIterator) {
    if (endIterator->chi2() >= chi2CutOff) {
      break;  // endIterator now points at the first track state with chi2
              // larger than our cutoff => defines the end of our returned
              // range
    }
  }
  ACTS_VERBOSE("Number of selected measurements: " << std::distance(
                   candidates.begin(), endIterator));

  // auto maxMeasIterator = std::next(candidates.begin() + nFinalCandidates);

  // Sort the initial measurement candidates based on chi2 in
  // ascending order
  // std::sort(measChi2.begin(), measChi2.begin() + nInitialCandidates,
  // [](const std::pair<size_t, double>& lchi2,
  // const std::pair<size_t, double>& rchi2) {
  // return lchi2.second < rchi2.second;
  // });

  // @TODO: Get rid of this
  // std::vector<MultiTrajectory::TrackStateProxy> trackStatesCopy;
  // trackStatesCopy = candidates;
  // candidates.clear();

  // Get only allowed number of measurement candidates, i.e. nFinalCandidates,
  // from the front and reset the values in the container
  // size_t nRecorded = 0;
  // for (const auto& [id, chi2] : measChi2) {
  // if (nRecorded >= nFinalCandidates) {
  // break;
  // }
  // candidates.push_back(trackStatesCopy.at(id));
  // nRecorded++;
  // }

  isOutlier = false;
  return std::pair{candidates.begin(), endIterator};
}

#else

Result<std::pair<std::vector<MultiTrajectory::TrackStateProxy>::iterator,
                 std::vector<MultiTrajectory::TrackStateProxy>::iterator>>
MeasurementSelector::select(
    std::vector<MultiTrajectory::TrackStateProxy>& candidates, bool& isOutlier,
    LoggerWrapper logger) const {
  using Result = Result<
      std::pair<std::vector<MultiTrajectory::TrackStateProxy>::iterator,
                std::vector<MultiTrajectory::TrackStateProxy>::iterator>>;

  ACTS_VERBOSE("Invoked MeasurementSelector");

  // Return error if no measurement
  if (candidates.empty()) {
    return CombinatorialKalmanFilterError::MeasurementSelectionFailed;
  }

  // Get geoID of this surface
  auto surface = &candidates.front().referenceSurface();
  auto geoID = surface->geometryId();

  // Find the appropriate cuts
  auto cuts = m_config.find(geoID);
  if (cuts == m_config.end()) {
    // for now we consider missing cuts an unrecoverable error
    // TODO consider other options e.g. do not add measurements at all (not
    // even as outliers)
    return CombinatorialKalmanFilterError::MeasurementSelectionFailed;
  }

  const double chi2CutOff = cuts->chi2CutOff;
  const size_t numMeasurementsCutOff = cuts->numMeasurementsCutOff;
  ACTS_VERBOSE("Allowed maximum chi2: " << chi2CutOff);
  ACTS_VERBOSE(
      "Allowed maximum number of measurements: " << numMeasurementsCutOff);

  std::vector<std::pair<size_t, double>> measChi2;
  measChi2.resize(candidates.size());
  double minChi2 = std::numeric_limits<double>::max();
  size_t minIndex = 0;
  size_t index = 0;
  size_t nInitialCandidates = 0;
  // Loop over all measurements to select the compatible measurements
  for (const auto& trackState : candidates) {
    // Take the parameter covariance
    const auto predicted = trackState.predicted();
    const auto predictedCovariance = trackState.predictedCovariance();
    visit_measurement(
        trackState.calibrated(), trackState.calibratedCovariance(),
        trackState.calibratedSize(),
        [&](const auto calibrated, const auto calibratedCovariance) {
          constexpr size_t kMeasurementSize =
              decltype(calibrated)::RowsAtCompileTime;

          using ParametersVector = ActsVector<kMeasurementSize>;

          // Take the projector (measurement mapping function)
          const auto H =
              trackState.projector()
                  .template topLeftCorner<kMeasurementSize, eBoundSize>()
                  .eval();

          // Get the residuals
          ParametersVector res;
          res = calibrated - H * predicted;
          // ACTS_VERBOSE("Residual: " << res.transpose());

          // Get the chi2
          double chi2 = (res.transpose() *
                         ((calibratedCovariance +
                           H * predictedCovariance * H.transpose()))
                             .inverse() *
                         res)
                            .eval()(0, 0);
          ACTS_VERBOSE("Chi2: " << chi2);
          // Push the measurement index and chi2 if satisfying the criteria
          if (chi2 < chi2CutOff) {
            measChi2.at(nInitialCandidates) = {index, chi2};
            nInitialCandidates++;
          }
          // Search for the measurement with the min chi2
          if (chi2 < minChi2) {
            minChi2 = chi2;
            minIndex = index;
          }
        });

    index++;
  }

  // Get the number of measurement candidates with provided constraint
  // considered
  size_t nFinalCandidates = std::min(nInitialCandidates, numMeasurementsCutOff);

  // If there is no selected measurement, return the measurement with the best
  // chi2 and tag it as an outlier
  if (nFinalCandidates == 0) {
    ACTS_VERBOSE("No measurement candidate. Return an outlier measurement.");
    isOutlier = true;
    // return single item range, no sorting necessary
    return Result::success(std::pair{candidates.begin() + minIndex,
                                     candidates.begin() + minIndex + 1});
  }

  ACTS_VERBOSE("Number of measurement candidates: " << nFinalCandidates);
  // Sort the initial measurement candidates based on chi2 in
  // ascending order
  std::sort(measChi2.begin(),
            // measChi2.end(),
            measChi2.begin() + nInitialCandidates,
            [](const std::pair<size_t, double>& lchi2,
               const std::pair<size_t, double>& rchi2) {
              return lchi2.second < rchi2.second;
            });

  // @TODO: Get rid of this
  std::vector<MultiTrajectory::TrackStateProxy> trackStatesCopy;
  trackStatesCopy = candidates;
  candidates.clear();

  // Get only allowed number of measurement candidates, i.e. nFinalCandidates,
  // from the front and reset the values in the container
  size_t nRecorded = 0;
  for (const auto& [id, chi2] : measChi2) {
    if (nRecorded >= nFinalCandidates) {
      break;
    }
    candidates.push_back(trackStatesCopy.at(id));
    trackStatesCopy.at(id).chi2() = 42.0;
    nRecorded++;
  }

  isOutlier = false;
  return std::pair{candidates.begin(), candidates.end()};
}
#endif

}  // namespace Acts
