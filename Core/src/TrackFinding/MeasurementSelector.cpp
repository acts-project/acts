// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/MeasurementSelector.hpp"

namespace Acts {

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

  const auto& chi2CutOff = cuts->chi2CutOff;
  const auto& numMeasurementsCutOff = cuts->numMeasurementsCutOff;

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
          const auto chi2Cut =
              VariableCut(predicted, cuts, chi2CutOff, logger);
          ACTS_VERBOSE("Chi2: " << chi2 << ", max: " << chi2Cut);

          // use track state chi2 storage for this
          if (chi2 < chi2Cut) {
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
    // return single item range, no sorting necessary
    return Result::success(std::pair{outlierIt, std::next(outlierIt, 1)});
  }

  std::sort(
      candidates.begin(), candidates.end(),
      [](const auto& tsa, const auto& tsb) { return tsa.chi2() < tsb.chi2(); });

  const size_t numMeasurementsCut =
      VariableCut(candidates.begin()->predicted(), cuts, numMeasurementsCutOff, logger);

  auto endIterator = candidates.begin();
  auto maxIterator = candidates.end();
  if (candidates.size() > numMeasurementsCut) {
    maxIterator = std::next(candidates.begin(), numMeasurementsCut);
  }

  for (; endIterator != maxIterator; ++endIterator) {
    const auto chi2Cut = VariableCut(endIterator->predicted(), cuts, chi2CutOff, logger);
    if (endIterator->chi2() >= chi2Cut) {
      break;  // endIterator now points at the first track state with chi2
              // larger than our cutoff => defines the end of our returned
              // range
    }
  }
  ACTS_VERBOSE("Number of selected measurements: "
               << std::distance(candidates.begin(), endIterator)
               << ", max: " << numMeasurementsCut);

  isOutlier = false;
  return std::pair{candidates.begin(), endIterator};
}

template <typename cut_value_t>
static cut_value_t VariableCut(
    const Acts::MultiTrajectory::TrackStateProxy::Parameters& predictedParams,
    const Acts::MeasurementSelector::Config::Iterator selector,
    const std::vector<cut_value_t>& cuts, LoggerWrapper logger) {
  const auto& etaBins = selector->etaBins;
  if (etaBins.empty()) {
    return cuts[0];  // shortcut if no etaBins
  }
  const auto eta = Acts::VectorHelpers::eta(predictedParams.unitDirection());
  const auto abseta = std::abs(eta);
  size_t bin = 0;
  for (auto etaBin : etaBins) {
    if (etaBin >= abseta) {
      break;
    }
    bin++;
  }
  if (bin >= cuts.size()) {
    bin = cuts.size() - 1;
  }
  ACTS_VERBOSE("Variable cut for eta=" << eta << ": " << cuts[bin]);
  return cuts[bin];
}

}  // namespace Acts
