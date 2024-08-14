// This file is part of the Acts project.
//
// Copyright (C) 2021-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/MeasurementSelector.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"

#include <algorithm>
#include <limits>

namespace Acts {

MeasurementSelector::MeasurementSelector()
    : MeasurementSelector({{GeometryIdentifier(), MeasurementSelectorCuts{}}}) {
}

MeasurementSelector::MeasurementSelector(const MeasurementSelectorCuts& cuts)
    : MeasurementSelector({{GeometryIdentifier(), cuts}}) {}

MeasurementSelector::MeasurementSelector(Config config)
    : m_config(std::move(config)) {
  for (const auto& cuts : m_config) {
    validateCuts(cuts);
  }
}

void MeasurementSelector::validateCuts(const MeasurementSelectorCuts& cuts) {
  if (!cuts.chi2CutOffOutlier.empty() &&
      cuts.chi2CutOff.size() != cuts.chi2CutOffOutlier.size()) {
    throw std::invalid_argument(
        "chi2CutOff and chi2CutOffOutlier must have the same size");
  }

  if (cuts.chi2CutOff.size() != cuts.numMeasurementsCutOff.size()) {
    throw std::invalid_argument(
        "chi2CutOff and numMeasurementsCutOff must have the same size");
  }

  if (cuts.chi2CutOff.size() != cuts.etaBins.size() + 1) {
    throw std::invalid_argument(
        "chi2CutOff must have one more element than etaBins");
  }
}

double MeasurementSelector::calculateChi2(
    const double* fullCalibrated, const double* fullCalibratedCovariance,
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Parameters predicted,
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Covariance predictedCovariance,
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Projector projector,
    unsigned int calibratedSize) const {
  return visit_measurement(
      calibratedSize,
      [&fullCalibrated, &fullCalibratedCovariance, &predicted,
       &predictedCovariance, &projector](auto N) -> double {
        constexpr std::size_t kMeasurementSize = decltype(N)::value;

        typename TrackStateTraits<kMeasurementSize, true>::Calibrated
            calibrated{fullCalibrated};

        typename TrackStateTraits<kMeasurementSize, true>::CalibratedCovariance
            calibratedCovariance{fullCalibratedCovariance};

        using ParametersVector = ActsVector<kMeasurementSize>;

        // Take the projector (measurement mapping function)
        const auto H =
            projector.template topLeftCorner<kMeasurementSize, eBoundSize>()
                .eval();

        // Get the residuals
        ParametersVector res = calibrated - H * predicted;

        // Get the chi2
        return (res.transpose() *
                (calibratedCovariance + H * predictedCovariance * H.transpose())
                    .inverse() *
                res)
            .eval()(0, 0);
      });
}

MeasurementSelector::Cuts MeasurementSelector::getCutsByEta(
    const MeasurementSelectorCuts& config, double eta) {
  const double etaAbs = std::abs(eta);
  std::size_t bin = 0;
  for (; bin < config.etaBins.size(); bin++) {
    if (config.etaBins[bin] >= etaAbs) {
      break;
    }
  }

  const double chi2CutOffMeasurement = config.chi2CutOff.at(bin);
  const double chi2CutOffOutlier = config.chi2CutOffOutlier.empty()
                                       ? std::numeric_limits<double>::infinity()
                                       : config.chi2CutOffOutlier.at(bin);
  const std::size_t numMeasurementsCutOff =
      config.numMeasurementsCutOff.at(bin);
  return {chi2CutOffMeasurement, chi2CutOffOutlier, numMeasurementsCutOff};
}

Result<MeasurementSelector::Cuts> MeasurementSelector::getCuts(
    const GeometryIdentifier& geoID, double eta) const {
  // Find the appropriate cuts
  auto cuts = m_config.find(geoID);
  if (cuts == m_config.end()) {
    // for now we consider missing cuts an unrecoverable error
    // TODO consider other options e.g. do not add measurements at all (not
    // even as outliers)
    return CombinatorialKalmanFilterError::MeasurementSelectionFailed;
  }
  assert(!cuts->chi2CutOff.empty());
  return getCutsByEta(*cuts, eta);
}

}  // namespace Acts
