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
    : m_config{{GeometryIdentifier(), MeasurementSelectorCuts{}}} {}

MeasurementSelector::MeasurementSelector(const MeasurementSelectorCuts& cuts)
    : m_config{{GeometryIdentifier(), cuts}} {}

MeasurementSelector::MeasurementSelector(Config config)
    : m_config(std::move(config)) {}

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

  auto select = [](const auto& vec, std::size_t idx) {
    return vec.at(std::clamp(idx, std::size_t{0}, vec.size() - 1));
  };

  const double chi2CutOffMeasurement = select(config.chi2CutOff, bin);
  const double chi2CutOffOutlier = config.chi2CutOffOutlier.empty()
                                       ? std::numeric_limits<double>::infinity()
                                       : select(config.chi2CutOffOutlier, bin);
  const std::size_t numMeasurementsCutOff =
      select(config.numMeasurementsCutOff, bin);
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
