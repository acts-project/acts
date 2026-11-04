// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/MeasurementSelector.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>

#include <boost/math/distributions/chi_squared.hpp>

namespace Acts {

MeasurementSelector::MeasurementSelector()
    : MeasurementSelector({{GeometryIdentifier(), MeasurementSelectorCuts{}}}) {
}

MeasurementSelector::MeasurementSelector(const MeasurementSelectorCuts& cuts)
    : MeasurementSelector({{GeometryIdentifier(), cuts}}) {}

MeasurementSelector::MeasurementSelector(const Config& config) {
  if (config.empty()) {
    throw std::invalid_argument(
        "MeasurementSelector: Configuration must not be empty");
  }

  std::vector<InternalConfig::InputElement> tmp;
  tmp.reserve(config.size());
  for (std::size_t i = 0; i < config.size(); ++i) {
    GeometryIdentifier geoID = config.idAt(i);
    MeasurementSelectorCuts cuts = config.valueAt(i);

    InternalCutBins internalCuts = convertCutBins(cuts);
    tmp.emplace_back(geoID, std::move(internalCuts));
  }
  m_config = InternalConfig(std::move(tmp));
}

namespace {

/// Convert the chi2 cut @p cut1D for one-dimensional measurements to the cut
/// with the same p-value for each measurement dimension. Without
/// @p perDimension every dimension uses @p cut1D.
std::array<double, kMeasurementSizeMax + 1> convertChi2Cut(double cut1D,
                                                           bool perDimension) {
  std::array<double, kMeasurementSizeMax + 1> result{};
  result.fill(cut1D);
  if (!perDimension || !std::isfinite(cut1D) ||
      cut1D == std::numeric_limits<double>::max() || cut1D <= 0) {
    return result;
  }

  const double pValue = boost::math::cdf(boost::math::complement(
      boost::math::chi_squared_distribution<>(1), cut1D));
  // A cut far in the tail has no representable p-value, so it stays as it is
  if (!(pValue > std::numeric_limits<double>::min())) {
    return result;
  }
  for (std::size_t dim = 2; dim < result.size(); ++dim) {
    result[dim] = boost::math::quantile(boost::math::complement(
        boost::math::chi_squared_distribution<>(static_cast<double>(dim)),
        pValue));
  }
  return result;
}

}  // namespace

MeasurementSelector::InternalCutBins MeasurementSelector::convertCutBins(
    const MeasurementSelectorCuts& config) {
  InternalCutBins cutBins;

  auto getEtaOrInf = [](const auto& vec, std::size_t bin) {
    if (bin >= vec.size()) {
      return std::numeric_limits<double>::infinity();
    }
    assert(vec[bin] >= 0 && "Eta bins must be positive");
    return vec[bin];
  };

  auto getBinOrBackOrMax = [](const auto& vec, std::size_t bin) {
    using Value = std::remove_reference_t<decltype(vec[0])>;
    static constexpr Value max = std::numeric_limits<Value>::max();
    return vec.empty() ? max : (bin < vec.size() ? vec[bin] : vec.back());
  };

  for (std::size_t bin = 0; bin < config.etaBins.size() + 1; ++bin) {
    InternalCutBin cuts;
    cuts.maxTheta = getEtaOrInf(config.etaBins, bin);
    cuts.maxNumMeasurements =
        getBinOrBackOrMax(config.numMeasurementsCutOff, bin);
    cuts.maxChi2Measurement = convertChi2Cut(
        getBinOrBackOrMax(config.chi2CutOff, bin), config.chi2CutPerDimension);
    cuts.maxChi2Outlier =
        convertChi2Cut(getBinOrBackOrMax(config.chi2CutOffOutlier, bin),
                       config.chi2CutPerDimension);
    cutBins.push_back(cuts);
  }

  return cutBins;
}

double MeasurementSelector::calculateChi2(
    const double* fullCalibrated, const double* fullCalibratedCovariance,
    TrackStateTraits<kMeasurementSizeMax, false>::Parameters predicted,
    TrackStateTraits<kMeasurementSizeMax, false>::Covariance
        predictedCovariance,
    BoundSubspaceIndices projector, unsigned int calibratedSize) const {
  return visit_measurement(
      calibratedSize,
      [&fullCalibrated, &fullCalibratedCovariance, &predicted,
       &predictedCovariance, &projector](auto N) -> double {
        constexpr std::size_t kMeasurementSize = decltype(N)::value;

        typename TrackStateTraits<kMeasurementSize, true>::Calibrated
            calibrated{fullCalibrated};

        typename TrackStateTraits<kMeasurementSize, true>::CalibratedCovariance
            calibratedCovariance{fullCalibratedCovariance};

        using ParametersVector = Vector<kMeasurementSize>;

        std::span<const std::uint8_t, kMeasurementSize> validSubspaceIndices(
            projector.begin(), projector.begin() + kMeasurementSize);
        FixedBoundSubspaceHelper<kMeasurementSize> subspaceHelper(
            validSubspaceIndices);

        // Get the residuals
        ParametersVector res =
            calibrated - subspaceHelper.projectVector(predicted);

        // Get the chi2
        return (res.transpose() *
                (calibratedCovariance +
                 subspaceHelper.projectMatrix(predictedCovariance))
                    .inverse() *
                res)
            .eval()(0, 0);
      });
}

MeasurementSelector::Cuts MeasurementSelector::getCutsByTheta(
    const InternalCutBins& config, double theta) {
  // since theta is in [0, pi] and we have a symmetric cut in eta, we can just
  // look at the positive half of the Z axis
  const double constrainedTheta = std::min(theta, std::numbers::pi - theta);

  auto it = std::ranges::find_if(
      config, [constrainedTheta](const InternalCutBin& cuts) {
        return constrainedTheta < cuts.maxTheta;
      });
  assert(it != config.end());
  return {it->maxNumMeasurements, it->maxChi2Measurement, it->maxChi2Outlier};
}

Result<MeasurementSelector::Cuts> MeasurementSelector::getCuts(
    const GeometryIdentifier& geoID, double theta) const {
  // Find the appropriate cuts
  auto cuts = m_config.find(geoID);
  if (cuts == m_config.end()) {
    // for now we consider missing cuts an unrecoverable error
    // TODO consider other options e.g. do not add measurements at all (not
    // even as outliers)
    return CombinatorialKalmanFilterError::MeasurementSelectionFailed;
  }
  return getCutsByTheta(*cuts, theta);
}

}  // namespace Acts
