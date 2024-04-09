// This file is part of the Acts project.
//
// Copyright (C) 2019-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterError.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <cassert>
#include <cstddef>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>
#include <chrono>
#include <iostream>

namespace Acts {

/// Selection cuts for associating measurements with predicted track
/// parameters on a surface.
///
/// The default configuration only takes the best matching measurement without a
/// cut on the local chi2.
struct MeasurementSelectorCuts {
  /// bins in |eta| to specify variable selections
  std::vector<double> etaBins{};
  /// Maximum local chi2 contribution.
  std::vector<double> chi2CutOff{std::numeric_limits<double>::max()};
  /// Maximum number of associated measurements on a single surface.
  std::vector<std::size_t> numMeasurementsCutOff{1};
};

/// @brief Measurement selection struct selecting those measurements compatible
/// with the given track parameter against provided criteria on one surface
///
/// The selection criteria could be allowed maximum chi2
/// and allowed maximum number of measurements on one surface
///
/// If there is no compatible measurement, the measurement with the minimum
/// chi2 will be selected and the status will be tagged as an outlier
///
class MeasurementSelector {
 public:
  /// Geometry-dependent cut configuration.
  ///
  /// Different components on the geometry can require different cut settings.
  /// The configuration must either contain explicit settings for all geometry
  /// components that are used or contain a global default.
  using Config = Acts::GeometryHierarchyMap<MeasurementSelectorCuts>;

  /// @brief Default constructor
  MeasurementSelector() = default;
  /// @brief Constructor with config and (non-owning) logger
  ///
  /// @param config a config instance
  MeasurementSelector(Config config);

  /// @brief Function that select the measurements compatible with
  /// the given track parameter on a surface
  ///
  /// @param candidates The track state candidates which already contain predicted parameters
  /// @param isOutlier The indicator for outlier or not
  /// @param logger The logger wrapper
  ///
  /// @return Pair of iterators into @a candidates marking the range of selected candidates
  ///

  template <typename traj_t>
  Result<std::pair<
      typename std::vector<typename traj_t::TrackStateProxy>::iterator,
      typename std::vector<typename traj_t::TrackStateProxy>::iterator>>
  select(std::vector<typename traj_t::TrackStateProxy>& candidates,
         bool& isOutlier, const Logger& logger) const;

private:
  template <typename traj_t, typename cut_value_t>
  static cut_value_t VariableCut(
      const typename traj_t::TrackStateProxy& trackState,
      const Acts::MeasurementSelector::Config::Iterator selector,
      const std::vector<cut_value_t>& cuts, const Logger& logger);

  double calculateChi2(
      const double* fullCalibrated, const double* fullCalibratedCovariance,
      const TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                       false>::Parameters predicted,
      const TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                       false>::Covariance predictedCovariance,
      const TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                       false>::Projector projector,
      const unsigned int calibratedSize) const;

  Config m_config;
};

}  // namespace Acts

#include "MeasurementSelector.ipp"
