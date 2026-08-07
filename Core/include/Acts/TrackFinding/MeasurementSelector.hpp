// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterError.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cassert>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace Acts {

/// Local chi2 contribution of a calibrated measurement with respect to
/// predicted bound track parameters.
///
/// The calibrated data is passed as raw pointers into an overallocated buffer
/// so that the Eigen operations can be performed with compile time sizes and
/// no dynamic memory allocation is needed.
///
/// @param calibrated pointer to @p calibratedSize contiguous measurement values
/// @param calibratedCovariance pointer to the contiguous measurement covariance
/// @param predicted the predicted bound parameters
/// @param predictedCovariance the predicted bound parameters covariance
/// @param projector the subspace indices of the measurement
/// @param calibratedSize the dimension of the measurement
///
/// @return the local chi2 contribution
double calculatePredictedChi2(const double* calibrated,
                              const double* calibratedCovariance,
                              Eigen::Ref<const BoundVector> predicted,
                              Eigen::Ref<const BoundMatrix> predictedCovariance,
                              BoundSubspaceIndices projector,
                              unsigned int calibratedSize);

/// Selection cuts for associating measurements with predicted track
/// parameters on a surface.
///
/// The default configuration only takes the best matching measurement without a
/// cut on the local chi2.
struct MeasurementSelectorCuts {
  /// bins in |eta| to specify variable selections
  std::vector<double> etaBins{};
  /// Maximum local chi2 contribution to classify as measurement.
  std::vector<double> chi2CutOff{15};
  /// Maximum number of associated measurements on a single surface.
  std::vector<std::size_t> numMeasurementsCutOff{1};
  /// Maximum local chi2 contribution to classify as outlier.
  std::vector<double> chi2CutOffOutlier{};
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
  ///
  /// This will use the default configuration for the cuts.
  MeasurementSelector();

  /// @brief Constructor with cuts
  ///
  /// @param cuts The cuts to use
  explicit MeasurementSelector(const MeasurementSelectorCuts& cuts);

  /// @brief Constructor with config
  ///
  /// @param config a config instance
  explicit MeasurementSelector(const Config& config);

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

  /// The selection cuts resolved for one surface and one polar angle.
  struct Cuts {
    /// Maximum number of associated measurements on this surface
    std::size_t numMeasurements{};
    /// Maximum local chi2 contribution to classify as measurement
    double chi2Measurement{};
    /// Maximum local chi2 contribution to classify as outlier
    double chi2Outlier{};
  };

  /// Resolve the selection cuts for a surface and a polar angle.
  ///
  /// @param geoID the geometry identifier of the surface
  /// @param theta the polar angle of the predicted parameters
  ///
  /// @return the cuts or `MeasurementSelectionFailed` if the geometry
  ///         identifier is not covered by the configuration
  Result<Cuts> getCuts(const GeometryIdentifier& geoID, double theta) const;

 private:
  struct InternalCutBin {
    double maxTheta{};
    std::size_t maxNumMeasurements{};
    double maxChi2Measurement{};
    double maxChi2Outlier{};
  };
  using InternalCutBins = std::vector<InternalCutBin>;
  using InternalConfig = Acts::GeometryHierarchyMap<InternalCutBins>;

  static InternalCutBins convertCutBins(const MeasurementSelectorCuts& config);

  static Cuts getCutsByTheta(const InternalCutBins& config, double theta);

  InternalConfig m_config;
};

}  // namespace Acts

#include "MeasurementSelector.ipp"
