// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
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

namespace Acts {

/// Selection cuts for associating measurements with predicted track
/// parameters on a surface.
///
/// The default configuration only takes the best matching measurement without a
/// cut on the local chi2.
struct MeasurementSelectorCuts {
  /// bins in |eta| to specify variable selections
  std::vector<double> etaBins;
  /// Maximum local chi2 contribution.
  std::vector<double> chi2CutOff{std::numeric_limits<double>::max()};
  /// Maximum number of associated measurements on a single surface.
  std::vector<size_t> numMeasurementsCutOff{1};
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
  MeasurementSelector(Config config) : m_config(std::move(config)) {}

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
         bool& isOutlier, const Logger& logger) const {
    using Result = Result<std::pair<
        typename std::vector<typename traj_t::TrackStateProxy>::iterator,
        typename std::vector<typename traj_t::TrackStateProxy>::iterator>>;

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

    assert(!cuts->chi2CutOff.empty());
    const auto& chi2CutOff = cuts->chi2CutOff;
    auto maxChi2Cut = *std::max_element(chi2CutOff.begin(), chi2CutOff.end());
    double minChi2 = std::numeric_limits<double>::max();
    size_t minIndex = 0;
    auto trackStateIterEnd = candidates.end();
    {
      auto trackStateIter = candidates.begin();
      // Loop over all measurements to select the compatible measurements
      // Sort track states which do not satisfy the chi2 cut to the end.
      // When done trackStateIterEnd will point to the first element that
      // does not satisfy the chi2 cut.
      assert(trackStateIter != trackStateIterEnd);
      for (;;) {
        double chi2 = calculateChi2(
            // This abuses an incorrectly sized vector / matrix to access the
            // data pointer! This works (don't use the matrix as is!), but be
            // careful!
            trackStateIter
                ->template calibrated<
                    MultiTrajectoryTraits::MeasurementSizeMax>()
                .data(),
            trackStateIter
                ->template calibratedCovariance<
                    MultiTrajectoryTraits::MeasurementSizeMax>()
                .data(),
            trackStateIter->predicted(), trackStateIter->predictedCovariance(),
            trackStateIter->projector(), trackStateIter->calibratedSize());

        trackStateIter->chi2() = chi2;

        // only consider track states which pass the chi2 cut
        if (chi2 >= maxChi2Cut ||
            chi2 >= VariableCut<traj_t>(*trackStateIter, cuts, chi2CutOff,
                                        logger)) {
          --trackStateIterEnd;
          // still check whether the element has the smallest chi2 in case an
          // outlier is returned.
          if (chi2 < minChi2) {
            minChi2 = chi2;
            // the current element will be swapped with the last unchecked
            // element if they are different
            minIndex = std::distance(candidates.begin(), trackStateIterEnd);
          }

          if (trackStateIter == trackStateIterEnd) {
            break;
          } else {
            // swap rejected element with last element in list
            std::swap(*trackStateIter, *trackStateIterEnd);
          }
        } else {
          // Search for the measurement with the min chi2
          // @if there is a track state which passes the cut-off there is
          // no need to remember the track state with the smallest chi2.
          ++trackStateIter;
          if (trackStateIter == trackStateIterEnd) {
            break;
          }
        }
      }
    }

    // If there are no measurements below the chi2 cut off, return the
    // measurement with the best chi2 and tag it as an outlier
    if (candidates.begin() == trackStateIterEnd) {
      const auto bestIt = std::next(candidates.begin(), minIndex);
      ACTS_VERBOSE(
          "No measurement candidate. Return an outlier measurement chi2="
          << bestIt->chi2());
      isOutlier = true;
      // return single item range, no sorting necessary
      return Result::success(std::pair{bestIt, std::next(bestIt, 1)});
    }

    std::sort(candidates.begin(), trackStateIterEnd,
              [](const auto& tsa, const auto& tsb) {
                return tsa.chi2() < tsb.chi2();
              });

    // use |eta| of best measurement to select numMeasurementsCut
    const auto numMeasurementsCut = VariableCut<traj_t>(
        *candidates.begin(), cuts, cuts->numMeasurementsCutOff, logger);

    if (static_cast<std::size_t>(std::distance(
            candidates.begin(), trackStateIterEnd)) > numMeasurementsCut &&
        numMeasurementsCut > 0) {
      trackStateIterEnd = std::next(candidates.begin(), numMeasurementsCut);
    }

    ACTS_VERBOSE("Number of selected measurements: "
                 << std::distance(candidates.begin(), trackStateIterEnd)
                 << ", max: " << numMeasurementsCut);

    isOutlier = false;
    return std::pair{candidates.begin(), trackStateIterEnd};
  }

 private:
  template <typename traj_t, typename cut_value_t>
  static cut_value_t VariableCut(
      const typename traj_t::TrackStateProxy& trackState,
      const Acts::MeasurementSelector::Config::Iterator selector,
      const std::vector<cut_value_t>& cuts, const Logger& logger) {
    const auto& etaBins = selector->etaBins;
    if (etaBins.empty()) {
      return cuts[0];  // shortcut if no etaBins
    }
    const auto eta = std::atanh(std::cos(trackState.predicted()[eBoundTheta]));
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

  double calculateChi2(
      double* fullCalibrated, double* fullCalibratedCovariance,
      TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                       false>::Parameters predicted,
      TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                       false>::Covariance predictedCovariance,
      TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                       false>::Projector projector,
      unsigned int calibratedSize) const;

  Config m_config;
};

}  // namespace Acts
