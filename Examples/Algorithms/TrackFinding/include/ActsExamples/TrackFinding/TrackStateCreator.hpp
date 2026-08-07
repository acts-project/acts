// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/TrackFinding/TrackStateCreatorBase.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <utility>

namespace ActsExamples {

/// Cuts applied when extending a track onto a surface.
struct TrackStateSelectionCuts {
  /// Maximum local chi2 contribution to accept a measurement
  double chi2CutOffMeasurement = 15;
  /// Maximum local chi2 contribution to accept an outlier. Below the
  /// measurement cut this has no effect, since a compatible measurement is
  /// always preferred over an outlier.
  double chi2CutOffOutlier = 25;
  /// Maximum number of branches created on a single surface
  std::size_t maxBranchesPerSurface = 1;
};

/// Creates track states directly from the measurement container, i.e. without
/// pushing every measurement on a surface through the track EDM first.
///
/// This carries the selection, which is the local chi2 with respect to the
/// predicted parameters cut on @ref TrackStateSelectionCuts. The derived class
/// says which measurements are on the current surface via `measurementRange`
/// and hands out the calibrated form of one via `calibrate`; those are the
/// only places the concrete measurement container shows up.
///
/// @tparam derived_t the concrete creator, see CRTP
template <typename derived_t>
class TrackStateCreatorBase
    : public Acts::TrackStateCreatorBase<derived_t, TrackContainer> {
 public:
  /// Type alias for the CRTP base
  using Base = Acts::TrackStateCreatorBase<derived_t, TrackContainer>;
  /// Type alias for the per surface state object
  using State = typename Base::State;
  /// Type alias for the track state proxy
  using TrackStateProxy = typename Base::TrackStateProxy;

  /// The selection cuts
  TrackStateSelectionCuts cuts;

  /// The uncalibrated source link of a measurement.
  ///
  /// @param state the state of the current surface
  /// @param measurement the measurement
  /// @return the source link
  Acts::SourceLink sourceLink(const State& state,
                              const IndexSourceLink& measurement) const {
    static_cast<void>(state);

    return Acts::SourceLink{measurement};
  }

  /// Whether @ref preselectMeasurement should be consulted on this surface.
  ///
  /// @param state the state of the current surface
  /// @return false, i.e. no preselection
  bool hasMeasurementPreselection(const State& state) const {
    static_cast<void>(state);

    return false;
  }

  /// Cheap filter applied before any measurement is calibrated.
  ///
  /// @note If this rejects every measurement on the surface it is ignored and
  ///       all measurements are considered.
  ///
  /// @param state the state of the current surface
  /// @param measurement the measurement to check
  /// @return true, i.e. accept everything
  bool preselectMeasurement(const State& state,
                            const IndexSourceLink& measurement) const {
    static_cast<void>(state);
    static_cast<void>(measurement);

    return true;
  }

  /// Select the measurements compatible with the predicted parameters.
  ///
  /// If none of them satisfies the chi2 measurement cut, the best one is
  /// selected as an outlier if it satisfies the outlier cut. Otherwise nothing
  /// is selected and the caller creates a hole. Measurements and outliers are
  /// never mixed.
  ///
  /// @tparam range_t the measurement range type
  /// @tparam candidates_t the candidate buffer type
  ///
  /// @param state the state of the current surface
  /// @param range the measurements on the current surface
  /// @param candidates the buffer to append the selected measurements to
  ///
  /// @return always success
  template <typename range_t, typename candidates_t>
  Acts::Result<void> selectMeasurements(const State& state, range_t&& range,
                                        candidates_t& candidates) const {
    const Acts::Logger& logger = *state.logger;

    if (derived().hasMeasurementPreselection(state)) {
      for (const IndexSourceLink& measurement : range) {
        if (derived().preselectMeasurement(state, measurement)) {
          candidates.push_back({measurement, 0., false});
        }
      }
    }
    // If the preselection is inactive, or rejected every measurement on this
    // surface, fall back to considering all of them.
    if (candidates.empty()) {
      for (const IndexSourceLink& measurement : range) {
        candidates.push_back({measurement, 0., false});
      }
    }

    // Compute the chi2 of every candidate and sort those which do not satisfy
    // the chi2 cut to the back. When done `passedCandidates` is the number of
    // candidates at the front which do satisfy it.
    double minChi2 = std::numeric_limits<double>::max();
    std::size_t minIndex = std::numeric_limits<std::size_t>::max();
    std::size_t passedCandidates = 0ul;
    for (std::size_t i = 0ul; i < candidates.size(); ++i) {
      const double chi2 = derived().calculatePredictedChi2(
          state, derived().calibrate(state, candidates[i].measurement));
      candidates[i].chi2 = chi2;

      if (chi2 < minChi2) {
        minChi2 = chi2;
        minIndex = i;
      }

      if (chi2 >= cuts.chi2CutOffMeasurement) {
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
      if (minChi2 < cuts.chi2CutOffOutlier) {
        ACTS_VERBOSE(
            "No measurement candidate. Create an outlier measurement chi2="
            << minChi2);
        candidates[minIndex].isOutlier = true;
        keepOnly(candidates, minIndex);
      } else {
        ACTS_VERBOSE("No measurement candidate. Create none chi2=" << minChi2);
        candidates.clear();
      }
      return Acts::Result<void>::success();
    }

    if (passedCandidates == 1ul) {
      // single candidate, no sorting necessary
      ACTS_VERBOSE("Creating only 1 track state chi2=" << minChi2);
      keepOnly(candidates, minIndex);
      return Acts::Result<void>::success();
    }

    std::sort(candidates.begin(), candidates.begin() + passedCandidates,
              [](const auto& a, const auto& b) { return a.chi2 < b.chi2; });

    ACTS_VERBOSE("Number of selected measurements: "
                 << passedCandidates
                 << ", max: " << cuts.maxBranchesPerSurface);

    candidates.erase(candidates.begin() +
                         std::min(cuts.maxBranchesPerSurface, passedCandidates),
                     candidates.end());
    return Acts::Result<void>::success();
  }

 private:
  const derived_t& derived() const {
    return static_cast<const derived_t&>(*this);
  }

  /// Reduce the candidates to the single one at @p index.
  template <typename candidates_t>
  static void keepOnly(candidates_t& candidates, std::size_t index) {
    if (index != 0ul) {
      candidates[0] = candidates[index];
    }
    candidates.erase(candidates.begin() + 1, candidates.end());
  }
};

}  // namespace ActsExamples
