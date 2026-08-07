// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/CalibratedBoundMeasurement.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterError.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterExtensions.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <optional>
#include <ranges>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @addtogroup track_finding
/// @{

/// Everything a @ref TrackStateCreatorBase customization point needs to know
/// about the surface it is currently working on.
///
/// This deliberately does not depend on the concrete creator type, so that
/// derived classes can implement their customization points as plain,
/// non-template member functions taking a `const State&`.
///
/// @tparam track_container_t the track container used by the track finder
template <typename track_container_t>
struct TrackStateCreatorState {
  /// Type alias for the track state container backend
  using TrackStateContainerBackend =
      typename track_container_t::TrackStateContainerBackend;
  /// Type alias for the track state proxy
  using TrackStateProxy = typename track_container_t::TrackStateProxy;

  /// The current geometry context
  const GeometryContext* geoContext{};
  /// The current calibration context
  const CalibrationContext* calibrationContext{};
  /// The surface the trajectory is being extended onto
  const Surface* surface{};
  /// The predicted bound state on @ref surface
  const CkfTypes::BoundState* boundState{};
  /// The tip of the trajectory which is being extended
  TrackIndexType prevTip{kTrackIndexInvalid};
  /// The trajectory new track states are added to
  TrackStateContainerBackend* trajectory{};
  /// A logger for messages
  const Logger* logger{};

  /// The selection cuts for this surface.
  /// @note only valid from `computeChi2` onwards
  MeasurementSelector::Cuts cuts{};

  /// The track state on this surface which owns the predicted parameters and
  /// the jacobian. All further track states on the same surface share them.
  std::optional<TrackStateProxy> firstTrackState{};

  /// The predicted bound parameters on the surface
  /// @return the bound track parameters
  const BoundTrackParameters& boundParameters() const {
    return std::get<0>(*boundState);
  }
  /// The predicted parameter vector
  /// @return the predicted parameters
  const BoundVector& predicted() const {
    return boundParameters().parameters();
  }
  /// The predicted parameter covariance
  /// @return the predicted covariance
  const BoundMatrix& predictedCovariance() const {
    assert(boundParameters().covariance().has_value() &&
           "Track state creation requires a predicted covariance");
    return *boundParameters().covariance();
  }
  /// The jacobian from the previous surface to this one
  /// @return the jacobian
  const BoundMatrix& jacobian() const { return std::get<1>(*boundState); }
  /// The path length accumulated up to this surface
  /// @return the path length
  double pathLength() const { return std::get<2>(*boundState); }
};

/// CRTP base for the creation of track states in the combinatorial Kalman
/// filter.
///
/// In contrast to @ref TrackStateCreator, which creates a temporary track
/// state for every measurement on a surface just to select a few of them, this
/// computes the compatibility chi2 directly from whatever the derived class
/// calls a measurement, and only materializes track states for the selected
/// ones. Measurements never have to pass through the track EDM.
///
/// The derived class has to provide
/// - `measurementRange(state)`: the measurements associated to `state.surface`
/// - `sourceLink(state, measurement)`: the uncalibrated source link
/// - `calibrate(state, measurement)`: an @ref CalibratedBoundMeasurement
///
/// and may override any of the defaulted customization points below. The
/// defaults reproduce the selection implemented by @ref MeasurementSelector.
///
/// @note This is expected to be connected to
///       @ref CombinatorialKalmanFilterExtensions::trackStateCreator.
/// @note The API is not stable yet and is expected to grow further
///       customization points, in particular for hole creation.
///
/// @tparam derived_t the concrete creator, see CRTP
/// @tparam track_container_t the track container used by the track finder
template <typename derived_t, typename track_container_t>
class TrackStateCreatorBase {
 public:
  /// Type alias for the track state container backend
  using TrackStateContainerBackend =
      typename track_container_t::TrackStateContainerBackend;
  /// Type alias for the track state proxy
  using TrackStateProxy = typename track_container_t::TrackStateProxy;
  /// Type alias for the per surface state object
  using State = TrackStateCreatorState<track_container_t>;
  /// Type alias for the selection cuts
  using Cuts = MeasurementSelector::Cuts;
  /// Type alias for the result of the track state creation
  using TrackStatesResult = Result<CkfTypes::BranchVector<TrackIndexType>>;

  /// Extend the trajectory onto the given surface.
  ///
  /// @param gctx The current geometry context
  /// @param cctx The current calibration context
  /// @param surface The surface the trajectory is extended onto
  /// @param boundState The predicted bound state on the given surface
  /// @param prevTip The tip of the trajectory which is to be extended
  /// @param trajectory The trajectory to be extended
  /// @param logger A logger for messages
  ///
  /// @return the indices of the newly created track states, ordered by
  ///         ascending chi2, or an error. An empty list means that no
  ///         compatible measurement was found and the caller is expected to
  ///         create a hole.
  TrackStatesResult createTrackStates(const GeometryContext& gctx,
                                      const CalibrationContext& cctx,
                                      const Surface& surface,
                                      const CkfTypes::BoundState& boundState,
                                      TrackIndexType prevTip,
                                      TrackStateContainerBackend& trajectory,
                                      const Logger& logger) const {
    State state{};
    state.geoContext = &gctx;
    state.calibrationContext = &cctx;
    state.surface = &surface;
    state.boundState = &boundState;
    state.prevTip = prevTip;
    state.trajectory = &trajectory;
    state.logger = &logger;

    CkfTypes::BranchVector<TrackIndexType> trackStates;

    auto range = derived().measurementRange(state);
    // Note this has to happen before the cuts are resolved. Resolving cuts can
    // fail for surfaces which are not covered by the configuration, and a
    // surface without measurements must not turn that into an error.
    if (std::ranges::empty(range)) {
      ACTS_VERBOSE("No measurements on surface " << surface.geometryId());
      return trackStates;
    }

    Result<Cuts> cutsResult = derived().selectionCuts(state);
    if (!cutsResult.ok()) {
      ACTS_DEBUG("Could not resolve selection cuts for surface "
                 << surface.geometryId() << ": "
                 << cutsResult.error().message());
      return cutsResult.error();
    }
    state.cuts = *cutsResult;
    if (state.cuts.numMeasurements == 0ul) {
      return CombinatorialKalmanFilterError::MeasurementSelectionFailed;
    }

    using Measurement = std::ranges::range_value_t<decltype(range)>;
    struct Candidate {
      Measurement measurement;
      double chi2{};
    };
    boost::container::small_vector<Candidate, s_maxCandidatesPerSurface>
        candidates;

    if (derived().hasMeasurementPreselection(state)) {
      for (const Measurement& measurement : range) {
        if (derived().preselectMeasurement(state, measurement)) {
          candidates.push_back({measurement, 0.});
        }
      }
    }
    // If the preselection is inactive, or rejected every measurement on this
    // surface, fall back to considering all of them.
    if (candidates.empty()) {
      for (const Measurement& measurement : range) {
        candidates.push_back({measurement, 0.});
      }
    }

    // Compute the chi2 of every candidate and sort those which do not satisfy
    // the chi2 cut to the back. When done `passedCandidates` is the number of
    // candidates at the front which do satisfy it.
    double minChi2 = std::numeric_limits<double>::max();
    std::size_t minIndex = std::numeric_limits<std::size_t>::max();
    std::size_t passedCandidates = 0ul;
    for (std::size_t i = 0ul; i < candidates.size(); ++i) {
      const double chi2 =
          derived().computeChi2(state, candidates[i].measurement);
      candidates[i].chi2 = chi2;

      if (chi2 < minChi2) {
        minChi2 = chi2;
        minIndex = i;
      }

      if (chi2 >= state.cuts.chi2Measurement) {
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
      if (minChi2 < state.cuts.chi2Outlier) {
        ACTS_VERBOSE(
            "No measurement candidate. Create an outlier measurement chi2="
            << minChi2);
        trackStates.push_back(createTrackState(state,
                                               candidates[minIndex].measurement,
                                               minChi2, /*isOutlier=*/true));
      } else {
        ACTS_VERBOSE("No measurement candidate. Create none chi2=" << minChi2);
      }
      return trackStates;
    }

    if (passedCandidates == 1ul) {
      // single candidate, no sorting necessary
      ACTS_VERBOSE("Creating only 1 track state chi2=" << minChi2);
      trackStates.push_back(createTrackState(state,
                                             candidates[minIndex].measurement,
                                             minChi2, /*isOutlier=*/false));
      return trackStates;
    }

    std::sort(
        candidates.begin(), candidates.begin() + passedCandidates,
        [](const Candidate& a, const Candidate& b) { return a.chi2 < b.chi2; });

    const std::size_t selected =
        std::min(state.cuts.numMeasurements, passedCandidates);
    ACTS_VERBOSE("Number of selected measurements: "
                 << passedCandidates
                 << ", max: " << state.cuts.numMeasurements);

    trackStates.reserve(selected);
    for (std::size_t i = 0ul; i < selected; ++i) {
      trackStates.push_back(createTrackState(state, candidates[i].measurement,
                                             candidates[i].chi2,
                                             /*isOutlier=*/false));
    }
    return trackStates;
  }

 protected:
  /// Expected maximum number of measurements on a single surface. Exceeding
  /// this only costs a dynamic allocation of the candidate buffer.
  static constexpr std::size_t s_maxCandidatesPerSurface = 16;

  /// The selection cuts to apply on this surface.
  ///
  /// @note This is only called when the surface has measurements.
  ///
  /// @param state the state of the current surface
  /// @return fully permissive cuts
  Result<Cuts> selectionCuts(const State& state) const {
    static_cast<void>(state);

    return Cuts{std::numeric_limits<std::size_t>::max(),
                std::numeric_limits<double>::infinity(),
                -std::numeric_limits<double>::infinity()};
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
  template <typename measurement_t>
  bool preselectMeasurement(const State& state,
                            const measurement_t& measurement) const {
    static_cast<void>(state);
    static_cast<void>(measurement);

    return true;
  }

  /// The chi2 of a measurement with respect to the predicted parameters.
  ///
  /// Override this together with @ref fillCalibrated to select on something
  /// cheaper than the full calibration.
  ///
  /// @param state the state of the current surface
  /// @param measurement the measurement to evaluate
  /// @return the local chi2 contribution
  template <typename measurement_t>
  double computeChi2(const State& state,
                     const measurement_t& measurement) const {
    const CalibratedBoundMeasurement calibrated =
        derived().calibrate(state, measurement);
    return calculatePredictedChi2(
        calibrated.parametersData(), calibrated.covarianceData(),
        state.predicted(), state.predictedCovariance(),
        calibrated.paddedSubspaceIndices(), calibrated.size());
  }

  /// Write the calibrated measurement into a selected track state.
  ///
  /// @param state the state of the current surface
  /// @param measurement the selected measurement
  /// @param trackState the track state to fill
  template <typename measurement_t>
  void fillCalibrated(const State& state, const measurement_t& measurement,
                      TrackStateProxy trackState) const {
    derived().calibrate(state, measurement).fill(trackState);
  }

  /// Create and fill a single track state for a selected measurement.
  ///
  /// @param state the state of the current surface
  /// @param measurement the selected measurement
  /// @param chi2 the chi2 of the measurement
  /// @param isOutlier whether the measurement is to be flagged as an outlier
  /// @return the index of the new track state
  template <typename measurement_t>
  TrackIndexType createTrackState(State& state,
                                  const measurement_t& measurement, double chi2,
                                  bool isOutlier) const {
    using PM = TrackStatePropMask;

    const bool isFirst = !state.firstTrackState.has_value();

    // The filtered parameters are deliberately not allocated here. Outliers
    // never get separate filtered parameters, and for measurements the caller
    // allocates them right before the Kalman update writes them, so that
    // nobody can observe allocated but uninitialized filtered parameters via
    // `parameters()`.
    PM mask = PM::Calibrated;
    if (isFirst) {
      mask |= PM::Predicted | PM::Jacobian;
    }

    TrackStateProxy trackState =
        state.trajectory->makeTrackState(mask, state.prevTip);

    const BoundTrackParameters& boundParams = state.boundParameters();
    if (isFirst) {
      trackState.predicted() = boundParams.parameters();
      if (boundParams.covariance().has_value()) {
        trackState.predictedCovariance() = *boundParams.covariance();
      }
      trackState.jacobian() = state.jacobian();
      state.firstTrackState = trackState;
    } else {
      // subsequent track states share these with the first one
      trackState.shareFrom(*state.firstTrackState, PM::Predicted);
      trackState.shareFrom(*state.firstTrackState, PM::Jacobian);
    }

    trackState.pathLength() = state.pathLength();
    trackState.chi2() = chi2;
    trackState.setReferenceSurface(
        boundParams.referenceSurface().getSharedPtr());
    trackState.setUncalibratedSourceLink(
        derived().sourceLink(state, measurement));
    derived().fillCalibrated(state, measurement, trackState);

    auto typeFlags = trackState.typeFlags();
    typeFlags.setHasParameters();
    typeFlags.setHasMeasurement();
    if (trackState.referenceSurface().hasMaterial()) {
      typeFlags.setHasMaterial();
    }
    if (isOutlier) {
      typeFlags.setIsOutlier();
    }

    return trackState.index();
  }

 private:
  const derived_t& derived() const {
    return static_cast<const derived_t&>(*this);
  }
};

/// A @ref TrackStateCreatorBase which sources its selection cuts from an
/// @ref MeasurementSelector, i.e. reproduces the geometry and eta binned cut
/// lookup of the classic track state creator.
///
/// @tparam derived_t the concrete creator, see CRTP
/// @tparam track_container_t the track container used by the track finder
template <typename derived_t, typename track_container_t>
class MeasurementSelectorTrackStateCreatorBase
    : public TrackStateCreatorBase<derived_t, track_container_t> {
 public:
  /// Type alias for the CRTP base
  using Base = TrackStateCreatorBase<derived_t, track_container_t>;
  /// Type alias for the per surface state object
  using State = typename Base::State;
  /// Type alias for the selection cuts
  using Cuts = typename Base::Cuts;

  /// The selector providing the geometry and eta binned cuts
  MeasurementSelector measurementSelector;

  /// Resolve the cuts for the current surface.
  ///
  /// @param state the state of the current surface
  /// @return the cuts or an error if the surface is not covered
  Result<Cuts> selectionCuts(const State& state) const {
    return measurementSelector.getCuts(state.surface->geometryId(),
                                       state.predicted()[eBoundTheta]);
  }
};

/// @}

}  // namespace Acts
