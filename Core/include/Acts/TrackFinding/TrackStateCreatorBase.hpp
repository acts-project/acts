// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/MeasurementConcept.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterExtensions.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cassert>
#include <cstddef>
#include <optional>
#include <ranges>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @addtogroup track_finding
/// @{

/// Write a calibrated measurement into a track state.
///
/// @note This passes exactly `measurement.size()` subspace indices. The
///       projector is serialized including its unused slots, so passing a
///       padded range here would not be equivalent.
///
/// @tparam measurement_t the calibrated measurement type
/// @tparam track_state_proxy_t the track state proxy type
///
/// @param measurement the calibrated measurement
/// @param trackState the track state to write to
template <MeasurementConcept measurement_t, typename track_state_proxy_t>
void fillTrackStateCalibrated(const measurement_t& measurement,
                              track_state_proxy_t trackState) {
  if constexpr (StaticMeasurementConcept<measurement_t>) {
    trackState.allocateCalibrated(measurement.parameters().eval(),
                                  measurement.covariance().eval());
  } else {
    visit_measurement(
        measurement.parameters(), measurement.covariance(), measurement.size(),
        [&trackState](const auto& parameters, const auto& covariance) {
          trackState.allocateCalibrated(parameters.eval(), covariance.eval());
        });
  }
  trackState.setProjectorSubspaceIndices(measurement.subspaceIndices());
}

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

/// A measurement which was selected to become a track state on a surface.
///
/// @tparam measurement_t whatever the creator calls a measurement
template <typename measurement_t>
struct TrackStateCandidate {
  /// The selected measurement
  measurement_t measurement;
  /// Its compatibility with the predicted parameters, written to the track
  /// state as is. Zero if the selection does not compute one.
  double chi2{};
  /// Whether the track state is to be flagged as an outlier
  bool isOutlier{false};
};

/// CRTP base for the creation of track states in the combinatorial Kalman
/// filter.
///
/// The derived class selects on whatever it calls a measurement, and only the
/// selected ones are materialized as track states. Measurements never have to
/// pass through the track EDM, which is what the classic creator had to do to
/// run its selection.
///
/// This class only provides the mechanism, i.e. correctly populating the track
/// EDM. It does not impose a selection; that is what `selectMeasurements` is
/// for.
///
/// The derived class has to provide
/// - `measurementRange(state)`: the measurements associated to `state.surface`
/// - `sourceLink(state, measurement)`: the uncalibrated source link
/// - `calibrate(state, measurement)`: a @ref MeasurementConcept value
///
/// and may override any of the defaulted customization points below.
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
  /// Type alias for the result of the track state creation
  using TrackStatesResult = Result<CkfTypes::BranchVector<TrackIndexType>>;

  /// Expected maximum number of measurements considered on a single surface.
  /// Exceeding this only costs a dynamic allocation of the candidate buffer.
  static constexpr std::size_t kMaxCandidatesPerSurface = 16;

  /// The buffer a selection fills with the measurements it selected
  /// @tparam measurement_t whatever the creator calls a measurement
  template <typename measurement_t>
  using Candidates =
      boost::container::small_vector<TrackStateCandidate<measurement_t>,
                                     kMaxCandidatesPerSurface>;

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
  /// @return the indices of the newly created track states, in the order the
  ///         selection produced them, or an error. An empty list means that no
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
    // Note this has to happen before the selection is consulted. Resolving
    // cuts can fail for surfaces which are not covered by the configuration,
    // and a surface without measurements must not turn that into an error.
    if (std::ranges::empty(range)) {
      ACTS_VERBOSE("No measurements on surface " << surface.geometryId());
      return trackStates;
    }

    using Measurement = std::ranges::range_value_t<decltype(range)>;
    Candidates<Measurement> candidates;

    if (Result<void> selection =
            derived().selectMeasurements(state, range, candidates);
        !selection.ok()) {
      return selection.error();
    }

    trackStates.reserve(candidates.size());
    for (const TrackStateCandidate<Measurement>& candidate : candidates) {
      trackStates.push_back(derived().createTrackState(state, candidate));
    }
    return trackStates;
  }

  /// Pick the measurements which are to become track states on this surface.
  ///
  /// This is the main customization point. Implementations append to
  /// @p candidates in the order the track states are to be created in, which
  /// for the combinatorial Kalman filter should be the most compatible one
  /// first. An empty result means that the caller creates a hole.
  ///
  /// @note This is only called when the surface has measurements.
  ///
  /// @tparam range_t the measurement range type
  /// @tparam candidates_t the candidate buffer type
  ///
  /// @param state the state of the current surface
  /// @param range the measurements on the current surface
  /// @param candidates the buffer to append the selected measurements to
  ///
  /// @return success, i.e. the default takes all measurements
  template <typename range_t, typename candidates_t>
  Result<void> selectMeasurements(const State& state, range_t&& range,
                                  candidates_t& candidates) const {
    static_cast<void>(state);

    for (const auto& measurement : range) {
      candidates.push_back({measurement, 0., false});
    }
    return Result<void>::success();
  }

  /// Local chi2 contribution of a calibrated measurement with respect to the
  /// predicted bound parameters on this surface.
  ///
  /// This is what a selection is expected to cut on. It is a customization
  /// point so that a derived class which knows more about its measurements can
  /// compute the compatibility differently, e.g. accounting for a non Gaussian
  /// measurement model.
  ///
  /// The default reads the measurement through @ref MeasurementConcept, which
  /// hands out Eigen maps into the original storage, so nothing is copied. If
  /// the measurement dimension is known at compile time the Eigen operations
  /// are statically sized; otherwise the runtime dimension is dispatched onto
  /// a statically sized instantiation first. Either way Eigen never allocates.
  ///
  /// @tparam measurement_t the calibrated measurement type
  ///
  /// @param state the state of the current surface
  /// @param measurement the calibrated measurement
  ///
  /// @return the local chi2 contribution
  template <MeasurementConcept measurement_t>
  double calculatePredictedChi2(const State& state,
                                const measurement_t& measurement) const {
    if constexpr (StaticMeasurementConcept<measurement_t>) {
      return predictedChi2(measurement.parameters(), measurement.covariance(),
                           measurement.subspaceIndices(), state.predicted(),
                           state.predictedCovariance());
    } else {
      return visit_measurement(
          measurement.parameters(), measurement.covariance(),
          measurement.size(),
          [&measurement, &state](const auto& calibrated,
                                 const auto& calibratedCovariance) {
            return predictedChi2(
                calibrated, calibratedCovariance, measurement.subspaceIndices(),
                state.predicted(), state.predictedCovariance());
          });
    }
  }

  /// Write the calibrated measurement into a selected track state.
  ///
  /// Override this together with the selection to calibrate lazily, i.e. to
  /// select on something cheaper than the full calibration.
  ///
  /// @tparam measurement_t whatever the creator calls a measurement
  ///
  /// @param state the state of the current surface
  /// @param measurement the selected measurement
  /// @param trackState the track state to fill
  template <typename measurement_t>
  void fillCalibrated(const State& state, const measurement_t& measurement,
                      TrackStateProxy trackState) const {
    fillTrackStateCalibrated(derived().calibrate(state, measurement),
                             trackState);
  }

  /// Create and fill a single track state for a selected measurement.
  ///
  /// @tparam measurement_t whatever the creator calls a measurement
  ///
  /// @param state the state of the current surface
  /// @param candidate the selected measurement
  ///
  /// @return the index of the new track state
  template <typename measurement_t>
  TrackIndexType createTrackState(
      State& state, const TrackStateCandidate<measurement_t>& candidate) const {
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
    trackState.chi2() = candidate.chi2;
    trackState.setReferenceSurface(
        boundParams.referenceSurface().getSharedPtr());
    trackState.setUncalibratedSourceLink(
        derived().sourceLink(state, candidate.measurement));
    derived().fillCalibrated(state, candidate.measurement, trackState);

    auto typeFlags = trackState.typeFlags();
    typeFlags.setHasParameters();
    typeFlags.setHasMeasurement();
    if (trackState.referenceSurface().hasMaterial()) {
      typeFlags.setHasMaterial();
    }
    if (candidate.isOutlier) {
      typeFlags.setIsOutlier();
    }

    return trackState.index();
  }

 protected:
  /// Access to the concrete creator
  /// @return the concrete creator
  const derived_t& derived() const {
    return static_cast<const derived_t&>(*this);
  }

  /// The Eigen math behind @ref calculatePredictedChi2, with the measurement
  /// dimension known at compile time.
  ///
  /// Exposed so that an override of @ref calculatePredictedChi2 can still
  /// reuse it after unpacking its measurement differently.
  ///
  /// @tparam calibrated_t the measured values type
  /// @tparam calibrated_covariance_t the measurement covariance type
  /// @tparam index_range_t the subspace index range type
  ///
  /// @param calibrated the measured values
  /// @param calibratedCovariance the covariance of the measured values
  /// @param subspaceIndices the bound indices the measurement constrains,
  ///        exactly as many as the measurement dimension
  /// @param predicted the predicted bound parameters
  /// @param predictedCovariance the predicted bound parameters covariance
  ///
  /// @return the local chi2 contribution
  template <typename calibrated_t, typename calibrated_covariance_t,
            std::ranges::sized_range index_range_t>
  static double predictedChi2(
      const Eigen::MatrixBase<calibrated_t>& calibrated,
      const Eigen::MatrixBase<calibrated_covariance_t>& calibratedCovariance,
      const index_range_t& subspaceIndices,
      Eigen::Ref<const BoundVector> predicted,
      Eigen::Ref<const BoundMatrix> predictedCovariance) {
    constexpr int kMeasurementSize = calibrated_t::RowsAtCompileTime;
    static_assert(kMeasurementSize != Eigen::Dynamic,
                  "Measurement dimension must be known at compile time");

    const FixedBoundSubspaceHelper<kMeasurementSize> subspaceHelper(
        subspaceIndices);

    // Get the residuals
    Vector<kMeasurementSize> res =
        calibrated - subspaceHelper.projectVector(predicted);

    // Get the chi2
    return (res.transpose() *
            (calibratedCovariance +
             subspaceHelper.projectMatrix(predictedCovariance))
                .inverse() *
            res)
        .eval()(0, 0);
  }
};

/// @}

}  // namespace Acts
