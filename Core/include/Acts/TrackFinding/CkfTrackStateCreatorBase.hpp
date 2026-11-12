// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
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
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <ranges>
#include <type_traits>
#include <utility>

#include <boost/container/small_vector.hpp>

namespace Acts {

namespace detail {

/// The parameter vector type of a measurement like type
/// @tparam T the measurement like type
template <typename T>
using MeasurementParametersType =
    std::decay_t<decltype(std::declval<const T&>().parameters())>;

/// The covariance matrix type of a measurement like type
/// @tparam T the measurement like type
template <typename T>
using MeasurementCovarianceType =
    std::decay_t<decltype(std::declval<const T&>().covariance())>;

/// Requirements on a calibrated measurement in a subspace of the bound
/// parameter space.
///
/// The track state creator uses such a value as is: it computes the chi2 from
/// it and writes it to the calibrated slot of the track state. The concept only
/// checks the shape, so it cannot tell a calibrated measurement from an
/// uncalibrated one. Calibration is the responsibility of `calibrate` in the
/// derived creator.
///
/// The dimension may be known at compile time or only at runtime, see
/// @c StaticCalibratedMeasurementConcept.
///
/// @tparam T the measurement like type
template <typename T>
concept CalibratedMeasurementConcept = requires(const T& measurement) {
  /// The measurement dimension
  { measurement.size() } -> std::integral;
  /// The bound parameters the measurement constrains, exactly `size()` of
  /// them
  { measurement.subspaceIndices() } -> std::ranges::sized_range;
  /// The measured values
  { measurement.parameters() };
  /// The covariance of the measured values
  { measurement.covariance() };

  requires std::convertible_to<
      std::ranges::range_value_t<decltype(measurement.subspaceIndices())>,
      std::uint8_t>;
  requires detail::MeasurementParametersType<T>::ColsAtCompileTime == 1;
  requires detail::MeasurementParametersType<T>::RowsAtCompileTime ==
               detail::MeasurementCovarianceType<T>::RowsAtCompileTime;
};

/// The dimension of a measurement if it is known at compile time,
/// `Eigen::Dynamic` otherwise.
///
/// @tparam measurement_t the measurement type
template <CalibratedMeasurementConcept measurement_t>
constexpr int CalibratedMeasurementSizeAtCompileTime =
    detail::MeasurementParametersType<measurement_t>::RowsAtCompileTime;

/// A calibrated measurement whose dimension is known at compile time.
///
/// Algorithms can use this to skip the runtime dispatch over the measurement
/// dimension.
///
/// @tparam measurement_t the measurement type
template <typename measurement_t>
concept StaticCalibratedMeasurementConcept =
    CalibratedMeasurementConcept<measurement_t> &&
    (CalibratedMeasurementSizeAtCompileTime<measurement_t> != Eigen::Dynamic);

/// Transport a bound state onto another surface along a straight line.
///
/// This ignores the magnetic field and the material between the two surfaces,
/// so it is only a good approximation for surfaces close to each other.
///
/// @param gctx the geometry context
/// @param boundState the bound state to transport, with a covariance
/// @param surface the target surface
/// @param propagationDirection the direction of the propagation with respect
///        to the momentum. The target must be ahead in this direction.
///
/// @return the bound state on @p surface. Its jacobian starts where the
///         jacobian of @p boundState starts, and its path length includes the
///         distance between the surfaces. An error if the surface is not ahead.
Result<CkfTypes::BoundState> transportBoundStateStraightLine(
    const GeometryContext& gctx, const CkfTypes::BoundState& boundState,
    const Surface& surface, Direction propagationDirection);

}  // namespace detail

/// @addtogroup track_finding
/// @{

/// Write a calibrated measurement into a track state.
///
/// @note The projector is serialized including its unused slots, so the
///       subspace indices have to be exactly `size()` and not padded.
///
/// @tparam measurement_t the calibrated measurement type
/// @tparam track_state_proxy_t the track state proxy type
///
/// @param measurement the calibrated measurement
/// @param trackState the track state to write to
template <detail::CalibratedMeasurementConcept measurement_t,
          typename track_state_proxy_t>
void fillTrackStateCalibrated(const measurement_t& measurement,
                              track_state_proxy_t trackState) {
  if constexpr (detail::StaticCalibratedMeasurementConcept<measurement_t>) {
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

/// Everything a @ref CkfTrackStateCreatorBase customization point needs to know
/// about the surface it is currently working on.
///
/// This is the input from the CKF, see @ref CkfTypes::TrackStateCreatorInput,
/// together with the scratch data of the creator. It does not depend on the
/// concrete creator type, so that customization points can be plain,
/// non-template member functions.
///
/// @tparam track_container_t the track container used by the track finder
template <typename track_container_t>
struct CkfTrackStateCreatorState
    : public CkfTypes::TrackStateCreatorInput<track_container_t> {
  /// Type alias for the input from the CKF
  using Input = CkfTypes::TrackStateCreatorInput<track_container_t>;
  /// Type alias for the track state proxy
  using TrackStateProxy = typename track_container_t::TrackStateProxy;

  /// Make the state for one surface from the input of the CKF
  /// @param input the input from the CKF
  explicit CkfTrackStateCreatorState(const Input& input) : Input(input) {}

  /// The track state on this surface which owns the predicted parameters and
  /// the jacobian. All further track states on the same surface share them.
  std::optional<TrackStateProxy> firstTrackState{};

  /// The predicted bound parameters on the surface
  /// @return the bound track parameters
  const BoundTrackParameters& boundParameters() const {
    return std::get<0>(*this->boundState);
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
  const BoundMatrix& jacobian() const { return std::get<1>(*this->boundState); }
  /// The path length accumulated up to this surface
  /// @return the path length
  double pathLength() const { return std::get<2>(*this->boundState); }
};

/// A measurement which was selected to become a track state on a surface.
///
/// @tparam measurement_t whatever the creator calls a measurement
template <typename measurement_t>
struct CkfTrackStateCandidate {
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
/// selected ones are materialized as track states, so measurements never have
/// to pass through the track EDM just to be selected on.
///
/// This class only provides the mechanism, i.e. correctly populating the track
/// EDM. The derived class has to provide
/// - `measurementRange(state)`: the measurements associated to `state.surface`
/// - `sourceLink(state, measurement)`: the uncalibrated source link
/// - `calibrate(state, measurement)`: the calibrated measurement, see
///   @c detail::CalibratedMeasurementConcept
///
/// and may override any of the defaulted customization points below.
///
/// @note This is expected to be connected to
///       @ref CombinatorialKalmanFilterExtensions::trackStateCreator.
/// @note The API is not stable yet and is expected to grow further
///       customization points.
///
/// @tparam derived_t the concrete creator, see CRTP
/// @tparam track_container_t the track container used by the track finder
template <typename derived_t, typename track_container_t>
class CkfTrackStateCreatorBase {
 public:
  /// Type alias for the track state container backend
  using TrackStateContainerBackend =
      typename track_container_t::TrackStateContainerBackend;
  /// Type alias for the track state proxy
  using TrackStateProxy = typename track_container_t::TrackStateProxy;
  /// Type alias for the per surface state object
  using State = CkfTrackStateCreatorState<track_container_t>;
  /// Type alias for the result of the track state creation
  using TrackStatesResult = Result<CkfTypes::BranchVector<TrackIndexType>>;

  /// Expected maximum number of measurements considered on a single surface.
  /// Exceeding this only costs a dynamic allocation of the candidate buffer.
  static constexpr std::size_t kMaxCandidatesPerSurface = 16;

  /// The buffer a selection fills with the measurements it selected
  /// @tparam measurement_t whatever the creator calls a measurement
  template <typename measurement_t>
  using Candidates =
      boost::container::small_vector<CkfTrackStateCandidate<measurement_t>,
                                     kMaxCandidatesPerSurface>;

  /// Type alias for the input from the CKF
  using Input = CkfTypes::TrackStateCreatorInput<track_container_t>;

  /// Extend the trajectory onto the given surface.
  ///
  /// This is the function to connect to
  /// @ref CombinatorialKalmanFilterExtensions::trackStateCreator.
  ///
  /// @param input the surface, the predicted bound state, the trajectory to
  ///        extend, and the contexts
  ///
  /// @return the indices of the newly created track states, in the order the
  ///         selection produced them, first on the current surface and then
  ///         on each of @ref neighbourSurfaces, or an error. If no selection
  ///         keeps a measurement, these are the track states from
  ///         @ref createHoleTrackStates.
  TrackStatesResult createTrackStates(const Input& input) const {
    assert(input.geoContext != nullptr && input.calibrationContext != nullptr &&
           input.surface != nullptr && input.boundState != nullptr &&
           input.trajectory != nullptr && input.logger != nullptr &&
           "Incomplete input for the track state creation");

    const Logger& logger = *input.logger;
    State state(input);

    CkfTypes::BranchVector<TrackIndexType> trackStates;

    if (Result<void> res = createMeasurementTrackStates(state, trackStates);
        !res.ok()) {
      return res.error();
    }

    for (const Surface* neighbour : derived().neighbourSurfaces(state)) {
      assert(neighbour != nullptr && neighbour != input.surface);

      Result<CkfTypes::BoundState> neighbourBoundState =
          derived().transportToSurface(state, *neighbour);
      if (!neighbourBoundState.ok()) {
        ACTS_VERBOSE("Skip neighbour surface " << neighbour->geometryId()
                                               << ": "
                                               << neighbourBoundState.error());
        continue;
      }

      Input neighbourInput = input;
      neighbourInput.surface = neighbour;
      neighbourInput.boundState = &*neighbourBoundState;
      State neighbourState(neighbourInput);

      if (Result<void> res =
              createMeasurementTrackStates(neighbourState, trackStates);
          !res.ok()) {
        return res.error();
      }
    }

    if (trackStates.empty()) {
      ACTS_VERBOSE("No measurement selected on surface "
                   << input.surface->geometryId());
      if (Result<void> res =
              derived().createHoleTrackStates(state, trackStates);
          !res.ok()) {
        return res.error();
      }
    }

    return trackStates;
  }

  /// The surfaces next to the current one which are also searched for
  /// measurements, e.g. overlapping modules.
  ///
  /// The creator selects on each of these surfaces as on the current one, and
  /// the track states on them come after the ones on the current surface. The
  /// creator makes a hole only if no surface has a track state.
  ///
  /// @note A surface must be ahead in the propagation direction, because the
  ///       CKF continues from there. @ref transportToSurface rejects one
  ///       which is not.
  /// @note If the CKF continues from a neighbour surface, it skips all
  ///       surfaces between the current surface and the neighbour, so these
  ///       get no track state, not even a hole.
  /// @note If the CKF continues from a neighbour surface, it skips the part of
  ///       the material of the current surface which comes after the update.
  ///       The default @ref transportToSurface also skips the part of the
  ///       material of the neighbour surface which comes before the update.
  ///
  /// @param state the state of the current surface
  ///
  /// @return a range of `const Surface*`, empty by default
  auto neighbourSurfaces(const State& state) const {
    static_cast<void>(state);

    return std::ranges::empty_view<const Surface*>{};
  }

  /// The predicted bound state on a neighbour surface.
  ///
  /// The default transports the bound state of the current surface along a
  /// straight line, see @ref detail::transportBoundStateStraightLine. It
  /// ignores the magnetic field and all material, so override this if the
  /// surfaces are far apart or carry a lot of material.
  ///
  /// @param state the state of the current surface
  /// @param surface the neighbour surface
  ///
  /// @return the bound state on @p surface, or an error to skip it
  Result<CkfTypes::BoundState> transportToSurface(
      const State& state, const Surface& surface) const {
    return detail::transportBoundStateStraightLine(*state.geoContext,
                                                   *state.boundState, surface,
                                                   state.propagationDirection);
  }

  /// Create the track states on a surface on which the selection kept no
  /// measurement, or which has no measurements at all, and on whose neighbour
  /// surfaces the selection kept none either.
  ///
  /// The default creates one track state without a measurement. It is a hole
  /// if @ref isMeasurementExpected, and otherwise a state flagged with no
  /// expected hit, which does not count as a hole.
  ///
  /// The CKF counts a track state flagged as a hole in the holes of the
  /// branch. Each track state makes its own branch. If this adds no track
  /// state, the CKF creates the hole itself.
  ///
  /// @param state the state of the current surface
  /// @param trackStates the list to append the indices of new track states to
  ///
  /// @return success, or an error which stops the track finding
  Result<void> createHoleTrackStates(
      State& state, CkfTypes::BranchVector<TrackIndexType>& trackStates) const {
    trackStates.push_back(derived().createTrackStateWithoutMeasurement(
        state, derived().isMeasurementExpected(state)));
    return Result<void>::success();
  }

  /// Whether a measurement is expected on this surface if the selection keeps
  /// none, e.g. false for a dead module or a prediction outside of the active
  /// area.
  ///
  /// @param state the state of the current surface
  ///
  /// @return true, i.e. a missing measurement is a hole
  bool isMeasurementExpected(const State& state) const {
    static_cast<void>(state);

    return true;
  }

  /// Pick the measurements which are to become track states on this surface.
  ///
  /// This is the main customization point. Implementations append to
  /// @p candidates in the order the track states are to be created in, which
  /// for the combinatorial Kalman filter should be the most compatible one
  /// first. If it selects nothing, the creator calls
  /// @ref createHoleTrackStates.
  ///
  /// The default selects the measurement with the smallest chi2, without a
  /// cut, as a measurement and not as an outlier. This is the default of the
  /// deprecated @ref MeasurementSelector. A measurement with a non-finite chi2
  /// is never selected, so if all of them have one, nothing is selected.
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
  /// @return always success
  template <typename range_t, typename candidates_t>
  Result<void> selectMeasurements(const State& state, range_t&& range,
                                  candidates_t& candidates) const {
    double minChi2 = std::numeric_limits<double>::infinity();
    for (const auto& measurement : range) {
      const double chi2 = derived().calculatePredictedChi2(
          state, derived().calibrate(state, measurement));
      // check finiteness first, because a comparison with NaN raises an FPE
      if (!std::isfinite(chi2) || chi2 >= minChi2) {
        continue;
      }
      minChi2 = chi2;
      candidates.clear();
      candidates.push_back({measurement, chi2, false});
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
  /// The default reads the measurement in place, dispatching a dimension which
  /// is only known at runtime onto a statically sized instantiation, so that
  /// nothing is copied and Eigen never allocates.
  ///
  /// @tparam measurement_t the calibrated measurement type
  ///
  /// @param state the state of the current surface
  /// @param measurement the calibrated measurement
  ///
  /// @return the local chi2 contribution
  template <detail::CalibratedMeasurementConcept measurement_t>
  double calculatePredictedChi2(const State& state,
                                const measurement_t& measurement) const {
    if constexpr (detail::StaticCalibratedMeasurementConcept<measurement_t>) {
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
      State& state,
      const CkfTrackStateCandidate<measurement_t>& candidate) const {
    using PM = TrackStatePropMask;

    // The filtered parameters are deliberately left unallocated. Outliers
    // never get their own, and the caller allocates them right before the
    // Kalman update writes them.
    TrackStateProxy trackState = makeTrackStateOnSurface(state, PM::Calibrated);

    trackState.chi2() = candidate.chi2;
    trackState.setUncalibratedSourceLink(
        derived().sourceLink(state, candidate.measurement));
    derived().fillCalibrated(state, candidate.measurement, trackState);

    auto typeFlags = trackState.typeFlags();
    typeFlags.setHasMeasurement();
    if (candidate.isOutlier) {
      typeFlags.setIsOutlier();
    }

    return trackState.index();
  }

  /// Create and fill a single track state without a measurement.
  ///
  /// The filtered parameters are the predicted ones, as for a hole which the
  /// CKF creates.
  ///
  /// @param state the state of the current surface
  /// @param isHole flag the track state as a hole, otherwise as a state with no
  ///        expected hit
  ///
  /// @return the index of the new track state
  TrackIndexType createTrackStateWithoutMeasurement(State& state,
                                                    bool isHole) const {
    using PM = TrackStatePropMask;

    TrackStateProxy trackState = makeTrackStateOnSurface(state, PM::None);
    trackState.shareFrom(PM::Predicted, PM::Filtered);

    auto typeFlags = trackState.typeFlags();
    if (isHole) {
      typeFlags.setIsHole();
    } else {
      typeFlags.setHasNoExpectedHit();
    }

    return trackState.index();
  }

 protected:
  /// Access to the concrete creator
  /// @return the concrete creator
  const derived_t& derived() const {
    return static_cast<const derived_t&>(*this);
  }

  /// Make a track state on the current surface with the predicted parameters,
  /// the jacobian, the path length, the surface and the flags which do not
  /// depend on a measurement.
  ///
  /// The first track state on the surface owns the predicted parameters and
  /// the jacobian. All further track states share them.
  ///
  /// @param state the state of the current surface
  /// @param mask the components to allocate beyond the shared ones
  ///
  /// @return the new track state
  TrackStateProxy makeTrackStateOnSurface(State& state,
                                          TrackStatePropMask mask) const {
    using PM = TrackStatePropMask;

    const bool isFirst = !state.firstTrackState.has_value();
    if (isFirst) {
      mask |= PM::Predicted | PM::Jacobian;
    }

    TrackStateProxy trackState =
        state.trajectory->makeTrackState(mask, state.prevTip);

    if (isFirst) {
      trackState.predicted() = state.predicted();
      trackState.predictedCovariance() = state.predictedCovariance();
      trackState.jacobian() = state.jacobian();
      state.firstTrackState = trackState;
    } else {
      // subsequent track states share these with the first one
      trackState.shareFrom(*state.firstTrackState, PM::Predicted);
      trackState.shareFrom(*state.firstTrackState, PM::Jacobian);
    }

    trackState.pathLength() = state.pathLength();
    trackState.setReferenceSurface(
        state.boundParameters().referenceSurface().getSharedPtr());

    auto typeFlags = trackState.typeFlags();
    typeFlags.setHasParameters();
    if (trackState.referenceSurface().hasMaterial()) {
      typeFlags.setHasMaterial();
    }

    return trackState;
  }

  /// Select the measurements on the surface of @p state and create their
  /// track states
  ///
  /// @param state the state of the surface
  /// @param trackStates the list to append the indices of new track states to
  ///
  /// @return success or the error of the selection
  Result<void> createMeasurementTrackStates(
      State& state, CkfTypes::BranchVector<TrackIndexType>& trackStates) const {
    const Logger& logger = *state.logger;

    auto range = derived().measurementRange(state);
    // Bail out before the selection is consulted, so that an ordinary surface
    // without measurements can never turn into a selection error.
    if (std::ranges::empty(range)) {
      ACTS_VERBOSE("No measurements on surface "
                   << state.surface->geometryId());
      return Result<void>::success();
    }

    using Measurement = std::ranges::range_value_t<decltype(range)>;
    Candidates<Measurement> candidates;

    if (Result<void> selection =
            derived().selectMeasurements(state, range, candidates);
        !selection.ok()) {
      return selection.error();
    }

    for (const CkfTrackStateCandidate<Measurement>& candidate : candidates) {
      trackStates.push_back(derived().createTrackState(state, candidate));
    }
    return Result<void>::success();
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
      const Eigen::Ref<const BoundVector>& predicted,
      const Eigen::Ref<const BoundMatrix>& predictedCovariance) {
    constexpr int kMeasurementSize = calibrated_t::RowsAtCompileTime;
    static_assert(kMeasurementSize != Eigen::Dynamic,
                  "Measurement dimension must be known at compile time");

    const FixedBoundSubspaceHelper<kMeasurementSize> subspaceHelper(
        subspaceIndices);

    const Vector<kMeasurementSize> residuals =
        calibrated - subspaceHelper.projectVector(predicted);

    return (residuals.transpose() *
            (calibratedCovariance +
             subspaceHelper.projectMatrix(predictedCovariance))
                .inverse() *
            residuals)
        .eval()(0, 0);
  }
};

/// @}

}  // namespace Acts
