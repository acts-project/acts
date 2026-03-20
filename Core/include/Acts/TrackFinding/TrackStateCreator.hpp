// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/TrackStateType.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterExtensions.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

namespace Experimental {

// TODO surface token / detector element token
// TODO state object
// TODO cuts object
template <typename derived_t>
class TrackStateCreatorBase {
 public:
  /// Type alias for result of track states creation operation
  using TrackStatesResult = Result<CkfTypes::BranchVector<TrackIndexType>>;
  /// Type alias for bound state tuple containing parameters, jacobian and path
  /// length
  using BoundState = std::tuple<BoundTrackParameters, BoundMatrix, double>;

  template <typename track_state_proxy_t,
            typename track_state_container_backend_t>
  TrackStatesResult createTrackStates(
      const GeometryContext& gctx, const CalibrationContext& cctx,
      const Surface& surface, const BoundState& boundState,
      const TrackIndexType prevTip,
      std::vector<track_state_proxy_t>& trackStateCandidates,
      track_state_container_backend_t& trajectory, const Logger& logger) const {
    static_cast<void>(trackStateCandidates);

    auto state = derived().makeState(gctx, cctx, surface, boundState, prevTip,
                                     trackStateCandidates, trajectory, logger);

    TrackStatesResult result = TrackStatesResult::success({});

    ACTS_VERBOSE("Perform measurement selection for surface "
                 << surface.geometryId());

    const auto measurementRange = derived().measurementRange(state);
    ACTS_VERBOSE("Found " << measurementRange.size()
                          << " measurements on surface "
                          << surface.geometryId());
    if (measurementRange.begin() == measurementRange.end()) {
      ACTS_VERBOSE("No measurements on surface "
                   << surface.geometryId() << ". No track states created.");
      return result;
    }

    auto selectedMeasurements =
        derived().selectMeasurements(state, measurementRange);

    ACTS_VERBOSE("Selected " << selectedMeasurements.size()
                             << " measurements for surface "
                             << surface.geometryId());

    derived().trimSelectedMeasurements(state, selectedMeasurements);

    ACTS_VERBOSE("Trimmed to " << selectedMeasurements.size()
                               << " measurements for surface "
                               << surface.geometryId());

    derived().createTrackStatesImpl(state, measurementRange,
                                    selectedMeasurements);

    return result;
  }

 protected:
  template <typename track_state_proxy_t,
            typename track_state_container_backend_t>
  struct State {
    using TrackStateProxy =
        typename track_state_container_backend_t::TrackStateProxy;

    State(const GeometryContext& gctx_, const CalibrationContext& cctx_,
          const Surface& surface_, const BoundState& boundState_,
          const TrackIndexType prevTip_,
          std::vector<track_state_proxy_t>& trackStateCandidates_,
          track_state_container_backend_t& trajectory_, const Logger& logger_)
        : gctx(&gctx_),
          cctx(&cctx_),
          surface(&surface_),
          boundState(&boundState_),
          prevTip(prevTip_),
          trackStateCandidates(&trackStateCandidates_),
          trajectory(&trajectory_),
          logger(&logger_) {}

    const GeometryContext* gctx{};
    const CalibrationContext* cctx{};
    const Surface* surface{};
    const BoundState* boundState{};
    TrackIndexType prevTip{};
    std::vector<track_state_proxy_t>* trackStateCandidates{};
    track_state_container_backend_t* trajectory{};
    const Logger* logger{};

    std::uint32_t maxNumSelectedMeasurements{};
    float maxChi2Compatible{};
    float maxChi2Outlier{};

    std::optional<TrackStateProxy> firstTrackState;
  };

  template <typename track_state_proxy_t,
            typename track_state_container_backend_t>
  auto makeState(const GeometryContext& gctx, const CalibrationContext& cctx,
                 const Surface& surface, const BoundState& boundState,
                 const TrackIndexType prevTip,
                 std::vector<track_state_proxy_t>& trackStateCandidates,
                 track_state_container_backend_t& trajectory,
                 const Logger& logger) const {
    State state(gctx, cctx, surface, boundState, prevTip, trackStateCandidates,
                trajectory, logger);

    state.maxNumSelectedMeasurements =
        derived().getMaxNumSelectedMeasurements(state);
    state.maxChi2Compatible = derived().getMaxChi2Compatible(state);
    state.maxChi2Outlier = derived().getMaxChi2Outlier(state);

    ACTS_VERBOSE("Max number of selected measurements: "
                 << state.maxNumSelectedMeasurements
                 << ", max chi2 compatible: " << state.maxChi2Compatible
                 << ", max chi2 outlier: " << state.maxChi2Outlier);

    return state;
  }

  template <typename state_t>
  std::uint32_t getMaxNumSelectedMeasurements(state_t& state) const {
    static_cast<void>(state);

    return std::numeric_limits<std::uint32_t>::max();
  }

  template <typename state_t>
  float getMaxChi2Compatible(state_t& state) const {
    static_cast<void>(state);

    return std::numeric_limits<float>::max();
  }

  template <typename state_t>
  float getMaxChi2Outlier(state_t& state) const {
    static_cast<void>(state);

    return std::numeric_limits<float>::max();
  }

  enum class MeasurementClassification : std::uint8_t {
    Incompatible,
    Outlier,
    Compatible
  };

  struct MeasurementCandidate {
    std::uint32_t index{std::numeric_limits<std::uint32_t>::max()};
    float chi2{std::numeric_limits<float>::max()};
    MeasurementClassification classification{
        MeasurementClassification::Incompatible};
  };

  struct SelectedMeasurements {
    MeasurementCandidate bestCandidate;
    CkfTypes::BranchVector<MeasurementCandidate> nonIncompatible;

    std::size_t size() const { return nonIncompatible.size(); }
    bool empty() const { return nonIncompatible.empty(); }

    auto begin() const { return nonIncompatible.begin(); }
    auto end() const { return nonIncompatible.end(); }

    void push_back(const MeasurementCandidate& candidate) {
      if (candidate.chi2 < bestCandidate.chi2) {
        bestCandidate = candidate;
      }
      if (candidate.classification != MeasurementClassification::Incompatible) {
        nonIncompatible.push_back(candidate);
      }
    }

    void emplace_back(std::uint32_t index, float chi2,
                      MeasurementClassification classification) {
      push_back({index, chi2, classification});
    }

    void sort(std::uint32_t maxNumSelectedMeasurements) {
      if (maxNumSelectedMeasurements < nonIncompatible.size()) {
        nonIncompatible.resize(maxNumSelectedMeasurements);
      }

      std::ranges::sort(
          nonIncompatible, {},
          [](const MeasurementCandidate& candidate) { return candidate.chi2; });
    }

    void resize(std::uint32_t maxNumSelectedMeasurements) {
      if (maxNumSelectedMeasurements < nonIncompatible.size()) {
        nonIncompatible.resize(maxNumSelectedMeasurements);
      }
    }
  };

  template <typename state_t, typename measurement_range_t>
  SelectedMeasurements selectMeasurements(
      state_t& state, const measurement_range_t& range) const {
    SelectedMeasurements result;

    for (const auto& [i, measurement] :
         enumerate<const measurement_range_t&, std::uint32_t>(range)) {
      const float chi2 = derived().computeChi2(state, measurement);
      const MeasurementClassification classification =
          derived().classifyMeasurement(state, measurement, chi2);
      result.emplace_back(i, chi2, classification);
    }

    return result;
  }

  template <typename state_t>
  void trimSelectedMeasurements(
      state_t& state, SelectedMeasurements& selectedMeasurements) const {
    selectedMeasurements.sort(state.maxNumSelectedMeasurements);

    selectedMeasurements.resize(state.maxNumSelectedMeasurements);
  }

  template <typename state_t, typename measurement_t>
  float computeChi2(state_t& state, const measurement_t& measurement) const {
    const auto& subspaceHelper =
        derived().measurementSubspace(state, measurement);
    return derived().computeChi2Impl(state, measurement, subspaceHelper);
  }

  template <typename state_t, typename measurement_t, std::size_t Dim>
  double computeChi2Impl(state_t& state, const measurement_t& measurement,
                         FixedBoundSubspaceHelper<Dim> subspaceHelper) const {
    const auto& [boundParams, jacobian, pathLength] = *state.boundState;

    const Vector<Dim> predictedParameters =
        subspaceHelper.projectVector(boundParams.parameters());
    const SquareMatrix<Dim> predictedCovariance =
        subspaceHelper.projectMatrix(*boundParams.covariance());

    const Vector<Dim> measuredParameters =
        derived().measuredParameters(state, measurement);
    const SquareMatrix<Dim> measurementCovariance =
        derived().measurementCovariance(state, measurement);

    const Vector<Dim> residualParameters =
        measuredParameters - predictedParameters;
    const SquareMatrix<Dim> residualCovariance =
        predictedCovariance + measurementCovariance;

    const double chi2 = (residualParameters.transpose() *
                         residualCovariance.inverse() * residualParameters)
                            .eval()(0, 0);
    return chi2;
  }

  template <typename state_t, typename measurement_t>
  double computeChi2Impl(state_t& state, const measurement_t& measurement,
                         VariableBoundSubspaceHelper subspaceHelper) const {
    return visit_measurement(subspaceHelper.size(), [&](auto N) -> double {
      constexpr std::size_t kDim = decltype(N)::value;
      const FixedBoundSubspaceHelper<kDim> fixedSubspaceHelper(subspaceHelper);
      return derived().computeChi2Impl(state, measurement, fixedSubspaceHelper);
    });
  }

  template <typename state_t, typename measurement_t>
  MeasurementClassification classifyMeasurement(
      state_t& state, const measurement_t& measurement,
      const float chi2) const {
    static_cast<void>(measurement);

    if (chi2 < state.maxChi2Compatible) {
      return MeasurementClassification::Compatible;
    }
    if (chi2 < state.maxChi2Outlier) {
      return MeasurementClassification::Outlier;
    }
    return MeasurementClassification::Incompatible;
  }

  template <typename state_t, typename measurement_range_t>
  auto selectedMeasurementSourceLink(
      state_t& state, const measurement_range_t& measurements,
      const MeasurementCandidate& selectedMeasurement) const {
    return derived().measurementSourceLink(
        state, measurements[selectedMeasurement.index]);
  }

  template <typename state_t, typename measurement_range_t>
  auto selectedMeasurementSubspace(
      state_t& state, const measurement_range_t& measurements,
      const MeasurementCandidate& selectedMeasurement) const {
    return derived().measurementSubspace(
        state, measurements[selectedMeasurement.index]);
  }

  template <typename state_t, typename measurement_range_t>
  auto selectedMeasurementParameters(
      state_t& state, const measurement_range_t& measurements,
      const MeasurementCandidate& selectedMeasurement) const {
    return derived().measuredParameters(
        state, measurements[selectedMeasurement.index]);
  }

  template <typename state_t, typename measurement_range_t>
  auto selectedMeasurementCovariance(
      state_t& state, const measurement_range_t& measurements,
      const MeasurementCandidate& selectedMeasurement) const {
    return derived().measurementCovariance(
        state, measurements[selectedMeasurement.index]);
  }

  template <typename state_t, typename measurement_range_t,
            typename selected_measurement_range_t>
  void createTrackStatesImpl(
      state_t& state, const measurement_range_t& measurements,
      const selected_measurement_range_t& selectedMeasurements) const {
    const Logger& logger = *state.logger;

    if (derived().hasHole(state, selectedMeasurements)) {
      ACTS_VERBOSE("No compatible measurements on surface "
                   << state.surface->geometryId()
                   << ". No track states created.");
      derived().createHoleState(state);
    }

    for (const auto& selectedMeasurement : selectedMeasurements) {
      derived().createMeasurementState(state, measurements,
                                       selectedMeasurement);
    }
  }

  template <typename state_t>
  bool hasHole(state_t& state,
               const SelectedMeasurements& selectedMeasurements) const {
    static_cast<void>(state);

    return selectedMeasurements.empty() &&
           (selectedMeasurements.bestCandidate.classification ==
            MeasurementClassification::Incompatible);
  }

  template <typename state_t>
  TrackStateType determineTrackStateType(state_t& state) const {
    TrackStateType result;

    result.setHasParameters();

    result.setHasMeasurement();

    if (state.surface->surfaceMaterial() != nullptr) {
      result.setHasMaterial();
    }

    return result;
  }

  template <typename state_t, typename measurement_range_t>
  TrackStateType determineTrackStateType(
      state_t& state, const measurement_range_t& measurements,
      const MeasurementCandidate& selectedMeasurement) const {
    static_cast<void>(measurements);

    TrackStateType result = derived().determineTrackStateType(state);

    if (selectedMeasurement.classification ==
        MeasurementClassification::Outlier) {
      result.setIsOutlier();
    }

    return result;
  }

  template <typename state_t>
  auto createTrackState(state_t& state,
                        const TrackStateType trackStateType) const {
    TrackStatePropMask mask = TrackStatePropMask::None;
    if (!state.firstTrackState.has_value()) {
      mask |= TrackStatePropMask::Predicted | TrackStatePropMask::Jacobian;
    }

    auto trackState = state.trajectory->makeTrackState(mask, state.prevTip);

    trackState.setReferenceSurface(state.surface->getSharedPtr());

    trackState.typeFlags() = trackStateType;

    if (!state.firstTrackState.has_value()) {
      const auto& [boundParams, jacobian, pathLength] = *state.boundState;
      trackState.predicted() = boundParams.parameters();
      trackState.predictedCovariance() = *boundParams.covariance();
      trackState.jacobian() = jacobian;
      trackState.pathLength() = pathLength;
    } else {
      trackState.shareFrom(*state.firstTrackState,
                           TrackStatePropMask::Predicted);
      trackState.shareFrom(*state.firstTrackState,
                           TrackStatePropMask::Jacobian);
    }

    if (!state.firstTrackState.has_value()) {
      state.firstTrackState = trackState;
    }
    return trackState;
  }

  template <typename state_t>
  auto createHoleState(state_t& state) const {
    TrackStateType trackStateType = derived().determineTrackStateType(state);
    trackStateType.setIsHole();

    auto trackState = derived().createTrackState(state, trackStateType);

    return trackState;
  }

  template <typename state_t, typename measurement_range_t,
            typename selected_measurement_t>
  auto createMeasurementState(
      state_t& state, const measurement_range_t& measurements,
      const selected_measurement_t& selectedMeasurement) const {
    const TrackStateType trackStateType = derived().determineTrackStateType(
        state, measurements, selectedMeasurement);

    auto trackState = derived().createTrackState(state, trackStateType);

    trackState.addComponents(TrackStatePropMask::Calibrated);

    derived().postCalibrateTrackState(state, measurements, selectedMeasurement,
                                      trackState);

    return trackState;
  }

  template <typename state_t, typename measurement_range_t,
            typename selected_measurement_t>
  void postCalibrateTrackState(
      state_t& state, const measurement_range_t& measurements,
      const selected_measurement_t& selectedMeasurement,
      typename state_t::TrackStateProxy& trackState) const {
    trackState.setUncalibratedSourceLink(
        derived().selectedMeasurementSourceLink(state, measurements,
                                                selectedMeasurement));

    const auto& subspaceHelper = derived().selectedMeasurementSubspace(
        state, measurements, selectedMeasurement);

    visit_measurement(subspaceHelper.size(), [&](auto N) -> void {
      constexpr std::size_t kDim = decltype(N)::value;
      const FixedBoundSubspaceHelper<kDim> fixedSubspaceHelper(subspaceHelper);
      const Vector<kDim> measuredParameters =
          derived().selectedMeasurementParameters(state, measurements,
                                                  selectedMeasurement);
      const SquareMatrix<kDim> measurementCovariance =
          derived().selectedMeasurementCovariance(state, measurements,
                                                  selectedMeasurement);
      trackState.allocateCalibrated(measuredParameters, measurementCovariance);
    });

    trackState.setProjectorSubspaceIndices(subspaceHelper.indices());
  }

 private:
  const derived_t& derived() const {
    return static_cast<const derived_t&>(*this);
  }
};

}  // namespace Experimental

/// @brief Create track states for selected measurements associated to a surface.
///
/// - First get a source link range covering relevant measurements associated to
///   the given surface. This task is delegated to a SourceLinkAccessor.
/// - Then create temporary track states for all measurements defined
///   by a source link range, calibrate the measurements and fill the
///   the calibrated data of these track states using a dedicated calibrator
/// - The measurement selection is delegated to a dedicated measurement
///   selector.
/// - Finally add branches to the given trajectory for the selected, temporary
///   track states. The track states of these branches still lack the filtered
///    data which is to be filled by the next stage e.g. the
///    CombinatorialKalmanFilter.
/// All track states, the temporary track states and track states for selected
/// measurements, are created in the given trajectory. The resulting container
/// may become big. Thus, it is advisable to copy selected tracks and their
/// track states to a separate container after each track finding step.
///
template <typename source_link_iterator_t, typename track_container_t>
struct TrackStateCreator {
  /// Type alias for result of track states creation operation
  using TrackStatesResult = Result<CkfTypes::BranchVector<TrackIndexType>>;
  /// Type alias for track state container backend from track container
  using TrackStateContainerBackend =
      typename track_container_t::TrackStateContainerBackend;
  /// Type alias for track proxy from track container
  using TrackProxy = typename track_container_t::TrackProxy;
  /// Type alias for track state proxy from track container
  using TrackStateProxy = typename track_container_t::TrackStateProxy;
  /// Type alias for bound state tuple containing parameters, jacobian and path
  /// length
  using BoundState = std::tuple<BoundTrackParameters, BoundMatrix, double>;
  /// Type alias for container of candidate track state proxies
  using candidate_container_t =
      typename std::vector<typename track_container_t::TrackStateProxy>;

  // delegate definition to get source link ranges for a surface
  /// Type alias for delegate to access source link ranges for a surface
  using SourceLinkAccessor =
      Delegate<std::pair<source_link_iterator_t, source_link_iterator_t>(
          const Surface&)>;

  // delegate to get calibrted measurements from a source link iterator
  /// Type alias for calibrator delegate to process measurements from source
  /// links
  using Calibrator =
      typename KalmanFitterExtensions<TrackStateContainerBackend>::Calibrator;

  // delegate to select measurements from a track state range
  /// Type alias for delegate to select measurements from track state candidates
  using MeasurementSelector =
      Delegate<Result<std::pair<typename candidate_container_t::iterator,
                                typename candidate_container_t::iterator>>(
          candidate_container_t& trackStates, bool&, const Logger&)>;

  /// The source link accessor will return an source link range for a surface
  /// which link to the associated measurements.
  SourceLinkAccessor sourceLinkAccessor;

  /// The Calibrator is a dedicated calibration algorithm that allows to
  /// calibrate measurements using track information, this could be e.g. sagging
  /// for wires, module deformations, etc.
  Calibrator calibrator{DelegateFuncTag<
      detail::voidFitterCalibrator<TrackStateContainerBackend>>{}};

  /// Delegate for measurement selection on surfaces
  MeasurementSelector measurementSelector{
      DelegateFuncTag<voidMeasurementSelector>{}};

  /// @brief extend the trajectory onto the given surface.
  ///
  /// @param gctx The geometry context to be used for this task
  /// @param calibrationContext The calibration context used to fill the calibrated data
  /// @param surface The surface onto which the trajectory is extended
  /// @param boundState the predicted bound state on the given surface
  /// @param prevTip the tip of the trajectory which is to be extended
  /// @param trackStateCandidates a temporary buffer which can be used to
  ///        to keep track of newly created temporary track states.
  /// @param trajectory the trajectory to be extended.
  /// @param logger a logger for messages.
  ///
  /// @return a list of indices of newly created track states which extend the
  ///    trajectory onto the given surface and match the bound state, or an
  ///    error.
  ///
  /// Extend or branch the trajectory onto the given surface. This may create
  /// new track states using measurements which match the predicted bound state.
  /// This may create multiple branches. The new track states still miss the
  /// "filtered" data.
  Result<CkfTypes::BranchVector<TrackIndexType>> createTrackStates(
      const GeometryContext& gctx, const CalibrationContext& calibrationContext,
      [[maybe_unused]] const Surface& surface, const BoundState& boundState,
      TrackIndexType prevTip,
      std::vector<TrackStateProxy>& trackStateCandidates,
      TrackStateContainerBackend& trajectory, const Logger& logger) const {
    TrackStatesResult tsRes = TrackStatesResult::success({});
    using SourceLinkRange = decltype(sourceLinkAccessor(surface));
    SourceLinkRange slRange = sourceLinkAccessor(surface);
    if (slRange.first != slRange.second) {
      auto [slBegin, slEnd] = slRange;
      tsRes = createSourceLinkTrackStates(
          gctx, calibrationContext, surface, boundState, slBegin, slEnd,
          prevTip, trackStateCandidates, trajectory, logger);
    }
    return tsRes;
  }

  /// Create track states for selected measurements given by the source links
  ///
  /// @param gctx The current geometry context
  /// @param calibrationContext pointer to the current calibration context
  /// @param surface the surface the sourceLinks are associated to
  /// @param boundState Bound state from the propagation on this surface
  /// @param slBegin Begin iterator for sourceLinks
  /// @param slEnd End iterator for sourceLinks
  /// @param prevTip Index pointing at previous trajectory state (i.e. tip)
  /// @param trackStateCandidates a temporary buffer which can be used to
  ///        to keep track of newly created temporary track states.
  /// @param trajectory the trajectory to which new track states for selected measurements will be added
  /// @param logger the logger for messages.
  /// @return Result containing vector of track state indices or error
  Result<CkfTypes::BranchVector<TrackIndexType>> createSourceLinkTrackStates(
      const GeometryContext& gctx, const CalibrationContext& calibrationContext,
      [[maybe_unused]] const Surface& surface, const BoundState& boundState,
      const source_link_iterator_t& slBegin,
      const source_link_iterator_t& slEnd, TrackIndexType prevTip,
      std::vector<TrackStateProxy>& trackStateCandidates,
      TrackStateContainerBackend& trajectory, const Logger& logger) const {
    using PM = TrackStatePropMask;

    using ResultTrackStateList = Result<CkfTypes::BranchVector<TrackIndexType>>;
    ResultTrackStateList resultTrackStateList{
        CkfTypes::BranchVector<TrackIndexType>()};
    const auto& [boundParams, jacobian, pathLength] = boundState;

    trackStateCandidates.clear();
    if constexpr (std::ranges::random_access_range<source_link_iterator_t>) {
      trackStateCandidates.reserve(std::distance(slBegin, slEnd));
    }

    // Calibrate all the source links on the surface since the selection has
    // to be done based on calibrated measurement
    for (auto it = slBegin; it != slEnd; ++it) {
      // get the source link
      const auto sourceLink = *it;

      // prepare the track state
      PM mask = PM::Predicted | PM::Jacobian | PM::Calibrated;
      if (it != slBegin) {
        // not the first TrackState, only need uncalibrated and calibrated
        mask = PM::Calibrated;
      }

      ACTS_VERBOSE("Create temp track state with mask: " << mask);
      // Temporary and final track states are created in the same
      // trajectory, which could lead to very large containers.

      // CAREFUL! This trackstate has a previous index that is not in this
      // MultiTrajectory Visiting backwards from this track state will
      // fail!
      auto ts = trajectory.makeTrackState(mask, prevTip);

      if (it == slBegin) {
        // only set these for first
        ts.predicted() = boundParams.parameters();
        if (boundParams.covariance()) {
          ts.predictedCovariance() = *boundParams.covariance();
        }
        ts.jacobian() = jacobian;
      } else {
        // subsequent track states can reuse
        auto& first = trackStateCandidates.front();
        ts.shareFrom(first, PM::Predicted);
        ts.shareFrom(first, PM::Jacobian);
      }

      ts.pathLength() = pathLength;
      ts.setReferenceSurface(boundParams.referenceSurface().getSharedPtr());

      // now calibrate the track state
      calibrator(gctx, calibrationContext, sourceLink, ts);

      trackStateCandidates.push_back(ts);
    }

    bool isOutlier = false;
    Result<std::pair<typename std::vector<TrackStateProxy>::iterator,
                     typename std::vector<TrackStateProxy>::iterator>>
        selectorResult =
            measurementSelector(trackStateCandidates, isOutlier, logger);
    if (!selectorResult.ok()) {
      ACTS_DEBUG("Selection of calibrated measurements failed: "
                 << selectorResult.error().message());
      resultTrackStateList =
          ResultTrackStateList::failure(selectorResult.error());
    } else {
      auto selectedTrackStateRange = *selectorResult;
      resultTrackStateList = processSelectedTrackStates(
          selectedTrackStateRange.first, selectedTrackStateRange.second,
          trajectory, isOutlier, logger);
    }

    return resultTrackStateList;
  }

  /// Create track states for the given trajectory from candidate track states
  ///
  /// @param begin begin iterator of the list of candidate track states
  /// @param end end iterator of the list of candidate track states
  /// @param trackStates the trajectory to which the new track states are added
  /// @param isOutlier true if the candidate(s) is(are) an outlier(s).
  /// @param logger the logger for messages
  /// @return Result containing vector of track state indices or error
  Result<CkfTypes::BranchVector<TrackIndexType>> processSelectedTrackStates(
      typename std::vector<TrackStateProxy>::const_iterator begin,
      typename std::vector<TrackStateProxy>::const_iterator end,
      TrackStateContainerBackend& trackStates, bool isOutlier,
      const Logger& logger) const {
    using PM = TrackStatePropMask;

    using ResultTrackStateList = Result<CkfTypes::BranchVector<TrackIndexType>>;
    ResultTrackStateList resultTrackStateList{
        CkfTypes::BranchVector<TrackIndexType>()};
    CkfTypes::BranchVector<TrackIndexType>& trackStateList =
        *resultTrackStateList;
    trackStateList.reserve(end - begin);

    std::optional<TrackStateProxy> firstTrackState{std::nullopt};
    for (auto it = begin; it != end; ++it) {
      auto& candidateTrackState = *it;

      PM mask = PM::Predicted | PM::Filtered | PM::Jacobian | PM::Calibrated;
      if (it != begin) {
        // subsequent track states don't need storage for these as they will
        // be shared
        mask &= ~PM::Predicted & ~PM::Jacobian;
      }
      if (isOutlier) {
        // outlier won't have separate filtered parameters
        mask &= ~PM::Filtered;
      }

      // copy this trackstate into fitted states MultiTrajectory
      auto trackState =
          trackStates.makeTrackState(mask, candidateTrackState.previous());
      ACTS_VERBOSE("Create SourceLink output track state #"
                   << trackState.index() << " with mask: " << mask);

      if (it != begin) {
        // assign indices pointing to first track state
        trackState.shareFrom(*firstTrackState, PM::Predicted);
        trackState.shareFrom(*firstTrackState, PM::Jacobian);
      } else {
        firstTrackState = trackState;
      }

      // either copy ALL or everything except for predicted and jacobian
      trackState.copyFrom(candidateTrackState, mask, false);

      auto typeFlags = trackState.typeFlags();
      typeFlags.setHasParameters();
      typeFlags.setHasMeasurement();
      if (trackState.referenceSurface().surfaceMaterial() != nullptr) {
        typeFlags.setHasMaterial();
      }
      if (isOutlier) {
        // propagate information that this is an outlier state
        ACTS_VERBOSE(
            "Creating outlier track state with tip = " << trackState.index());
        typeFlags.setIsOutlier();
      }

      trackStateList.push_back(trackState.index());
    }
    return resultTrackStateList;
  }

  /// Default measurement selector which will return all measurements
  /// @param candidates Measurement track state candidates
  /// @return Iterator pair representing the range of all candidates
  static Result<std::pair<typename std::vector<TrackStateProxy>::iterator,
                          typename std::vector<TrackStateProxy>::iterator>>
  voidMeasurementSelector(typename std::vector<TrackStateProxy>& candidates,
                          bool& /*isOutlier*/, const Logger& /*logger*/) {
    return std::pair{candidates.begin(), candidates.end()};
  };
};

}  // namespace Acts
