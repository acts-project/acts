// This file is part of the Acts project.
//
// Copyright (C) 2016-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterError.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/VoidFitterComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <functional>
#include <limits>
#include <memory>
#include <type_traits>
#include <unordered_map>

namespace Acts {

/// Track quality summary for one trajectory.
///
/// This could be used to decide if a track is to be recorded when the
/// filtering is done or to be terminated due to its bad quality
/// @todo: add other useful info, e.g. chi2
struct CombinatorialKalmanFilterTipState {
  // Number of passed sensitive surfaces
  std::size_t nSensitiveSurfaces = 0;
  // Number of track states
  std::size_t nStates = 0;
  // Number of (non-outlier) measurements
  std::size_t nMeasurements = 0;
  // Number of outliers
  std::size_t nOutliers = 0;
  // Number of holes
  std::size_t nHoles = 0;
};

enum class CombinatorialKalmanFilterBranchStopperResult {
  Continue,
  StopAndDrop,
  StopAndKeep,
};

/// Extension struct which holds the delegates to customize the CKF behavior
template <typename traj_t>
struct CombinatorialKalmanFilterExtensions {
  using candidate_container_t =
      typename std::vector<typename traj_t::TrackStateProxy>;
  using BranchStopperResult = CombinatorialKalmanFilterBranchStopperResult;

  using Calibrator = typename KalmanFitterExtensions<traj_t>::Calibrator;
  using Updater = typename KalmanFitterExtensions<traj_t>::Updater;
  using MeasurementSelector =
      Delegate<Result<std::pair<typename candidate_container_t::iterator,
                                typename candidate_container_t::iterator>>(
          candidate_container_t& trackStates, bool&, const Logger&)>;
  using BranchStopper =
      Delegate<BranchStopperResult(const CombinatorialKalmanFilterTipState&,
                                   typename traj_t::TrackStateProxy&)>;

  /// The Calibrator is a dedicated calibration algorithm that allows to
  /// calibrate measurements using track information, this could be e.g. sagging
  /// for wires, module deformations, etc.
  Calibrator calibrator{
      DelegateFuncTag<detail::voidFitterCalibrator<traj_t>>{}};

  /// The updater incorporates measurement information into the track parameters
  Updater updater{DelegateFuncTag<detail::voidFitterUpdater<traj_t>>{}};

  /// The measurement selector is called during the filtering by the Actor.
  MeasurementSelector measurementSelector{
      DelegateFuncTag<voidMeasurementSelector>{}};

  /// The branch stopper is called during the filtering by the Actor.
  BranchStopper branchStopper{DelegateFuncTag<voidBranchStopper>{}};

 private:
  /// Default measurement selector which will return all measurements
  /// @param candidates Measurement track state candidates
  static Result<std::pair<
      typename std::vector<typename traj_t::TrackStateProxy>::iterator,
      typename std::vector<typename traj_t::TrackStateProxy>::iterator>>
  voidMeasurementSelector(
      typename std::vector<typename traj_t::TrackStateProxy>& candidates,
      bool& /*isOutlier*/, const Logger& /*logger*/) {
    return std::pair{candidates.begin(), candidates.end()};
  };

  /// Default branch stopper which will never stop
  /// @return false
  static BranchStopperResult voidBranchStopper(
      const CombinatorialKalmanFilterTipState& /*tipState*/,
      typename traj_t::TrackStateProxy& /*trackState*/) {
    return BranchStopperResult::Continue;
  }
};

/// Delegate type that retrieves a range of source links to for a given surface
/// to be processed by the CKF
template <typename source_link_iterator_t>
using SourceLinkAccessorDelegate =
    Delegate<std::pair<source_link_iterator_t, source_link_iterator_t>(
        const Surface&)>;

/// expected max number of track states that are expected to be added by
/// stateCandidateCreator
/// @note if the number of states exceeds this number dynamic memory allocation will occur.
//.       the number is chosen to yield a container size of 64 bytes.
static constexpr std::size_t s_maxBranchesPerSurface = 10;

/// Combined options for the combinatorial Kalman filter.
///
/// @tparam source_link_iterator_t Type of the source link iterator
/// @tparam traj_t Type of the trajectory
template <typename source_link_iterator_t, typename traj_t>
struct CombinatorialKalmanFilterOptions {
  using SourceLinkIterator = source_link_iterator_t;
  using SourceLinkAccessor = SourceLinkAccessorDelegate<source_link_iterator_t>;

  /// PropagatorOptions with context
  ///
  /// @param gctx The geometry context for this track finding/fitting
  /// @param mctx The magnetic context for this track finding/fitting
  /// @param cctx The calibration context for this track finding/fitting
  /// @param accessor_ The source link accessor
  /// @param extensions_ The extension struct
  /// @param pOptions The plain propagator options
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  CombinatorialKalmanFilterOptions(
      const GeometryContext& gctx, const MagneticFieldContext& mctx,
      std::reference_wrapper<const CalibrationContext> cctx,
      SourceLinkAccessor accessor_,
      CombinatorialKalmanFilterExtensions<traj_t> extensions_,
      const PropagatorPlainOptions& pOptions, bool mScattering = true,
      bool eLoss = true)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        sourcelinkAccessor(std::move(accessor_)),
        extensions(extensions_),
        propagatorPlainOptions(pOptions),
        multipleScattering(mScattering),
        energyLoss(eLoss) {}

  /// Contexts are required and the options must not be default-constructible.
  CombinatorialKalmanFilterOptions() = delete;

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  /// The source link accessor
  SourceLinkAccessor sourcelinkAccessor;

  /// The filter extensions
  CombinatorialKalmanFilterExtensions<traj_t> extensions;

  /// The trivial propagator options
  PropagatorPlainOptions propagatorPlainOptions;

  /// The target surface
  /// @note This is useful if the filtering should be terminated at a
  ///       certain surface
  const Surface* targetSurface = nullptr;

  using BoundState = std::tuple<BoundTrackParameters, BoundMatrix, double>;

  /// Delegate definition to create track states for selected measurements
  /// @note expected to iterator over the given sourcelink range,
  ///       select measurements, and create track states for
  ///       which new tips are to be created, more over the outlier
  ///       flag should be set for states that are outlier.
  ///
  /// @param geoContext The current geometry context
  /// @param calibrationContext pointer to the current calibration context
  /// @param surface the surface at which new track states are to be created
  /// @param boundState the current bound state of the trajectory
  /// @param slBegin Begin iterator for sourcelinks
  /// @param slEnd End iterator for sourcelinks
  /// @param prevTip Index pointing at previous trajectory state (i.e. tip)
  /// @param bufferTrajectory a temporary trajectory which can be used to create temporary track states
  /// @param trackStateCandidates a temporary buffer that can be used to collect track states
  /// @param trajectory the trajectory to which the new states are to be added
  /// @param logger a logger for messages
  using TrackStateCandidateCreator =
      Delegate<Result<boost::container::small_vector<
          typename traj_t::TrackStateProxy::IndexType,
          s_maxBranchesPerSurface>>(
          const GeometryContext& geoContext,
          const CalibrationContext& calibrationContext, const Surface& surface,
          const BoundState& boundState, source_link_iterator_t slBegin,
          source_link_iterator_t slEnd, std::size_t prevTip,
          traj_t& bufferTrajectory,
          std::vector<typename traj_t::TrackStateProxy>& trackStateCandidates,
          traj_t& trajectory, const Logger& logger)>;

  /// The delegate to create new track states.
  TrackStateCandidateCreator trackStateCandidateCreator;

  /// Whether to consider multiple scattering.
  bool multipleScattering = true;

  /// Whether to consider energy loss.
  bool energyLoss = true;
};

template <typename traj_t>
struct CombinatorialKalmanFilterResult {
  /// Fitted states that the actor has handled.
  traj_t* fittedStates{nullptr};

  /// This is used internally to store candidate trackstates
  std::shared_ptr<traj_t> stateBuffer;
  std::vector<typename traj_t::TrackStateProxy> trackStateCandidates;

  /// This is the indices of the 'tip' of the tracks stored in multitrajectory.
  /// This corresponds to the last measurement state in the multitrajectory.
  std::vector<MultiTrajectoryTraits::IndexType> lastMeasurementIndices;

  /// This is the indices of the 'tip' of the tracks stored in multitrajectory.
  /// This corresponds to the last state in the multitrajectory.
  std::vector<MultiTrajectoryTraits::IndexType> lastTrackIndices;

  /// The Parameters at the provided surface for separate tracks
  std::unordered_map<MultiTrajectoryTraits::IndexType, BoundTrackParameters>
      fittedParameters;

  /// The indices of the 'tip' of the unfinished tracks
  std::vector<std::pair<MultiTrajectoryTraits::IndexType,
                        CombinatorialKalmanFilterTipState>>
      activeTips;

  /// The indices of track states and corresponding source links on different
  /// surfaces
  std::unordered_map<const Surface*,
                     std::unordered_map<std::size_t, std::size_t>>
      sourcelinkTips;

  /// Indicator if track finding has been done
  bool finished = false;

  /// Last encountered error
  Result<void> lastError{Result<void>::success()};

  /// Path limit aborter
  PathLimitReached pathLimitReached;
};

/// Combinatorial Kalman filter to find tracks.
///
/// @tparam propagator_t Type of the propagator
///
/// The CombinatorialKalmanFilter contains an Actor and a Sequencer sub-class.
/// The Sequencer has to be part of the Navigator of the Propagator in order to
/// initialize and provide the measurement surfaces.
///
/// The Actor is part of the Propagation call and does the Kalman update and
/// eventually the smoothing. Updater and Calibrator are given to the Actor for
/// further use:
/// - The Updater is the implemented kalman updater formalism, it
///   runs via a visitor pattern through the measurements.
///
/// Measurements are not required to be ordered for the
/// CombinatorialKalmanFilter, measurement ordering needs to be figured out by
/// the navigation of the propagator.
///
/// The void components are provided mainly for unit testing.
///
template <typename propagator_t, typename traj_t>
class CombinatorialKalmanFilter {
 public:
  /// Default constructor is deleted
  CombinatorialKalmanFilter() = delete;
  /// Constructor from arguments
  CombinatorialKalmanFilter(propagator_t pPropagator,
                            std::unique_ptr<const Logger> _logger =
                                getDefaultLogger("CKF", Logging::INFO))
      : m_propagator(std::move(pPropagator)),
        m_logger(std::move(_logger)),
        m_actorLogger{m_logger->cloneWithSuffix("Actor")},
        m_updaterLogger{m_logger->cloneWithSuffix("Updater")} {}

 private:
  using KalmanNavigator = typename propagator_t::Navigator;

  using BoundState = std::tuple<BoundTrackParameters, BoundMatrix, double>;

  /// The propagator for the transport and material update
  propagator_t m_propagator;

  std::unique_ptr<const Logger> m_logger;
  std::shared_ptr<const Logger> m_actorLogger;
  std::shared_ptr<const Logger> m_updaterLogger;

  const Logger& logger() const { return *m_logger; }

  struct DefaultTrackStateCreator {
    typename CombinatorialKalmanFilterExtensions<traj_t>::Calibrator calibrator;
    typename CombinatorialKalmanFilterExtensions<traj_t>::MeasurementSelector
        measurementSelector;

    /// Create track states for selected measurements given by the source links
    ///
    /// @param gctx The current geometry context
    /// @param calibrationContext pointer to the current calibration context
    /// @param surface the surface the sourceLinks are associated to
    /// @param boundState Bound state from the propagation on this surface
    /// @param slBegin Begin iterator for sourcelinks
    /// @param slEnd End iterator for sourcelinks
    /// @param prevTip Index pointing at previous trajectory state (i.e. tip)
    /// @param bufferTrajectory a buffer for temporary candidate track states
    /// @param trackStateCandidates a buffer for temporary track state proxies for candidates
    /// @param trajectory the trajectory to which new track states for selected measurements will be added
    /// @param logger the logger for messages.
    template <typename source_link_iterator_t>
    Result<boost::container::small_vector<
        typename traj_t::TrackStateProxy::IndexType, s_maxBranchesPerSurface>>
    createSourceLinkTrackStates(
        const GeometryContext& gctx,
        const CalibrationContext& calibrationContext,
        [[maybe_unused]] const Surface& surface, const BoundState& boundState,
        source_link_iterator_t slBegin, source_link_iterator_t slEnd,
        std::size_t prevTip, traj_t& bufferTrajectory,
        std::vector<typename traj_t::TrackStateProxy>& trackStateCandidates,
        traj_t& trajectory, const Logger& logger) const {
      using ResultTrackStateList = Acts::Result<boost::container::small_vector<
          typename traj_t::TrackStateProxy::IndexType,
          s_maxBranchesPerSurface>>;
      ResultTrackStateList resultTrackStateList{boost::container::small_vector<
          typename traj_t::TrackStateProxy::IndexType,
          s_maxBranchesPerSurface>()};
      const auto& [boundParams, jacobian, pathLength] = boundState;

      trackStateCandidates.clear();
      if constexpr (std::is_same_v<
                        typename std::iterator_traits<
                            source_link_iterator_t>::iterator_category,
                        std::random_access_iterator_tag>) {
        trackStateCandidates.reserve(std::distance(slBegin, slEnd));
      }

      bufferTrajectory.clear();

      using PM = TrackStatePropMask;

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
        // CAREFUL! This trackstate has a previous index that is not in this
        // MultiTrajectory Visiting brackwards from this track state will
        // fail!
        auto ts = bufferTrajectory.makeTrackState(mask, prevTip);

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
      Result<std::pair<
          typename std::vector<typename traj_t::TrackStateProxy>::iterator,
          typename std::vector<typename traj_t::TrackStateProxy>::iterator>>
          selectorResult =
              measurementSelector(trackStateCandidates, isOutlier, logger);
      if (!selectorResult.ok()) {
        ACTS_ERROR("Selection of calibrated measurements failed: "
                   << selectorResult.error());
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
    /// @param fittedStates the trajectory to which the new track states are added
    /// @param isOutlier true if the candidate(s) is(are) an outlier(s).
    /// @param logger the logger for messages
    Result<boost::container::small_vector<
        typename traj_t::TrackStateProxy::IndexType, s_maxBranchesPerSurface>>
    processSelectedTrackStates(
        typename std::vector<typename traj_t::TrackStateProxy>::const_iterator
            begin,
        typename std::vector<typename traj_t::TrackStateProxy>::const_iterator
            end,
        traj_t& fittedStates, bool isOutlier, const Logger& logger) const {
      Acts::Result<boost::container::small_vector<
          typename traj_t::TrackStateProxy::IndexType, s_maxBranchesPerSurface>>
          resultTrackStateList{boost::container::small_vector<
              typename traj_t::TrackStateProxy::IndexType,
              s_maxBranchesPerSurface>()};
      boost::container::small_vector<
          typename traj_t::TrackStateProxy::IndexType, s_maxBranchesPerSurface>&
          trackStateList = *resultTrackStateList;
      trackStateList.reserve(end - begin);
      using PM = TrackStatePropMask;

      std::optional<typename traj_t::TrackStateProxy> firstTrackState{
          std::nullopt};
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
        typename traj_t::TrackStateProxy trackState =
            fittedStates.makeTrackState(mask, candidateTrackState.previous());
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
        trackState.allocateCalibrated(candidateTrackState.calibratedSize());
        trackState.copyFrom(candidateTrackState, mask, false);

        auto typeFlags = trackState.typeFlags();
        if (trackState.referenceSurface().surfaceMaterial() != nullptr) {
          typeFlags.set(TrackStateFlag::MaterialFlag);
        }
        typeFlags.set(TrackStateFlag::ParameterFlag);

        if (isOutlier) {
          // propagate information that this is an outlier state
          ACTS_VERBOSE(
              "Creating outlier track state with tip = " << trackState.index());
          // Set the outlier flag needed by the next step to identify outlier
          // states
          typeFlags.set(TrackStateFlag::OutlierFlag);
        }
        trackStateList.push_back(trackState.index());
      }
      return resultTrackStateList;
    }
  };

  /// @brief Propagator Actor plugin for the CombinatorialKalmanFilter
  ///
  /// @tparam source_link_accessor_t The type of source link accessor
  /// @tparam parameters_t The type of parameters used for "local" parameters.
  ///
  /// The CombinatorialKalmanFilter Actor does not rely on the measurements to
  /// be sorted along the track.
  template <typename source_link_accessor_t, typename parameters_t>
  class Actor {
   public:
    using TipState = CombinatorialKalmanFilterTipState;
    using BoundState = std::tuple<parameters_t, BoundMatrix, double>;
    using CurvilinearState =
        std::tuple<CurvilinearTrackParameters, BoundMatrix, double>;
    /// Broadcast the result_type
    using result_type = CombinatorialKalmanFilterResult<traj_t>;

    /// The target surface aborter
    SurfaceReached targetReached{std::numeric_limits<double>::lowest()};

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// Whether to consider energy loss.
    bool energyLoss = true;

    /// Calibration context for the finding run
    const CalibrationContext* calibrationContextPtr{nullptr};

    /// @brief CombinatorialKalmanFilter actor operation
    ///
    /// @tparam propagator_state_t Type of the Propagator state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper is the stepper in use
    /// @param navigator is the navigator in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void operator()(propagator_state_t& state, const stepper_t& stepper,
                    const navigator_t& navigator, result_type& result,
                    const Logger& /*logger*/) const {
      assert(result.fittedStates && "No MultiTrajectory set");

      if (result.finished) {
        return;
      }

      ACTS_VERBOSE("CombinatorialKalmanFilter step");

      // Initialize path limit reached aborter
      if (result.pathLimitReached.internalLimit ==
          std::numeric_limits<double>::max()) {
        detail::setupLoopProtection(state, stepper, result.pathLimitReached,
                                    true, logger());
      }

      // Update:
      // - Waiting for a current surface
      if (auto surface = navigator.currentSurface(state.navigation);
          surface != nullptr) {
        // There are three scenarios:
        // 1) The surface is in the measurement map
        // -> Select source links
        // -> Perform the kalman update for selected non-outlier source links
        // -> Add track states in multitrajectory. Multiple states mean branch
        // splitting.
        // -> Call branch stopper to justify each branch
        // -> If there is non-outlier state, update stepper information
        // 2) The surface is not in the measurement map but with material or is
        // an active surface
        // -> Add a hole or passive material state in multitrajectory
        // -> Call branch stopper to justify the branch
        // 3) The surface is neither in the measurement map nor with material
        // -> Do nothing
        ACTS_VERBOSE("Perform filter step");
        auto res = filter(surface, state, stepper, navigator, result);
        if (!res.ok()) {
          ACTS_ERROR("Error in filter: " << res.error());
          result.lastError = res.error();
        }
      }

      const bool isEndOfWorldReached =
          endOfWorldReached(state, stepper, navigator, logger());
      const bool isPathLimitReached =
          result.pathLimitReached(state, stepper, navigator, logger());
      const bool isTargetReached =
          targetReached(state, stepper, navigator, logger());
      if (isEndOfWorldReached || isPathLimitReached || isTargetReached) {
        if (isEndOfWorldReached) {
          ACTS_VERBOSE("End of world reached");
        } else if (isPathLimitReached) {
          ACTS_VERBOSE("Path limit reached");
        } else if (isTargetReached) {
          ACTS_VERBOSE("Target surface reached");

          // Bind the parameter to the target surface
          auto res = stepper.boundState(state.stepping, *targetReached.surface);
          if (!res.ok()) {
            ACTS_ERROR("Error while acquiring bound state for target surface: "
                       << res.error() << " " << res.error().message());
            result.lastError = res.error();
          } else if (!result.activeTips.empty()) {
            const auto& fittedState = *res;
            std::size_t currentTip = result.activeTips.back().first;
            // Assign the fitted parameters
            result.fittedParameters.emplace(
                currentTip, std::get<BoundTrackParameters>(fittedState));
          }

          stepper.releaseStepSize(state.stepping, ConstrainedStep::actor);
        }

        if (!result.activeTips.empty()) {
          // Record the active tip as trajectory entry indices and remove it
          // from the list
          storeLastActiveTip(result);
          // Remove the tip from list of active tips
          result.activeTips.erase(result.activeTips.end() - 1);
        }
        // If no more active tip, done with filtering; Otherwise, reset
        // propagation state to track state at last tip of active tips
        if (result.activeTips.empty()) {
          ACTS_VERBOSE("Kalman filtering finds "
                       << result.lastTrackIndices.size() << " tracks");
          result.finished = true;
        } else {
          ACTS_VERBOSE("Propagation jumps to branch with tip = "
                       << result.activeTips.back().first);
          reset(state, stepper, navigator, result);
        }
      }
    }

    /// @brief CombinatorialKalmanFilter actor operation: reset propagation
    ///
    /// @tparam propagator_state_t Type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper is the stepper in use
    /// @param navigator is the navigator in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void reset(propagator_state_t& state, const stepper_t& stepper,
               const navigator_t& navigator, result_type& result) const {
      auto currentState =
          result.fittedStates->getTrackState(result.activeTips.back().first);

      // Reset the stepping state
      stepper.resetState(state.stepping, currentState.filtered(),
                         currentState.filteredCovariance(),
                         currentState.referenceSurface(),
                         state.options.maxStepSize);

      // Reset the navigation state
      // Set targetSurface to nullptr for forward filtering; it's only needed
      // after smoothing
      state.navigation =
          navigator.makeState(&currentState.referenceSurface(), nullptr);
      navigator.initialize(state, stepper);

      // No Kalman filtering for the starting surface, but still need
      // to consider the material effects here
      materialInteractor(navigator.currentSurface(state.navigation), state,
                         stepper, navigator, MaterialUpdateStage::PostUpdate);

      detail::setupLoopProtection(state, stepper, result.pathLimitReached, true,
                                  logger());
    }

    /// @brief CombinatorialKalmanFilter actor operation:
    /// - filtering for all measurement(s) on surface
    /// - store selected track states in multiTrajectory
    /// - update propagator state to the (last) selected track state
    ///
    /// @tparam propagator_state_t Type of the Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param result The mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    Result<void> filter(const Surface* surface, propagator_state_t& state,
                        const stepper_t& stepper, const navigator_t& navigator,
                        result_type& result) const {
      std::size_t nBranchesOnSurface = 0;
      // Count the number of source links on the surface
      auto [slBegin, slEnd] = m_sourcelinkAccessor(*surface);
      if (slBegin != slEnd) {
        // Screen output message
        ACTS_VERBOSE("Measurement surface " << surface->geometryId()
                                            << " detected.");

        // Transport the covariance to the surface
        stepper.transportCovarianceToBound(state.stepping, *surface);

        // Update state and stepper with pre material effects
        materialInteractor(surface, state, stepper, navigator,
                           MaterialUpdateStage::PreUpdate);

        // Bind the transported state to the current surface
        auto boundStateRes =
            stepper.boundState(state.stepping, *surface, false);
        if (!boundStateRes.ok()) {
          return boundStateRes.error();
        }
        auto& boundState = *boundStateRes;
        auto& [boundParams, jacobian, pathLength] = boundState;
        boundParams.covariance() = state.stepping.cov;

        // Retrieve the previous tip and its state
        // The states created on this surface will have the common previous tip
        std::size_t prevTip = kTrackIndexInvalid;
        TipState prevTipState;
        if (!result.activeTips.empty()) {
          prevTip = result.activeTips.back().first;
          prevTipState = result.activeTips.back().second;
          // New state is to be added. Remove the last tip from active tips
          result.activeTips.erase(result.activeTips.end() - 1);
        }

        // Create trackstates for all source links (will be filtered later)
        // Results are stored in result => no return value

        using TrackStatesResult = Acts::Result<boost::container::small_vector<
            typename traj_t::TrackStateProxy::IndexType,
            s_maxBranchesPerSurface>>;

        TrackStatesResult tsRes = trackStateCandidateCreator(
            state.geoContext, *calibrationContextPtr, *surface, boundState,
            slBegin, slEnd, prevTip, *result.stateBuffer,
            result.trackStateCandidates, *result.fittedStates, logger());

        if (!tsRes.ok()) {
          ACTS_ERROR(
              "Processing of selected track states failed: " << tsRes.error())
          return tsRes.error();
        }
        Result<std::tuple<unsigned int, bool>> procRes = processNewTrackStates(
            state.geoContext, prevTipState, *tsRes, result);
        if (!procRes.ok()) {
          ACTS_ERROR(
              "Processing of selected track states failed: " << procRes.error())
          return procRes.error();
        }
        auto [nNewBranchesOnSurface, isOutlier] = *procRes;
        nBranchesOnSurface = nNewBranchesOnSurface;

        if (nBranchesOnSurface > 0 && !isOutlier) {
          // If there are measurement track states on this surface
          ACTS_VERBOSE("Filtering step successful with " << nBranchesOnSurface
                                                         << " branches");
          // Update stepping state using filtered parameters of last track
          // state on this surface
          auto ts = result.fittedStates->getTrackState(
              result.activeTips.back().first);
          stepper.update(state.stepping,
                         MultiTrajectoryHelpers::freeFiltered(
                             state.options.geoContext, ts),
                         ts.filtered(), ts.filteredCovariance(), *surface);
          ACTS_VERBOSE("Stepping state is updated with filtered parameter:");
          ACTS_VERBOSE("-> " << ts.filtered().transpose()
                             << " of track state with tip = "
                             << result.activeTips.back().first);
        }

        // Update state and stepper with post material effects
        materialInteractor(surface, state, stepper, navigator,
                           MaterialUpdateStage::PostUpdate);
      } else if (surface->associatedDetectorElement() != nullptr ||
                 surface->surfaceMaterial() != nullptr) {
        // No splitting on the surface without source links. Set it to one
        // first, but could be changed later
        nBranchesOnSurface = 1;

        // Retrieve the previous tip and its state
        std::size_t prevTip = kTrackIndexInvalid;
        TipState tipState;
        if (!result.activeTips.empty()) {
          prevTip = result.activeTips.back().first;
          tipState = result.activeTips.back().second;
        }

        // The surface could be either sensitive or passive
        bool isSensitive = (surface->associatedDetectorElement() != nullptr);
        bool isMaterial = (surface->surfaceMaterial() != nullptr);
        ACTS_VERBOSE("Detected " << (isSensitive ? "sensitive" : "passive")
                                 << " surface: " << surface->geometryId());
        if (isSensitive) {
          // Increment of number of passed sensitive surfaces
          tipState.nSensitiveSurfaces++;
        }
        // Add state if there is already measurement detected on this branch
        if (tipState.nMeasurements > 0 || isMaterial) {
          // New state is to be added. Remove the last tip from active tips now
          if (!result.activeTips.empty()) {
            result.activeTips.erase(result.activeTips.end() - 1);
          }
          // No source links on surface, add either hole or passive material
          // TrackState. No storage allocation for uncalibrated/calibrated
          // measurement and filtered parameter
          auto stateMask =
              TrackStatePropMask::Predicted | TrackStatePropMask::Jacobian;

          // Increment of number of processed states
          tipState.nStates++;
          if (isSensitive) {
            // Increment of number of holes
            tipState.nHoles++;
          }

          // Transport the covariance to a curvilinear surface
          stepper.transportCovarianceToCurvilinear(state.stepping);

          // Update state and stepper with pre material effects
          materialInteractor(surface, state, stepper, navigator,
                             MaterialUpdateStage::PreUpdate);

          // Transport & bind the state to the current surface
          auto boundStateRes =
              stepper.boundState(state.stepping, *surface, false);
          if (!boundStateRes.ok()) {
            return boundStateRes.error();
          }
          auto& boundState = *boundStateRes;
          auto& [boundParams, jacobian, pathLength] = boundState;
          boundParams.covariance() = state.stepping.cov;

          // Add a hole or material track state to the multitrajectory
          std::size_t currentTip = addNonSourcelinkState(
              stateMask, boundState, result, isSensitive, prevTip);
          result.activeTips.emplace_back(currentTip, tipState);

          auto nonSourcelinkState =
              result.fittedStates->getTrackState(currentTip);

          using BranchStopperResult =
              CombinatorialKalmanFilterBranchStopperResult;
          BranchStopperResult branchStopperResult =
              m_extensions.branchStopper(tipState, nonSourcelinkState);

          // Check the branch
          if (branchStopperResult == BranchStopperResult::Continue) {
            // Remembered the active tip and its state
          } else {
            // No branch on this surface
            nBranchesOnSurface = 0;

            if (branchStopperResult == BranchStopperResult::StopAndKeep) {
              storeLastActiveTip(result);
            }

            // Remove the tip from list of active tips
            result.activeTips.erase(result.activeTips.end() - 1);
          }

          // Update state and stepper with post material effects
          materialInteractor(surface, state, stepper, navigator,
                             MaterialUpdateStage::PostUpdate);
        }
      } else {
        // Neither measurement nor material on surface, this branch is still
        // valid. Count the branch on current surface
        nBranchesOnSurface = 1;
      }

      // Reset current tip if there is no branch on current surface
      if (nBranchesOnSurface == 0) {
        ACTS_DEBUG("Branch on surface " << surface->geometryId()
                                        << " is stopped");
        if (!result.activeTips.empty()) {
          ACTS_VERBOSE("Propagation jumps to branch with tip = "
                       << result.activeTips.back().first);
          reset(state, stepper, navigator, result);
        } else {
          ACTS_VERBOSE("Stop Kalman filtering with "
                       << result.lastMeasurementIndices.size()
                       << " found tracks");
          result.finished = true;
        }
      }

      return Result<void>::success();
    }

    /// process new, incompomplete track states and set the filtered state
    /// @note will process the given list of new states, run the updater
    ///     or share the predicted state for states flagged as outliers
    ///     and add them to the list of active tips
    ///
    /// @param gctx The geometry context for this track finding/fitting
    /// @param prevTipState the previous tip state
    /// @param newTrackStateList index list of new track states
    /// @param result which contains among others the new states, and the list of active tips
    /// @return tuple of the number of newly added tips and outlier flag or an error
    Result<std::tuple<unsigned int, bool>> processNewTrackStates(
        const Acts::GeometryContext& gctx, const TipState& prevTipState,
        const boost::container::small_vector<
            typename traj_t::TrackStateProxy::IndexType,
            s_maxBranchesPerSurface>& newTrackStateList,
        result_type& result) const {
      unsigned int nBranchesOnSurface = 0;
      bool isOutlier = false;
      for (typename traj_t::TrackStateProxy::IndexType tipIndex :
           newTrackStateList) {
        // Inherit the tip state from the previous and will be updated
        // later
        typename traj_t::TrackStateProxy trackState(
            result.fittedStates->getTrackState(tipIndex));
        TipState tipState = prevTipState;

        // Increment of number of processedState and passed sensitive surfaces
        tipState.nSensitiveSurfaces++;
        tipState.nStates++;

        using PM = Acts::TrackStatePropMask;
        TrackStateType typeFlags(trackState.typeFlags());
        if (typeFlags.test(TrackStateFlag::OutlierFlag)) {
          // Increment number of outliers
          isOutlier = true;
          tipState.nOutliers++;
          // No Kalman update for outlier
          // Set the filtered parameter index to be the same with predicted
          // parameter
          trackState.shareFrom(PM::Predicted, PM::Filtered);
        } else {
          // Kalman update
          auto updateRes = m_extensions.updater(
              gctx, trackState, Direction::Forward, *updaterLogger);
          if (!updateRes.ok()) {
            ACTS_ERROR("Update step failed: " << updateRes.error());
            return updateRes.error();
          }
          ACTS_VERBOSE(
              "Creating measurement track state with tip = " << tipIndex);
          // Set the measurement flag
          typeFlags.set(TrackStateFlag::MeasurementFlag);
          // Increment number of measurements
          tipState.nMeasurements++;
        }

        // Put tipstate back into active tips to continue with it
        result.activeTips.emplace_back(tipIndex, tipState);

        using BranchStopperResult =
            CombinatorialKalmanFilterBranchStopperResult;
        BranchStopperResult branchStopperResult =
            m_extensions.branchStopper(tipState, trackState);

        // Check if need to stop this branch
        if (branchStopperResult == BranchStopperResult::Continue) {
          // Record the number of branches on surface
          nBranchesOnSurface++;
        } else {
          if (branchStopperResult == BranchStopperResult::StopAndKeep) {
            storeLastActiveTip(result);
          }

          // Remove the tip from list of active tips
          result.activeTips.erase(result.activeTips.end() - 1);
        }
      }
      return std::make_tuple(nBranchesOnSurface, isOutlier);
    }

    /// @brief CombinatorialKalmanFilter actor operation: add a hole or material track state
    ///
    /// @param stateMask The bitmask that instructs which components to allocate
    /// @param boundState The bound state on current surface
    /// @param result is the mutable result state object and which to leave invalid
    /// @param isSensitive The surface is sensitive or passive
    /// @param prevTip The index of the previous state
    ///
    /// @return The tip of added state
    std::size_t addNonSourcelinkState(TrackStatePropMask stateMask,
                                      const BoundState& boundState,
                                      result_type& result, bool isSensitive,
                                      std::size_t prevTip) const {
      // Add a track state
      auto trackStateProxy =
          result.fittedStates->makeTrackState(stateMask, prevTip);
      ACTS_VERBOSE("Create " << (isSensitive ? "Hole" : "Material")
                             << " output track state #"
                             << trackStateProxy.index()
                             << " with mask: " << stateMask);

      const auto& [boundParams, jacobian, pathLength] = boundState;
      // Fill the track state
      trackStateProxy.predicted() = boundParams.parameters();
      trackStateProxy.predictedCovariance() = boundParams.covariance().value();
      trackStateProxy.jacobian() = jacobian;
      trackStateProxy.pathLength() = pathLength;
      // Set the surface
      trackStateProxy.setReferenceSurface(
          boundParams.referenceSurface().getSharedPtr());

      // Set the track state flags
      auto typeFlags = trackStateProxy.typeFlags();
      if (trackStateProxy.referenceSurface().surfaceMaterial() != nullptr) {
        typeFlags.set(TrackStateFlag::MaterialFlag);
      }
      typeFlags.set(TrackStateFlag::ParameterFlag);
      if (isSensitive) {
        typeFlags.set(TrackStateFlag::HoleFlag);
      }

      // Set the filtered parameter index to be the same with predicted
      // parameter
      trackStateProxy.shareFrom(TrackStatePropMask::Predicted,
                                TrackStatePropMask::Filtered);

      return trackStateProxy.index();
    }

    /// @brief CombinatorialKalmanFilter actor operation: material interaction
    ///
    /// @tparam propagator_state_t is the type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param surface The surface where the material interaction happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param updateStage The material update stage
    ///
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void materialInteractor(const Surface* surface, propagator_state_t& state,
                            const stepper_t& stepper,
                            const navigator_t& navigator,
                            const MaterialUpdateStage& updateStage) const {
      if (surface == nullptr) {
        return;
      }

      // Indicator if having material
      bool hasMaterial = false;

      if (surface->surfaceMaterial() != nullptr) {
        // Prepare relevant input particle properties
        detail::PointwiseMaterialInteraction interaction(surface, state,
                                                         stepper);
        // Evaluate the material properties
        if (interaction.evaluateMaterialSlab(state, navigator, updateStage)) {
          // Surface has material at this stage
          hasMaterial = true;

          // Evaluate the material effects
          interaction.evaluatePointwiseMaterialInteraction(multipleScattering,
                                                           energyLoss);

          // Screen out material effects info
          ACTS_VERBOSE("Material effects on surface: "
                       << surface->geometryId()
                       << " at update stage: " << updateStage << " are :");
          ACTS_VERBOSE("eLoss = "
                       << interaction.Eloss * interaction.navDir << ", "
                       << "variancePhi = " << interaction.variancePhi << ", "
                       << "varianceTheta = " << interaction.varianceTheta
                       << ", "
                       << "varianceQoverP = " << interaction.varianceQoverP);

          // Update the state and stepper with material effects
          interaction.updateState(state, stepper, addNoise);
        }
      }

      if (!hasMaterial) {
        // Screen out message
        ACTS_VERBOSE("No material effects on surface: " << surface->geometryId()
                                                        << " at update stage: "
                                                        << updateStage);
      }
    }

    void storeLastActiveTip(result_type& result) const {
      const auto& [currentTip, tipState] = result.activeTips.back();
      // @TODO: Keep information on tip state around so we don't have to
      //        recalculate it later
      ACTS_VERBOSE("Find track with entry index = "
                   << currentTip << " and there are nMeasurements = "
                   << tipState.nMeasurements
                   << ", nOutliers = " << tipState.nOutliers
                   << ", nHoles = " << tipState.nHoles << " on track");
      result.lastTrackIndices.emplace_back(currentTip);

      std::optional<MultiTrajectoryTraits::IndexType>
          lastMeasurementIndexCandidate;
      result.fittedStates->visitBackwards(
          currentTip, [&](const auto& trackState) {
            bool isMeasurement =
                trackState.typeFlags().test(TrackStateFlag::MeasurementFlag);
            if (isMeasurement) {
              lastMeasurementIndexCandidate = trackState.index();
              return false;
            }
            return true;
          });
      if (lastMeasurementIndexCandidate.has_value()) {
        ACTS_VERBOSE("Last measurement found on track with entry index = "
                     << currentTip << " and measurement index = "
                     << lastMeasurementIndexCandidate.value());
        result.lastMeasurementIndices.emplace_back(
            lastMeasurementIndexCandidate.value());
      } else {
        ACTS_VERBOSE(
            "No measurement found on track with entry index = " << currentTip);
      }
    }

    CombinatorialKalmanFilterExtensions<traj_t> m_extensions;

    /// The source link accessor
    source_link_accessor_t m_sourcelinkAccessor;

    using source_link_iterator_t =
        decltype(std::declval<decltype(m_sourcelinkAccessor(
                     *static_cast<const Surface*>(nullptr)))>()
                     .first);

    using TrackStateCandidateCreator =
        typename CombinatorialKalmanFilterOptions<
            source_link_iterator_t, traj_t>::TrackStateCandidateCreator;

    /// the stateCandidator to be used
    /// @note will be set to a default trackStateCandidateCreator or the one
    //        provided via the extension
    TrackStateCandidateCreator trackStateCandidateCreator;

    /// End of world aborter
    EndOfWorldReached endOfWorldReached;

    /// Actor logger instance
    const Logger* actorLogger{nullptr};
    /// Updater logger instance
    const Logger* updaterLogger{nullptr};

    const Logger& logger() const { return *actorLogger; }
  };

  template <typename source_link_accessor_t, typename parameters_t>
  class Aborter {
   public:
    /// Broadcast the action type
    using action_type = Actor<source_link_accessor_t, parameters_t>;

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t, typename result_t>
    bool operator()(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, const result_t& result,
                    const Logger& /*logger*/) const {
      if (result.finished) {
        return true;
      }
      return false;
    }
  };

  /// Void path limit reached aborter to replace the default since the path
  /// limit is handled in the CKF actor internally.
  struct StubPathLimitReached {
    double internalLimit{};

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    bool operator()(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/,
                    const Logger& /*logger*/) const {
      return false;
    }
  };

 public:
  /// Combinatorial Kalman Filter implementation, calls the Kalman filter
  ///
  /// @tparam source_link_iterator_t Type of the source link iterator
  /// @tparam start_parameters_container_t Type of the initial parameters
  ///                                      container
  /// @tparam calibrator_t Type of the source link calibrator
  /// @tparam measurement_selector_t Type of the measurement selector
  /// @tparam track_container_t Type of the track container backend
  /// @tparam holder_t Type defining track container backend ownership
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param initialParameters The initial track parameters
  /// @param tfOptions CombinatorialKalmanFilterOptions steering the track
  ///                  finding
  /// @param trackContainer Input track container to use
  /// @note The input measurements are given in the form of @c SourceLinks.
  ///       It's @c calibrator_t's job to turn them into calibrated measurements
  ///       used in the track finding.
  ///
  /// @return a container of track finding result for all the initial track
  /// parameters
  template <typename source_link_iterator_t, typename start_parameters_t,
            typename track_container_t, template <typename> class holder_t,
            typename parameters_t = BoundTrackParameters>
  auto findTracks(
      const start_parameters_t& initialParameters,
      const CombinatorialKalmanFilterOptions<source_link_iterator_t, traj_t>&
          tfOptions,
      TrackContainer<track_container_t, traj_t, holder_t>& trackContainer) const
      -> Result<std::vector<
          typename std::decay_t<decltype(trackContainer)>::TrackProxy>> {
    using TrackContainer = typename std::decay_t<decltype(trackContainer)>;
    using SourceLinkAccessor =
        SourceLinkAccessorDelegate<source_link_iterator_t>;

    // Create the ActionList and AbortList
    using CombinatorialKalmanFilterAborter =
        Aborter<SourceLinkAccessor, parameters_t>;
    using CombinatorialKalmanFilterActor =
        Actor<SourceLinkAccessor, parameters_t>;
    using Actors = ActionList<CombinatorialKalmanFilterActor>;
    using Aborters = AbortList<CombinatorialKalmanFilterAborter>;

    // Create relevant options for the propagation options
    PropagatorOptions<Actors, Aborters> propOptions(tfOptions.geoContext,
                                                    tfOptions.magFieldContext);

    // Set the trivial propagator options
    propOptions.setPlainOptions(tfOptions.propagatorPlainOptions);

    // Catch the actor
    auto& combKalmanActor =
        propOptions.actionList.template get<CombinatorialKalmanFilterActor>();
    combKalmanActor.targetReached.surface = tfOptions.targetSurface;
    combKalmanActor.multipleScattering = tfOptions.multipleScattering;
    combKalmanActor.energyLoss = tfOptions.energyLoss;
    combKalmanActor.actorLogger = m_actorLogger.get();
    combKalmanActor.updaterLogger = m_updaterLogger.get();
    combKalmanActor.calibrationContextPtr = &tfOptions.calibrationContext.get();

    // copy source link accessor, calibrator and measurement selector
    combKalmanActor.m_sourcelinkAccessor = tfOptions.sourcelinkAccessor;
    combKalmanActor.m_extensions = tfOptions.extensions;
    combKalmanActor.trackStateCandidateCreator =
        tfOptions.trackStateCandidateCreator;
    DefaultTrackStateCreator defaultTrackStateCreator;
    // connect a default state candidate creator if no state candidate creator
    // was provided via the extension
    if (!combKalmanActor.trackStateCandidateCreator.connected()) {
      defaultTrackStateCreator.calibrator = tfOptions.extensions.calibrator;
      defaultTrackStateCreator.measurementSelector =
          tfOptions.extensions.measurementSelector;
      combKalmanActor.trackStateCandidateCreator.template connect<
          &DefaultTrackStateCreator::template createSourceLinkTrackStates<
              source_link_iterator_t>>(&defaultTrackStateCreator);
    }

    auto propState =
        m_propagator.template makeState(initialParameters, propOptions);

    auto& r = propState.template get<CombinatorialKalmanFilterResult<traj_t>>();
    r.fittedStates = &trackContainer.trackStateContainer();
    r.stateBuffer = std::make_shared<traj_t>();

    auto propagationResult = m_propagator.propagate(propState);

    auto result = m_propagator.makeResult(
        std::move(propState), propagationResult, propOptions, false);

    if (!result.ok()) {
      ACTS_ERROR("Propagation failed: " << result.error() << " "
                                        << result.error().message()
                                        << " with the initial parameters: \n"
                                        << initialParameters.parameters());
      return result.error();
    }

    auto& propRes = *result;

    /// Get the result of the CombinatorialKalmanFilter
    auto combKalmanResult = std::move(
        propRes.template get<CombinatorialKalmanFilterResult<traj_t>>());

    /// The propagation could already reach max step size
    /// before the track finding is finished during two phases:
    // -> filtering for track finding;
    // -> surface targeting to get fitted parameters at target surface.
    // This is regarded as a failure.
    // @TODO: Implement distinguishment between the above two cases if
    // necessary
    if (combKalmanResult.lastError.ok() && !combKalmanResult.finished) {
      combKalmanResult.lastError = Result<void>(
          CombinatorialKalmanFilterError::PropagationReachesMaxSteps);
    }

    if (!combKalmanResult.lastError.ok()) {
      ACTS_ERROR("CombinatorialKalmanFilter failed: "
                 << combKalmanResult.lastError.error() << " "
                 << combKalmanResult.lastError.error().message()
                 << " with the initial parameters: \n"
                 << initialParameters.parameters());
    }

    std::vector<typename TrackContainer::TrackProxy> tracks;
    tracks.reserve(combKalmanResult.lastMeasurementIndices.size());

    for (auto tip : combKalmanResult.lastMeasurementIndices) {
      auto track = trackContainer.makeTrack();
      track.tipIndex() = tip;

      // Set fitted track parameters if available. This will only be the case if
      // a target surface is set. Without a target surface there cannot be
      // fitted parameters and the user will have to extrapolate the track to a
      // target surface themselves.
      if (auto it = combKalmanResult.fittedParameters.find(tip);
          it != combKalmanResult.fittedParameters.end()) {
        const BoundTrackParameters& parameters = it->second;
        track.parameters() = parameters.parameters();
        track.covariance() = *parameters.covariance();
        track.setReferenceSurface(parameters.referenceSurface().getSharedPtr());
      }

      calculateTrackQuantities(track);

      tracks.push_back(std::move(track));
    }

    return tracks;
  }
};

}  // namespace Acts
