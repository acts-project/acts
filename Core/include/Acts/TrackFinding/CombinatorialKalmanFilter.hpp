// This file is part of the Acts project.
//
// Copyright (C) 2016-2021 CERN for the benefit of the Acts project
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
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterError.hpp"
#include "Acts/TrackFinding/SourceLinkAccessorConcept.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/VoidKalmanComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <functional>
#include <memory>
#include <unordered_map>

namespace Acts {

/// Track quality summary for one trajectory.
///
/// This could be used to decide if a track is to be recorded when the
/// filtering is done or to be terminated due to its bad quality
/// @todo: add other useful info, e.g. chi2
struct CombinatorialKalmanFilterTipState {
  // Number of passed sensitive surfaces
  size_t nSensitiveSurfaces = 0;
  // Number of track states
  size_t nStates = 0;
  // Number of (non-outlier) measurements
  size_t nMeasurements = 0;
  // Number of outliers
  size_t nOutliers = 0;
  // Number of holes
  size_t nHoles = 0;
};

/// Extension struct which holds the delegates to customize the CKF behavior
template <typename traj_t>
struct CombinatorialKalmanFilterExtensions {
  using candidate_container_t =
      typename std::vector<typename MultiTrajectory<traj_t>::TrackStateProxy>;
  using MeasurementSelector =
      Delegate<Result<std::pair<typename candidate_container_t::iterator,
                                typename candidate_container_t::iterator>>(
          candidate_container_t& trackStates, bool&, const Logger&)>;
  using BranchStopper =
      Delegate<bool(const CombinatorialKalmanFilterTipState&)>;

  /// The Calibrator is a dedicated calibration algorithm that allows
  /// to calibrate measurements using track information, this could be
  /// e.g. sagging for wires, module deformations, etc.
  typename KalmanFitterExtensions<traj_t>::Calibrator calibrator;

  /// The updater incorporates measurement information into the track parameters
  typename KalmanFitterExtensions<traj_t>::Updater updater;

  /// The smoother back-propagates measurement information along the track
  typename KalmanFitterExtensions<traj_t>::Smoother smoother;

  /// The measurement selector is called during the filtering by the Actor.
  MeasurementSelector measurementSelector;

  BranchStopper branchStopper;

  /// Default constructor which connects the default void components
  CombinatorialKalmanFilterExtensions() {
    calibrator.template connect<&voidKalmanCalibrator<traj_t>>();
    updater.template connect<&voidKalmanUpdater<traj_t>>();
    smoother.template connect<&voidKalmanSmoother<traj_t>>();
    branchStopper.connect<voidBranchStopper>();
    measurementSelector.template connect<voidMeasurementSelector>();
  }

 private:
  /// Default measurement selector which will return all measurements
  /// @param candidates Measurement track state candidates
  /// @param isOutlier Output variable indicating whether the returned state is an outlier (unused)
  /// @param logger A logger instance
  static Result<std::pair<typename std::vector<typename MultiTrajectory<
                              traj_t>::TrackStateProxy>::iterator,
                          typename std::vector<typename MultiTrajectory<
                              traj_t>::TrackStateProxy>::iterator>>
  voidMeasurementSelector(
      typename std::vector<typename MultiTrajectory<traj_t>::TrackStateProxy>&
          candidates,
      bool& isOutlier, const Logger& logger) {
    (void)isOutlier;
    (void)logger;
    return std::pair{candidates.begin(), candidates.end()};
  };

  /// Default branch stopper which will never stop
  /// @param tipState The tip state to decide whether to stop (unused)
  /// @return false
  static bool voidBranchStopper(
      const CombinatorialKalmanFilterTipState& tipState) {
    (void)tipState;
    return false;
  }
};

/// Delegate type that retrieves a range of source links to for a given surface
/// to be processed by the CKF
template <typename source_link_iterator_t>
using SourceLinkAccessorDelegate =
    Delegate<std::pair<source_link_iterator_t, source_link_iterator_t>(
        const Surface&)>;

/// Combined options for the combinatorial Kalman filter.
///
/// @tparam source_link_accessor_t Source link accessor type, should be
/// semiregular.
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
  /// @param rSurface The reference surface for the eventual track fitting to be
  /// expressed at
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  /// @param rSmoothing Whether to run smoothing to get fitted parameter
  CombinatorialKalmanFilterOptions(
      const GeometryContext& gctx, const MagneticFieldContext& mctx,
      std::reference_wrapper<const CalibrationContext> cctx,
      SourceLinkAccessor accessor_,
      CombinatorialKalmanFilterExtensions<traj_t> extensions_,
      const PropagatorPlainOptions& pOptions, const Surface* rSurface = nullptr,
      bool mScattering = true, bool eLoss = true, bool rSmoothing = true)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        sourcelinkAccessor(std::move(accessor_)),
        extensions(extensions_),
        propagatorPlainOptions(pOptions),
        referenceSurface(rSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        smoothing(rSmoothing) {}

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

  /// The reference Surface
  const Surface* referenceSurface = nullptr;

  /// Whether to consider multiple scattering.
  bool multipleScattering = true;

  /// Whether to consider energy loss.
  bool energyLoss = true;

  /// Whether to run smoothing to get fitted parameter
  bool smoothing = true;
};

template <typename traj_t>
struct CombinatorialKalmanFilterResult {
  // Fitted states that the actor has handled.
  traj_t* fittedStates{nullptr};

  // These is used internally to store candidate trackstates
  std::shared_ptr<traj_t> stateBuffer;
  std::vector<typename MultiTrajectory<traj_t>::TrackStateProxy>
      trackStateCandidates;

  // This is the indices of the 'tip' of the tracks stored in multitrajectory.
  // This correspond to the last measurment state in the multitrajectory.
  std::vector<MultiTrajectoryTraits::IndexType> lastMeasurementIndices;

  // This is the indices of the 'tip' of the tracks stored in multitrajectory.
  // This correspond to the last state in the multitrajectory.
  std::vector<MultiTrajectoryTraits::IndexType> lastTrackIndices;

  // The Parameters at the provided surface for separate tracks
  std::unordered_map<MultiTrajectoryTraits::IndexType, BoundTrackParameters>
      fittedParameters;

  // The indices of the 'tip' of the unfinished tracks
  std::vector<std::pair<MultiTrajectoryTraits::IndexType,
                        CombinatorialKalmanFilterTipState>>
      activeTips;

  // The indices of track states and corresponding source links on different
  // surfaces
  std::unordered_map<const Surface*, std::unordered_map<size_t, size_t>>
      sourcelinkTips;

  // Indicator if filtering has been done
  bool filtered = false;

  // Indicator if smoothing has been done.
  bool smoothed = false;

  // The index for the current smoothing track
  MultiTrajectoryTraits::IndexType iSmoothed = 0;

  // Indicator if track finding has been done
  bool finished = false;

  Result<void> result{Result<void>::success()};

  // TODO place into options and make them accessible?
  AbortList<PathLimitReached, EndOfWorldReached, ParticleStopped> abortList;
};

/// Combinatorial Kalman filter to find tracks.
///
///
/// @tparam propagator_t Type of the propagator
///
/// The CombinatorialKalmanFilter contains an Actor and a Sequencer sub-class.
/// The Sequencer has to be part of the Navigator of the Propagator
/// in order to initialize and provide the measurement surfaces.
///
/// The Actor is part of the Propagation call and does the Kalman update
/// and eventually the smoothing.  Updater, Smoother and Calibrator are
/// given to the Actor for further use:
/// - The Updater is the implemented kalman updater formalism, it
///   runs via a visitor pattern through the measurements.
/// - The Smoother is called at the end of the filtering (track finding) by the
/// Actor.
///
/// Measurements are not required to be ordered for the
/// CombinatorialKalmanFilter, measurement ordering needs to be figured out by
/// the navigation of the propagator.
///
/// The void components are provided mainly for unit testing.
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
        m_actorLogger{m_logger->cloneWithSuffix("Actor")} {}

 private:
  using KalmanNavigator = typename propagator_t::Navigator;

  /// The propgator for the transport and material update
  propagator_t m_propagator;

  std::unique_ptr<const Logger> m_logger;
  std::shared_ptr<const Logger> m_actorLogger;

  const Logger& logger() const { return *m_logger; }

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
    // The source link container type
    /// Broadcast the result_type
    using result_type = CombinatorialKalmanFilterResult<traj_t>;

    /// The target surface
    const Surface* targetSurface = nullptr;

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// Whether to consider energy loss.
    bool energyLoss = true;

    /// Whether to run smoothing to get fitted parameter
    bool smoothing = true;

    /// @brief CombinatorialKalmanFilter actor operation
    ///
    /// @tparam propagator_state_t Type of the Propagagor state
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

      // Update:
      // - Waiting for a current surface
      auto surface = navigator.currentSurface(state.navigation);
      if (surface != nullptr and not result.filtered) {
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
          result.result = res.error();
        }
      }

      // Reset propagation state:
      // - When navigation breaks and there is stil active tip present after
      // recording&removing track tips on current surface
      if (navigator.navigationBreak(state.navigation) and not result.filtered) {
        // Record the tips on current surface as trajectory entry indices
        // (taking advantage of fact that those tips are consecutive in list of
        // active tips) and remove those tips from active tips
        if (not result.activeTips.empty()) {
          // The last active tip
          const auto& lastActiveTip = result.activeTips.back().first;
          // Get the index of previous state
          const auto& iprevious =
              result.fittedStates->getTrackState(lastActiveTip).previous();
          // Find the track states which have the same previous state and remove
          // them from active tips
          while (not result.activeTips.empty()) {
            const auto& [currentTip, tipState] = result.activeTips.back();
            if (result.fittedStates->getTrackState(currentTip).previous() !=
                iprevious) {
              break;
            }
            // Record the tips if there are measurements on the track
            if (tipState.nMeasurements > 0) {
              ACTS_VERBOSE("Find track with entry index = "
                           << currentTip << " and there are nMeasurements = "
                           << tipState.nMeasurements
                           << ", nOutliers = " << tipState.nOutliers
                           << ", nHoles = " << tipState.nHoles << " on track");
              result.lastTrackIndices.emplace_back(currentTip);
              // Set the lastMeasurementIndex to the last measurement
              // to ignore the states after it in the rest of the algorithm
              auto lastMeasurementIndex = currentTip;
              auto lastMeasurementState =
                  result.fittedStates->getTrackState(lastMeasurementIndex);
              bool isMeasurement = lastMeasurementState.typeFlags().test(
                  TrackStateFlag::MeasurementFlag);
              while (!isMeasurement) {
                lastMeasurementIndex = lastMeasurementState.previous();
                lastMeasurementState =
                    result.fittedStates->getTrackState(lastMeasurementIndex);
                isMeasurement = lastMeasurementState.typeFlags().test(
                    TrackStateFlag::MeasurementFlag);
              }
              result.lastMeasurementIndices.emplace_back(lastMeasurementIndex);
              // @TODO: Keep information on tip state around so we don't have to
              //        recalculate it later
            }
            // Remove the tip from list of active tips
            result.activeTips.erase(result.activeTips.end() - 1);
          }
        }
        // If no more active tip, done with filtering; Otherwise, reset
        // propagation state to track state at last tip of active tips
        if (result.activeTips.empty()) {
          ACTS_VERBOSE("Kalman filtering finds "
                       << result.lastTrackIndices.size() << " tracks");
          result.filtered = true;
        } else {
          ACTS_VERBOSE("Propagation jumps to branch with tip = "
                       << result.activeTips.back().first);
          reset(state, stepper, navigator, result);
        }
      }

      if (result.abortList(state, stepper, navigator, result, logger())) {
        navigator.targetReached(state.navigation, false);
        if (result.activeTips.empty()) {
          // we are already done
        } else if (result.activeTips.size() == 1) {
          // this was the last track - we are done
          ACTS_VERBOSE("Kalman filtering finds "
                       << result.lastTrackIndices.size() << " tracks");
          result.filtered = true;
        } else {
          // remove the active tip and continue with the next
          result.activeTips.erase(result.activeTips.end() - 1);
          reset(state, stepper, navigator, result);
        }
      }

      // Post-processing after filtering phase
      if (result.filtered) {
        // Return error if filtering finds no tracks
        if (result.lastTrackIndices.empty()) {
          ACTS_WARNING("No tracks found");
          result.finished = true;
        } else {
          if (not smoothing) {
            ACTS_VERBOSE("Finish Kalman filtering");
            // Remember that track finding is done
            result.finished = true;
          } else {
            // Iterate over the found tracks for smoothing and getting the
            // fitted parameter. This needs to be accomplished in different
            // propagation steps:
            // -> first run smoothing for found track indexed with iSmoothed
            if (not result.smoothed) {
              ACTS_VERBOSE(
                  "Finalize/run smoothing for track with last measurement "
                  "index = "
                  << result.lastMeasurementIndices.at(result.iSmoothed));
              // --> Search the starting state to run the smoothing
              // --> Call the smoothing
              // --> Set a stop condition when all track states have been
              // handled
              auto res = finalize(state, stepper, result);
              if (!res.ok()) {
                ACTS_ERROR("Error in finalize: " << res.error());
                result.result = res.error();
              }
              result.smoothed = true;
            }
            // -> then progress to target/reference surface and built the final
            // track parameters for found track indexed with iSmoothed
            if (result.smoothed and targetReached(state, stepper, navigator,
                                                  *targetSurface, logger())) {
              ACTS_VERBOSE(
                  "Completing the track with last measurement index = "
                  << result.lastMeasurementIndices.at(result.iSmoothed));
              // Transport & bind the parameter to the final surface
              auto res = stepper.boundState(state.stepping, *targetSurface);
              if (!res.ok()) {
                ACTS_ERROR("Error in finalize: " << res.error());
                result.result = res.error();
                return;
              }

              auto fittedState = *res;
              // Assign the fitted parameters
              result.fittedParameters.emplace(
                  result.lastMeasurementIndices.at(result.iSmoothed),
                  std::get<BoundTrackParameters>(fittedState));
              // If there are more trajectories to handle:
              // -> set the targetReached status to false
              // -> set the smoothed status to false
              // -> update the index of track to be smoothed
              if (result.iSmoothed < result.lastMeasurementIndices.size() - 1) {
                navigator.targetReached(state.navigation, false);
                result.smoothed = false;
                result.iSmoothed++;
                // Reverse navigation direction to start targeting for the rest
                // tracks
                state.stepping.navDir =
                    (state.stepping.navDir == NavigationDirection::Backward)
                        ? NavigationDirection::Forward
                        : NavigationDirection::Backward;
                // To avoid meaningless navigation target call
                state.stepping.stepSize =
                    ConstrainedStep(state.stepping.navDir *
                                    std::abs(state.options.maxStepSize));
              } else {
                ACTS_VERBOSE("Finish Kalman filtering and smoothing");
                // Remember that track finding is done
                result.finished = true;
              }
            }
          }  // if run smoothing
        }    // if there are found tracks
      }      // if filtering is done
    }

    /// @brief Kalman actor operation : reset propagation
    ///
    /// @tparam propagator_state_t Type of Propagagor state
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

      // Update the stepping state
      stepper.resetState(state.stepping, currentState.filtered(),
                         currentState.filteredCovariance(),
                         currentState.referenceSurface(), state.stepping.navDir,
                         state.options.maxStepSize);

      // Reset the navigation state
      // Set targetSurface to nullptr for forward filtering; it's only needed
      // after smoothing
      navigator.resetState(
          state.navigation, state.geoContext, stepper.position(state.stepping),
          stepper.direction(state.stepping), state.stepping.navDir,
          &currentState.referenceSurface(), nullptr);

      // No Kalman filtering for the starting surface, but still need
      // to consider the material effects here
      materialInteractor(navigator.currentSurface(state.navigation), state,
                         stepper, navigator, MaterialUpdateStage::FullUpdate);

      detail::setupLoopProtection(
          state, stepper, result.abortList.template get<PathLimitReached>(),
          logger());
    }

    /// @brief CombinatorialKalmanFilter actor operation :
    /// - filtering for all measurement(s) on surface
    /// - store selected track states in multiTrajectory
    /// - update propagator state to the (last) selected track state
    ///
    /// @tparam propagator_state_t Type of the Propagagor state
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
      // Initialize the number of branches on current surface
      size_t nBranchesOnSurface = 0;

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
        auto boundState = *boundStateRes;

        // Retrieve the previous tip and its state
        // The states created on this surface will have the common previous tip
        size_t prevTip = SIZE_MAX;
        TipState prevTipState;
        if (not result.activeTips.empty()) {
          prevTip = result.activeTips.back().first;
          prevTipState = result.activeTips.back().second;
          // New state is to be added. Remove the last tip from active tips
          result.activeTips.erase(result.activeTips.end() - 1);
        }

        // Create trackstates for all source links (will be filtered later)
        // Results are stored in result => no return value
        createSourceLinkTrackStates(state.geoContext, result, boundState,
                                    prevTip, slBegin, slEnd);

        // Invoke the measurement selector to select compatible measurements
        // with the predicted track parameter.
        // It can modify the trackStateCandidates vector, and will return a pair
        // of iterators marking the range of accepted measurements (track
        // states)
        bool isOutlier = false;
        auto selectorResult = m_extensions.measurementSelector(
            result.trackStateCandidates, isOutlier, logger());

        if (!selectorResult.ok()) {
          ACTS_ERROR("Selection of calibrated measurements failed: "
                     << selectorResult.error());
          return selectorResult.error();
        }
        auto selectedTrackStateRange = *selectorResult;

        auto procRes = processSelectedTrackStates(
            state.geoContext, selectedTrackStateRange.first,
            selectedTrackStateRange.second, result, isOutlier, prevTipState,
            nBranchesOnSurface);

        if (!procRes.ok()) {
          ACTS_ERROR(
              "Processing of selected track states failed: " << procRes.error())
          return procRes.error();
        }

        if (nBranchesOnSurface > 0 and not isOutlier) {
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
        size_t prevTip = SIZE_MAX;
        TipState tipState;
        if (not result.activeTips.empty()) {
          prevTip = result.activeTips.back().first;
          tipState = result.activeTips.back().second;
        }

        // The surface could be either sensitive or passive
        bool isSensitive = (surface->associatedDetectorElement() != nullptr);
        bool isMaterial = (surface->surfaceMaterial() != nullptr);
        std::string type = isSensitive ? "sensitive" : "passive";
        ACTS_VERBOSE("Detected " << type
                                 << " surface: " << surface->geometryId());
        if (isSensitive) {
          // Increment of number of passed sensitive surfaces
          tipState.nSensitiveSurfaces++;
        }
        // Add state if there is already measurement detected on this branch
        // For in-sensitive surface, only add state when smoothing is
        // required
        bool createState = false;
        if (smoothing) {
          createState = (tipState.nMeasurements > 0 or isMaterial);
        } else {
          createState = (tipState.nMeasurements > 0 and isSensitive);
        }
        if (createState) {
          // New state is to be added. Remove the last tip from active tips now
          if (not result.activeTips.empty()) {
            result.activeTips.erase(result.activeTips.end() - 1);
          }
          // No source links on surface, add either hole or passive material
          // TrackState. No storage allocation for uncalibrated/calibrated
          // measurement and filtered parameter
          auto stateMask =
              ~(TrackStatePropMask::Calibrated | TrackStatePropMask::Filtered);

          // Increment of number of processed states
          tipState.nStates++;
          size_t currentTip = SIZE_MAX;
          if (isSensitive) {
            // Incremet of number of holes
            tipState.nHoles++;
          }

          // Transport & bind the state to the current surface
          auto res = stepper.boundState(state.stepping, *surface);
          if (!res.ok()) {
            ACTS_ERROR("Error in filter: " << res.error());
            return res.error();
          }
          const auto boundState = *res;
          // Add a hole or material track state to the multitrajectory
          currentTip = addNonSourcelinkState(stateMask, boundState, result,
                                             isSensitive, prevTip);

          // Check the branch
          if (not m_extensions.branchStopper(tipState)) {
            // Remember the active tip and its state
            result.activeTips.emplace_back(currentTip, tipState);
          } else {
            // No branch on this surface
            nBranchesOnSurface = 0;
          }
        }
        if (surface->surfaceMaterial() != nullptr) {
          // Update state and stepper with material effects
          materialInteractor(surface, state, stepper, navigator,
                             MaterialUpdateStage::FullUpdate);
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
        if (not result.activeTips.empty()) {
          ACTS_VERBOSE("Propagation jumps to branch with tip = "
                       << result.activeTips.back().first);
          reset(state, stepper, navigator, result);
        } else {
          ACTS_VERBOSE("Stop Kalman filtering with "
                       << result.lastMeasurementIndices.size()
                       << " found tracks");
          result.filtered = true;
        }
      }

      return Result<void>::success();
    }

    /// Create and fill track states for all source links
    /// @param gctx The current geometry context
    /// @param result Reference to the result struct of the actor
    /// @param boundState Bound state from the propagation on this surface
    /// @param prevTip Index pointing at previous trajectory state (i.e. tip)
    /// @param slBegin Begin iterator for sourcelinks
    /// @param slEnd End iterator for sourcelinks
    template <typename source_link_iterator_t>
    void createSourceLinkTrackStates(const Acts::GeometryContext& gctx,
                                     result_type& result,
                                     const BoundState& boundState,
                                     size_t prevTip,
                                     source_link_iterator_t slBegin,
                                     source_link_iterator_t slEnd) const {
      const auto& [boundParams, jacobian, pathLength] = boundState;

      result.trackStateCandidates.clear();
      if constexpr (std::is_same_v<
                        typename std::iterator_traits<
                            source_link_iterator_t>::iterator_category,
                        std::random_access_iterator_tag>) {
        result.trackStateCandidates.reserve(std::distance(slBegin, slEnd));
      }

      result.stateBuffer->clear();

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

        size_t tsi = result.stateBuffer->addTrackState(mask, prevTip);
        // CAREFUL! This trackstate has a previous index that is not in this
        // MultiTrajectory Visiting brackwards from this track state will
        // fail!
        auto ts = result.stateBuffer->getTrackState(tsi);

        if (it == slBegin) {
          // only set these for first
          ts.predicted() = boundParams.parameters();
          if (boundParams.covariance()) {
            ts.predictedCovariance() = *boundParams.covariance();
          }
          ts.jacobian() = jacobian;
        } else {
          // subsequent track states can reuse
          auto& first = result.trackStateCandidates.front();
          ts.shareFrom(first, PM::Predicted);
          ts.shareFrom(first, PM::Jacobian);
        }

        ts.pathLength() = pathLength;

        ts.setReferenceSurface(boundParams.referenceSurface().getSharedPtr());

        ts.setUncalibratedSourceLink(sourceLink);

        // now calibrate the track state
        m_extensions.calibrator(gctx, ts);

        result.trackStateCandidates.push_back(ts);
      }
    }

    /// Handle the list of selected track states
    /// @param gctx The current geometry context
    /// @param begin The start iterator for selected track states
    /// @param end The end iterator for selected track states
    /// @param result Reference to the actor result struct
    /// @param isOutlier If this track state is a single outlier one
    /// @param prevTipState Tip state prior to this surface
    /// @param [in,out] nBranchesOnSurface Number of branches on surface, will be updated
    Result<void> processSelectedTrackStates(
        const Acts::GeometryContext& gctx,
        typename std::vector<
            typename MultiTrajectory<traj_t>::TrackStateProxy>::const_iterator
            begin,
        typename std::vector<
            typename MultiTrajectory<traj_t>::TrackStateProxy>::const_iterator
            end,
        result_type& result, bool isOutlier, const TipState& prevTipState,
        size_t& nBranchesOnSurface) const {
      using PM = TrackStatePropMask;

      std::optional<typename MultiTrajectory<traj_t>::TrackStateProxy>
          firstTrackState{std::nullopt};
      for (auto it = begin; it != end; ++it) {
        auto& candidateTrackState = *it;

        PM mask = PM::All;

        if (it != begin) {
          // subsequent track states don't need storage for these
          mask = ~PM::Predicted & ~PM::Jacobian;
        }

        if (isOutlier) {
          mask &= ~PM::Filtered;  // outlier won't have separate filtered
                                  // parameters
        }

        // copy this trackstate into fitted states MultiTrajectory
        typename MultiTrajectory<traj_t>::TrackStateProxy trackState =
            result.fittedStates->getTrackState(
                result.fittedStates->addTrackState(
                    mask, candidateTrackState.previous()));

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

        auto& typeFlags = trackState.typeFlags();
        if (trackState.referenceSurface().surfaceMaterial() != nullptr) {
          typeFlags.set(TrackStateFlag::MaterialFlag);
        }
        typeFlags.set(TrackStateFlag::ParameterFlag);

        // Inherit the tip state from the previous and will be updated
        // later
        TipState tipState = prevTipState;
        size_t currentTip = trackState.index();

        // Increment of number of processedState and passed sensitive surfaces
        tipState.nSensitiveSurfaces++;
        tipState.nStates++;

        if (isOutlier) {
          ACTS_VERBOSE(
              "Creating outlier track state with tip = " << currentTip);
          // Set the outlier flag
          typeFlags.set(TrackStateFlag::OutlierFlag);
          // Increment number of outliers
          tipState.nOutliers++;
          // No Kalman update for outlier
          // Set the filtered parameter index to be the same with predicted
          // parameter
          trackState.shareFrom(PM::Predicted, PM::Filtered);

        } else {
          // Kalman update
          auto updateRes = m_extensions.updater(
              gctx, trackState, NavigationDirection::Forward, getDummyLogger());
          if (!updateRes.ok()) {
            ACTS_ERROR("Update step failed: " << updateRes.error());
            return updateRes.error();
          }
          ACTS_VERBOSE(
              "Creating measurement track state with tip = " << currentTip);
          // Set the measurement flag
          typeFlags.set(TrackStateFlag::MeasurementFlag);
          // Increment number of measurements
          tipState.nMeasurements++;
        }

        // Check if need to stop this branch
        if (not m_extensions.branchStopper(tipState)) {
          // Put tipstate back into active tips to continue with it
          result.activeTips.emplace_back(currentTip, tipState);
          // Record the number of branches on surface
          nBranchesOnSurface++;
        }
      }
      return Result<void>::success();
    }

    /// @brief CombinatorialKalmanFilter actor operation : add hole or material track state
    ///
    /// @param stateMask The bitmask that instructs which components to allocate
    /// @param boundState The bound state on current surface
    /// @param result is the mutable result state object
    /// and which to leave invalid
    /// @param isSensitive The surface is sensitive or passive
    /// @param prevTip The index of the previous state
    ///
    /// @return The tip of added state
    size_t addNonSourcelinkState(const TrackStatePropMask& stateMask,
                                 const BoundState& boundState,
                                 result_type& result, bool isSensitive,
                                 size_t prevTip) const {
      // Add a track state
      auto currentTip = result.fittedStates->addTrackState(stateMask, prevTip);
      if (isSensitive) {
        ACTS_VERBOSE("Creating Hole track state with tip = " << currentTip);
      } else {
        ACTS_VERBOSE("Creating Material track state with tip = " << currentTip);
      }
      // now get track state proxy back
      auto trackStateProxy = result.fittedStates->getTrackState(currentTip);

      const auto& [boundParams, jacobian, pathLength] = boundState;
      // Fill the track state
      trackStateProxy.predicted() = boundParams.parameters();
      if (boundParams.covariance().has_value()) {
        trackStateProxy.predictedCovariance() = *boundParams.covariance();
      }
      trackStateProxy.jacobian() = jacobian;
      trackStateProxy.pathLength() = pathLength;
      // Set the surface
      trackStateProxy.setReferenceSurface(
          boundParams.referenceSurface().getSharedPtr());
      // Set the filtered parameter index to be the same with predicted
      // parameter

      // Set the track state flags
      auto& typeFlags = trackStateProxy.typeFlags();
      if (trackStateProxy.referenceSurface().surfaceMaterial() != nullptr) {
        typeFlags.set(TrackStateFlag::MaterialFlag);
      }
      typeFlags.set(TrackStateFlag::ParameterFlag);
      if (isSensitive) {
        typeFlags.set(TrackStateFlag::HoleFlag);
      }

      trackStateProxy.shareFrom(TrackStatePropMask::Predicted,
                                TrackStatePropMask::Filtered);

      return currentTip;
    }

    /// @brief CombinatorialKalmanFilter actor operation : material interaction
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param surface The surface where the material interaction happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param updateStage The materal update stage
    ///
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void materialInteractor(const Surface* surface, propagator_state_t& state,
                            const stepper_t& stepper,
                            const navigator_t& navigator,
                            const MaterialUpdateStage& updateStage) const {
      // Indicator if having material
      bool hasMaterial = false;

      if (surface and surface->surfaceMaterial()) {
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
                       << interaction.Eloss << ", "
                       << "variancePhi = " << interaction.variancePhi << ", "
                       << "varianceTheta = " << interaction.varianceTheta
                       << ", "
                       << "varianceQoverP = " << interaction.varianceQoverP);

          // Update the state and stepper with material effects
          interaction.updateState(state, stepper, addNoise);
        }
      }

      if (not hasMaterial) {
        // Screen out message
        ACTS_VERBOSE("No material effects on surface: " << surface->geometryId()
                                                        << " at update stage: "
                                                        << updateStage);
      }
    }

    /// @brief Kalman actor operation : finalize
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    Result<void> finalize(propagator_state_t& state, const stepper_t& stepper,
                          result_type& result) const {
      // The measurement tip of the track being smoothed
      const auto& lastMeasurementIndex =
          result.lastMeasurementIndices.at(result.iSmoothed);

      // Get the indices of the first states (can be either a measurement or
      // material);
      size_t firstStateIndex = lastMeasurementIndex;
      // Count track states to be smoothed
      size_t nStates = 0;
      result.fittedStates->applyBackwards(lastMeasurementIndex, [&](auto st) {
        bool isMeasurement =
            st.typeFlags().test(TrackStateFlag::MeasurementFlag);
        bool isMaterial = st.typeFlags().test(TrackStateFlag::MaterialFlag);
        if (isMeasurement || isMaterial) {
          firstStateIndex = st.index();
        }
        nStates++;
      });
      // Return error if the track has no measurement states (but this should
      // not happen)
      if (nStates == 0) {
        ACTS_ERROR("Smoothing for a track without measurements.");
        return CombinatorialKalmanFilterError::SmoothFailed;
      }
      // Screen output for debugging
      ACTS_VERBOSE("Apply smoothing on " << nStates
                                         << " filtered track states.");
      // Smooth the track states
      auto smoothRes =
          m_extensions.smoother(state.geoContext, *result.fittedStates,
                                lastMeasurementIndex, getDummyLogger());
      if (!smoothRes.ok()) {
        ACTS_ERROR("Smoothing step failed: " << smoothRes.error());
        return smoothRes.error();
      }

      // Return in case no target surface
      if (targetSurface == nullptr) {
        return Result<void>::success();
      }

      // Obtain the smoothed parameters at first/last measurement state.
      // The first state can also be a material state
      auto firstCreatedState =
          result.fittedStates->getTrackState(firstStateIndex);
      auto lastCreatedMeasurement =
          result.fittedStates->getTrackState(lastMeasurementIndex);

      // Lambda to get the intersection of the free params on the target surface
      auto target = [&](const FreeVector& freeVector) -> SurfaceIntersection {
        return targetSurface->intersect(
            state.geoContext, freeVector.segment<3>(eFreePos0),
            state.stepping.navDir * freeVector.segment<3>(eFreeDir0), true);
      };

      // The smoothed free params at the first/last measurement state
      auto firstParams = MultiTrajectoryHelpers::freeSmoothed(
          state.options.geoContext, firstCreatedState);
      auto lastParams = MultiTrajectoryHelpers::freeSmoothed(
          state.options.geoContext, lastCreatedMeasurement);
      // Get the intersections of the smoothed free parameters with the target
      // surface
      const auto firstIntersection = target(firstParams);
      const auto lastIntersection = target(lastParams);

      // Update the stepping parameters - in order to progress to destination.
      // At the same time, reverse navigation direction for further
      // stepping if necessary.
      // @note The stepping parameters is updated to the smoothed parameters at
      // either the first measurement state or the last measurement state. It
      // assumes the target surface is not within the first and the last
      // smoothed measurement state. Also, whether the intersection is on
      // surface is not checked here.
      bool reverseDirection = false;
      bool closerToFirstCreatedState =
          (std::abs(firstIntersection.intersection.pathLength) <=
           std::abs(lastIntersection.intersection.pathLength));
      if (closerToFirstCreatedState) {
        stepper.update(state.stepping, firstParams,
                       firstCreatedState.smoothed(),
                       firstCreatedState.smoothedCovariance(),
                       firstCreatedState.referenceSurface());
        reverseDirection = (firstIntersection.intersection.pathLength < 0);
      } else {
        stepper.update(state.stepping, lastParams,
                       lastCreatedMeasurement.smoothed(),
                       lastCreatedMeasurement.smoothedCovariance(),
                       lastCreatedMeasurement.referenceSurface());
        reverseDirection = (lastIntersection.intersection.pathLength < 0);
      }
      const auto& surface = closerToFirstCreatedState
                                ? firstCreatedState.referenceSurface()
                                : lastCreatedMeasurement.referenceSurface();
      ACTS_VERBOSE(
          "Smoothing successful, updating stepping state to smoothed "
          "parameters at surface "
          << surface.geometryId() << ". Prepared to reach the target surface.");

      // Reverse the navigation direction if necessary
      if (reverseDirection) {
        ACTS_VERBOSE(
            "Reverse navigation direction after smoothing for reaching the "
            "target surface");
        state.stepping.navDir =
            (state.stepping.navDir == NavigationDirection::Forward)
                ? NavigationDirection::Backward
                : NavigationDirection::Forward;
      }
      // Reinitialize the stepping jacobian
      state.stepping.jacobian = BoundMatrix::Identity();
      state.stepping.jacTransport = FreeMatrix::Identity();
      state.stepping.derivative = FreeVector::Zero();
      // Reset the step size
      state.stepping.stepSize = ConstrainedStep(
          state.stepping.navDir * std::abs(state.options.maxStepSize));
      // Set accumulatd path to zero before targeting surface
      state.stepping.pathAccumulated = 0.;

      return Result<void>::success();
    }

    CombinatorialKalmanFilterExtensions<traj_t> m_extensions;

    /// The source link accesor
    source_link_accessor_t m_sourcelinkAccessor;

    /// The Surface being targeted
    SurfaceReached targetReached;

    /// Logger instance
    const Logger* actorLogger{nullptr};

    const Logger& logger() const { return *actorLogger; }
  };

  template <typename source_link_accessor_t, typename parameters_t>
  class Aborter {
   public:
    /// Broadcast the result_type
    using action_type = Actor<source_link_accessor_t, parameters_t>;

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t, typename result_t>
    bool operator()(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, const result_t& result,
                    const Logger& /*logger*/) const {
      if (!result.result.ok() or result.finished) {
        return true;
      }
      return false;
    }
  };

 public:
  /// Combinatorial Kalman Filter implementation, calls the Kalman filter
  /// and smoother
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
    combKalmanActor.targetSurface = tfOptions.referenceSurface;
    combKalmanActor.multipleScattering = tfOptions.multipleScattering;
    combKalmanActor.energyLoss = tfOptions.energyLoss;
    combKalmanActor.smoothing = tfOptions.smoothing;
    combKalmanActor.actorLogger = m_actorLogger.get();

    // copy source link accessor, calibrator and measurement selector
    combKalmanActor.m_sourcelinkAccessor = tfOptions.sourcelinkAccessor;
    combKalmanActor.m_extensions = tfOptions.extensions;

    // Run the CombinatorialKalmanFilter.
    auto stateBuffer = std::make_shared<traj_t>();

    typename propagator_t::template action_list_t_result_t<
        CurvilinearTrackParameters, Actors>
        inputResult;

    auto& r =
        inputResult.template get<CombinatorialKalmanFilterResult<traj_t>>();

    r.fittedStates = &trackContainer.trackStateContainer();
    r.stateBuffer = stateBuffer;
    r.stateBuffer->clear();

    auto result = m_propagator.template propagate(
        initialParameters, propOptions, std::move(inputResult));

    if (!result.ok()) {
      ACTS_ERROR("Propapation failed: " << result.error() << " "
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
    if (combKalmanResult.result.ok() and not combKalmanResult.finished) {
      combKalmanResult.result = Result<void>(
          CombinatorialKalmanFilterError::PropagationReachesMaxSteps);
    }

    if (!combKalmanResult.result.ok()) {
      ACTS_ERROR("CombinatorialKalmanFilter failed: "
                 << combKalmanResult.result.error() << " "
                 << combKalmanResult.result.error().message()
                 << " with the initial parameters: \n"
                 << initialParameters.parameters());
      return combKalmanResult.result.error();
    }

    std::vector<typename TrackContainer::TrackProxy> tracks;

    for (auto tip : combKalmanResult.lastMeasurementIndices) {
      auto track = trackContainer.getTrack(trackContainer.addTrack());
      track.tipIndex() = tip;
      const BoundTrackParameters& parameters =
          combKalmanResult.fittedParameters.find(tip)->second;
      track.parameters() = parameters.parameters();
      track.covariance() = *parameters.covariance();
      track.setReferenceSurface(parameters.referenceSurface().getSharedPtr());
      tracks.push_back(track);
    }

    return tracks;
  }
};

}  // namespace Acts
