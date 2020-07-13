// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/Fitter/detail/VoidKalmanComponents.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFinder/CombinatorialKalmanFilterError.hpp"
#include "Acts/TrackFinder/detail/VoidTrackFinderComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <functional>
#include <map>
#include <memory>
#include <unordered_map>

namespace Acts {

/// @brief struct to keep record of the track quality
///
/// This could be used to decide if a track is to be recorded when the forward
/// filtering is done or to be terminated due to its bad quality
/// @Todo: add other useful info, e.g. chi2
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

/// @brief Options struct how the CombinatorialKalmanFilter (CKF) is called
///
/// @tparam source_link_selector_t The source link selector type
///
/// It contains the context of the CKF call, the source link selector
/// config, the optional surface where to express the track finding/fitting
/// result, config for material effects and whether to run smoothing to get
/// fitted parameters
///
/// @note the context objects must be provided
template <typename source_link_selector_t>
struct CombinatorialKalmanFilterOptions {
  // Broadcast the source link selector type
  using SourceLinkSelector = source_link_selector_t;

  // Broadcast the source link selector config type
  using SourceLinkSelectorConfig = typename SourceLinkSelector::Config;

  /// Deleted default constructor
  CombinatorialKalmanFilterOptions() = delete;

  /// PropagatorOptions with context
  ///
  /// @param gctx The goemetry context for this track finding/fitting
  /// @param mctx The magnetic context for this track finding/fitting
  /// @param cctx The calibration context for this track finding/fitting
  /// @param slsCfg The config for the source link selector for this track
  /// finding/fitting
  /// @param rSurface The reference surface for the eventual track fitting to be
  /// expressed at
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  /// @param rSmoothing Whether to run smoothing to get fitted parameter
  CombinatorialKalmanFilterOptions(
      std::reference_wrapper<const GeometryContext> gctx,
      std::reference_wrapper<const MagneticFieldContext> mctx,
      std::reference_wrapper<const CalibrationContext> cctx,
      const SourceLinkSelectorConfig& slsCfg, const Surface* rSurface = nullptr,
      bool mScattering = true, bool eLoss = true, bool rSmoothing = true)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        sourcelinkSelectorConfig(slsCfg),
        referenceSurface(rSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        smoothing(rSmoothing) {}

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  /// The config for the source link selector
  SourceLinkSelectorConfig sourcelinkSelectorConfig;

  /// The reference Surface
  const Surface* referenceSurface = nullptr;

  /// Whether to consider multiple scattering.
  bool multipleScattering = true;

  /// Whether to consider energy loss.
  bool energyLoss = true;

  /// Whether to run smoothing to get fitted parameter
  bool smoothing = true;
};

template <typename source_link_t>
struct CombinatorialKalmanFilterResult {
  // Fitted states that the actor has handled.
  MultiTrajectory<source_link_t> fittedStates;

  // The indices of the 'tip' of the tracks stored in multitrajectory.
  std::vector<size_t> trackTips;

  // The Parameters at the provided surface for separate tracks
  std::unordered_map<size_t, BoundParameters> fittedParameters;

  // The indices of the 'tip' of the unfinished tracks
  std::vector<std::pair<size_t, CombinatorialKalmanFilterTipState>> activeTips;

  // The indices of source links in multitrajectory
  std::unordered_map<const Surface*, std::unordered_map<size_t, size_t>>
      sourcelinkTips;

  // Indicator if forward filtering has been done
  bool forwardFiltered = false;

  // Indicator if smoothing has been done.
  bool smoothed = false;

  // The index for the current smoothing track
  size_t iSmoothed = 0;

  // Indicator if the propagation state has been reset
  bool reset = false;

  // Indicator if track finding has been done
  bool finished = false;

  // Temporary container for index and chi2 of intermediate source link
  // candidates
  std::vector<std::pair<size_t, double>> sourcelinkChi2;

  // Temporary container for index of final source link candidates
  std::vector<size_t> sourcelinkCandidateIndices;

  Result<void> result{Result<void>::success()};
};

/// @brief CombinatorialKalmanFilter implementation of Acts as a plugin
///
/// to the Propgator
///
/// @tparam propagator_t Type of the propagation class
/// @tparam updater_t Type of the kalman updater class
/// @tparam smoother_t Type of the kalman smoother class
/// @tparam source_link_selector_t Type of the source link selector class
/// @tparam branch_stopper_t Type of the branch stopper class
/// @tparam calibrator_t Type of the calibrator class
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
/// - The Smoother is called at the end of the forward track finding by the
/// Actor.
/// - The Sourcelink selector is called during the filtering by the Actor.
/// - The Calibrator is a dedicated calibration algorithm that allows
///   to calibrate measurements using track information, this could be
///    e.g. sagging for wires, module deformations, etc.
///
/// Measurements are not required to be ordered for the
/// CombinatorialKalmanFilter, measurement ordering needs to be figured out by
/// the navigation of the propagator.
///
/// The void components are provided mainly for unit testing.
template <typename propagator_t, typename updater_t = VoidKalmanUpdater,
          typename smoother_t = VoidKalmanSmoother,
          typename source_link_selector_t = CKFSourceLinkSelector,
          typename branch_stopper_t = VoidBranchStopper,
          typename calibrator_t = VoidMeasurementCalibrator>
class CombinatorialKalmanFilter {
 public:
  /// Shorthand definition
  using MeasurementSurfaces = std::multimap<const Layer*, const Surface*>;
  /// The navigator type
  using KalmanNavigator = typename propagator_t::Navigator;

  /// Default constructor is deleted
  CombinatorialKalmanFilter() = delete;

  /// Constructor from arguments
  CombinatorialKalmanFilter(propagator_t pPropagator,
                            std::unique_ptr<const Logger> logger =
                                getDefaultLogger("CombinatorialKalmanFilter",
                                                 Logging::INFO))
      : m_propagator(std::move(pPropagator)), m_logger(logger.release()) {}

 private:
  /// The propgator for the transport and material update
  propagator_t m_propagator;

  /// Logger getter to support macros
  const Logger& logger() const { return *m_logger; }

  /// Owned logging instance
  std::shared_ptr<const Logger> m_logger;

  /// @brief Propagator Actor plugin for the CombinatorialKalmanFilter
  ///
  /// @tparam source_link_t is an type fulfilling the @c SourceLinkConcept
  /// @tparam parameters_t The type of parameters used for "local" paremeters.
  ///
  /// The CombinatorialKalmanFilterActor does not rely on the measurements to be
  /// sorted along the track.
  template <typename source_link_t, typename parameters_t>
  class Actor {
   public:
    using TipState = CombinatorialKalmanFilterTipState;
    using BoundState = std::tuple<BoundParameters, BoundMatrix, double>;
    using CurvilinearState =
        std::tuple<CurvilinearParameters, BoundMatrix, double>;
    /// Broadcast the result_type
    using result_type = CombinatorialKalmanFilterResult<source_link_t>;

    /// The target surface
    const Surface* targetSurface = nullptr;

    /// Allows retrieving measurements for a surface
    std::unordered_map<const Surface*, std::vector<source_link_t>>
        inputMeasurements;

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// Whether to consider energy loss.
    bool energyLoss = true;

    /// Whether to run smoothing to get fitted parameter
    bool smoothing = true;

    /// @brief CombinatorialKalmanFilter actor operation
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    void operator()(propagator_state_t& state, const stepper_t& stepper,
                    result_type& result) const {
      ACTS_VERBOSE("CombinatorialKalmanFilter step");

      // This following is added due to the fact that the navigation
      // reinitialization in reset call cannot guarantee the navigator to target
      // for extra layers in the reset volume.
      // Currently, manually set navigation stage to allow for targeting layers
      // after all the surfaces on the reset layer has been processed.
      // Otherwise, the navigation stage will be Stage::boundaryTarget after
      // navigator status call which means the extra layers on the reset volume
      // won't be targeted.
      // @Todo: Let the navigator do all the re-initialization
      if (result.reset and state.navigation.navSurfaceIter ==
                               state.navigation.navSurfaces.end()) {
        // So the navigator target call will target layers
        state.navigation.navigationStage = KalmanNavigator::Stage::layerTarget;
        // We only do this after the reset layer has been processed
        result.reset = false;
      }

      // Update:
      // - Waiting for a current surface
      auto surface = state.navigation.currentSurface;
      if (surface != nullptr and not result.forwardFiltered) {
        // There are three scenarios:
        // 1) The surface is in the measurement map
        // -> Select source links
        // -> Perform the kalman update for selected non-outlier source links
        // -> Add track states in multitrajectory. Multiple states mean branch
        // splitting.
        // -> Call branch stopper to justify each branch
        // -> If there is non-outlier state, update stepper information
        // 2) The surface is not in the measurement map but with material
        // -> Add a hole or passive material state in multitrajectory
        // -> Call branch stopper to justify the branch
        // 3) The surface is neither in the measurement map nor with material
        // -> Do nothing
        ACTS_VERBOSE("Perform filter step");
        auto res = filter(surface, state, stepper, result);
        if (!res.ok()) {
          ACTS_ERROR("Error in filter: " << res.error());
          result.result = res.error();
        }
      }

      // Reset propagation state:
      // - When navigation breaks and there is stil active tip present after
      // recording&removing track tips on current surface
      if (state.navigation.navigationBreak and not result.forwardFiltered) {
        // Record the tips on current surface as trajectory entry indices
        // (taking advantage of fact that those tips are consecutive in list of
        // active tips) and remove those tips from active tips
        if (not result.activeTips.empty()) {
          // The last active tip
          const auto& lastActiveTip = result.activeTips.back().first;
          // Get the index of previous state
          const auto& iprevious =
              result.fittedStates.getTrackState(lastActiveTip).previous();
          // Find the track states which have the same previous state and remove
          // them from active tips
          while (not result.activeTips.empty()) {
            const auto& [currentTip, tipState] = result.activeTips.back();
            if (result.fittedStates.getTrackState(currentTip).previous() !=
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
              result.trackTips.emplace_back(currentTip);
            }
            // Remove the tip from list of active tips
            result.activeTips.erase(result.activeTips.end() - 1);
          }
        }
        // If no more active tip, done with forward filtering; Otherwise, reset
        // propagation state to track state at last tip of active tips
        if (result.activeTips.empty()) {
          ACTS_VERBOSE("Forward Kalman filtering finds "
                       << result.trackTips.size() << " tracks");
          result.forwardFiltered = true;
        } else {
          ACTS_VERBOSE("Propagation jumps to branch with tip = "
                       << result.activeTips.back().first);
          reset(state, stepper, result);
        }
      }

      // Post-processing after forward filtering
      if (result.forwardFiltered and not result.finished) {
        // Return error if forward filtering finds no tracks
        if (result.trackTips.empty()) {
          result.result =
              Result<void>(CombinatorialKalmanFilterError::NoTracksFound);
        } else {
          if (not smoothing) {
            // Manually set the targetReached to abort the propagation
            ACTS_VERBOSE("Finish forward Kalman filtering");
            // Remember that track finding is done
            result.finished = true;
            state.navigation.targetReached = true;
          } else {
            // Iterate over the found tracks for smoothing and getting the
            // fitted parameter. This needs to be accomplished in different
            // propagation steps:
            // -> first run smoothing for found track indexed with iSmoothed
            if (not result.smoothed) {
              ACTS_VERBOSE(
                  "Finalize/run smoothing for track with entry index = "
                  << result.trackTips.at(result.iSmoothed));
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
            if (result.smoothed and
                targetReached(state, stepper, *targetSurface)) {
              ACTS_VERBOSE("Completing the track with entry index = "
                           << result.trackTips.at(result.iSmoothed));
              // Transport & bind the parameter to the final surface
              auto fittedState =
                  stepper.boundState(state.stepping, *targetSurface);
              // Assign the fitted parameters
              result.fittedParameters.emplace(
                  result.trackTips.at(result.iSmoothed),
                  std::get<BoundParameters>(fittedState));
              // If there are more trajectories to handle:
              // -> set the targetReached status to false
              // -> set the smoothed status to false
              // -> update the index of track to be smoothed
              if (result.iSmoothed < result.trackTips.size() - 1) {
                state.navigation.targetReached = false;
                result.smoothed = false;
                result.iSmoothed++;
                // To avoid meaningless navigation target call
                state.stepping.stepSize =
                    ConstrainedStep(state.options.maxStepSize);
                // Need to go back to start targeting for the rest tracks
                state.stepping.navDir = forward;
              } else {
                ACTS_VERBOSE(
                    "Finish forward Kalman filtering and backward smoothing");
                // Remember that track finding is done
                result.finished = true;
              }
            }
          }  // if run smoothing
        }    // if there are found tracks
      }      // if forward filtering is done
    }

    /// @brief Kalman actor operation : reset propagation
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    void reset(propagator_state_t& state, stepper_t& stepper,
               result_type& result) const {
      // Remember the propagation state has been reset
      result.reset = true;

      auto currentState =
          result.fittedStates.getTrackState(result.activeTips.back().first);
      // Reset the navigation state
      state.navigation = typename propagator_t::NavigatorState();
      state.navigation.startSurface = &currentState.referenceSurface();
      if (state.navigation.startSurface->associatedLayer() != nullptr) {
        state.navigation.startLayer =
            state.navigation.startSurface->associatedLayer();
      }
      state.navigation.startVolume =
          state.navigation.startLayer->trackingVolume();
      state.navigation.targetSurface = targetSurface;
      state.navigation.currentSurface = state.navigation.startSurface;
      state.navigation.currentVolume = state.navigation.startVolume;

      // Update the stepping state
      stepper.update(state.stepping,
                     MultiTrajectoryHelpers::freeFiltered(
                         state.options.geoContext, currentState),
                     currentState.filteredCovariance());
      // Reinitialize the stepping jacobian
      currentState.referenceSurface().initJacobianToGlobal(
          state.options.geoContext, state.stepping.jacToGlobal,
          state.stepping.pos, state.stepping.dir, currentState.filtered());
      state.stepping.jacobian = BoundMatrix::Identity();
      state.stepping.jacTransport = FreeMatrix::Identity();
      state.stepping.derivative = FreeVector::Zero();
      // Reset step size and accumulated path
      state.stepping.stepSize = ConstrainedStep(state.options.maxStepSize);
      state.stepping.pathAccumulated = currentState.pathLength();

      // No Kalman filtering for the starting surface, but still need
      // to consider the material effects here
      materialInteractor(state.navigation.startSurface, state, stepper);
    }

    /// @brief CombinatorialKalmanFilter actor operation :
    /// - filtering for all measurement(s) on surface
    /// - store selected track states in multiTrajectory
    /// - update propagator state to the (last) selected track state
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result The mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    Result<void> filter(const Surface* surface, propagator_state_t& state,
                        const stepper_t& stepper, result_type& result) const {
      // Initialize the number of branches on current surface
      size_t nBranchesOnSurface = 0;

      // Try to find the surface in the measurement surfaces
      auto sourcelink_it = inputMeasurements.find(surface);
      if (sourcelink_it != inputMeasurements.end()) {
        // Screen output message
        ACTS_VERBOSE("Measurement surface " << surface->geoID()
                                            << " detected.");

        // Update state and stepper with pre material effects
        materialInteractor(surface, state, stepper, preUpdate);

        // Transport & bind the state to the current surface
        auto boundState = stepper.boundState(state.stepping, *surface);
        auto boundParams = std::get<BoundParameters>(boundState);

        // Get all source links on surface
        auto& sourcelinks = sourcelink_it->second;

        // Invoke the source link selector to select source links for either
        // measurements or outlier.
        // Calibrator is passed to the selector because
        // selection has to be done based on calibrated measurement
        bool isOutlier = false;
        auto sourcelinkSelectionRes = m_sourcelinkSelector(
            m_calibrator, boundParams, sourcelinks, result.sourcelinkChi2,
            result.sourcelinkCandidateIndices, isOutlier);
        if (!sourcelinkSelectionRes.ok()) {
          ACTS_ERROR("Selection of source links failed: "
                     << sourcelinkSelectionRes.error());
          return sourcelinkSelectionRes.error();
        }

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

        // Remember the tip of the neighbor state on this surface
        size_t neighborTip = SIZE_MAX;
        // Loop over the selected source links
        for (const auto& index : result.sourcelinkCandidateIndices) {
          // Determine if predicted parameter is already contained in
          // neighboring state
          bool isPredictedShared = (neighborTip != SIZE_MAX);

          // Determine if uncalibrated measurement are already
          // contained in other track state
          bool isSourcelinkShared = false;
          size_t sharedTip = SIZE_MAX;
          auto sourcelinkTips_it = result.sourcelinkTips.find(surface);
          if (sourcelinkTips_it != result.sourcelinkTips.end()) {
            auto& sourcelinkTipsOnSurface = sourcelinkTips_it->second;
            auto index_it = sourcelinkTipsOnSurface.find(index);
            if (index_it != sourcelinkTipsOnSurface.end()) {
              isSourcelinkShared = true;
              sharedTip = index_it->second;
            }
          }

          // The mask for adding a state in the multitrajectory
          // No storage allocation for:
          // -> predicted parameter and uncalibrated measurement if
          // already stored
          // -> filtered parameter for outlier
          auto stateMask =
              (isPredictedShared ? ~TrackStatePropMask::Predicted
                                 : TrackStatePropMask::All) &
              (isSourcelinkShared ? ~TrackStatePropMask::Uncalibrated
                                  : TrackStatePropMask::All) &
              (isOutlier ? ~TrackStatePropMask::Filtered
                         : TrackStatePropMask::All);

          // Add measurement/outlier track state to the multitrajectory
          auto addStateRes = addSourcelinkState(
              stateMask, boundState, sourcelinks.at(index), isOutlier, result,
              state.geoContext, prevTip, prevTipState, neighborTip, sharedTip);
          if (addStateRes.ok()) {
            const auto& [currentTip, tipState] = addStateRes.value();
            // Remember the track state tip for this stored source link
            if (not isSourcelinkShared) {
              auto& sourcelinkTipsOnSurface = result.sourcelinkTips[surface];
              sourcelinkTipsOnSurface.emplace(index, currentTip);
            }
            // Remember the tip of neighbor state on this surface
            neighborTip = currentTip;

            // Check if need to stop this branch
            if (not m_branchStopper(tipState)) {
              // Remember the active tip and its state
              result.activeTips.emplace_back(std::move(currentTip),
                                             std::move(tipState));
              // Record the number of branches on surface
              nBranchesOnSurface++;
            }
          }
        }  // end of loop for all selected source links on this surface

        if (nBranchesOnSurface > 0 and not isOutlier) {
          // If there are measurement track states on this surface
          ACTS_VERBOSE("Filtering step successful with " << nBranchesOnSurface
                                                         << " branches");
          // Update stepping state using filtered parameters of last track
          // state on this surface
          auto ts =
              result.fittedStates.getTrackState(result.activeTips.back().first);
          stepper.update(state.stepping,
                         MultiTrajectoryHelpers::freeFiltered(
                             state.options.geoContext, ts),
                         ts.filteredCovariance());
          ACTS_VERBOSE("Stepping state is updated with filtered parameter: \n"
                       << ts.filtered().transpose()
                       << " of track state with tip = "
                       << result.activeTips.back().first);
        }
        // Update state and stepper with post material effects
        materialInteractor(surface, state, stepper, postUpdate);
      } else if (surface->surfaceMaterial() != nullptr) {
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
        std::string type = isSensitive ? "sensitive" : "passive";
        ACTS_VERBOSE("Detected " << type << " surface: " << surface->geoID());
        if (isSensitive) {
          // Increment of number of passed sensitive surfaces
          tipState.nSensitiveSurfaces++;
        }
        // Add state if there is already measurement detected on this branch
        // For in-sensitive surface, only add state when smoothing is
        // required
        if (tipState.nMeasurements > 0 and
            (isSensitive or (not isSensitive and smoothing))) {
          // New state is to be added. Remove the last tip from active tips now
          result.activeTips.erase(result.activeTips.end() - 1);

          // No source links on surface, add either hole or passive material
          // TrackState. No storage allocation for uncalibrated/calibrated
          // measurement and filtered parameter
          auto stateMask =
              ~(TrackStatePropMask::Uncalibrated |
                TrackStatePropMask::Calibrated | TrackStatePropMask::Filtered);

          // Increment of number of processed states
          tipState.nStates++;
          size_t currentTip = SIZE_MAX;
          if (isSensitive) {
            // Transport & bind the state to the current surface
            auto boundState = stepper.boundState(state.stepping, *surface);
            // Add a hole track state to the multitrajectory
            currentTip = addHoleState(stateMask, boundState, result, prevTip);
            // Incremet of number of holes
            tipState.nHoles++;
          } else {
            // Transport & get curvilinear state instead of bound state
            auto curvilinearState = stepper.curvilinearState(state.stepping);
            // Add a passive material track state to the multitrajectory
            currentTip =
                addPassiveState(stateMask, curvilinearState, result, prevTip);
          }

          // Check the branch
          if (not m_branchStopper(tipState)) {
            // Remember the active tip and its state
            result.activeTips.emplace_back(std::move(currentTip),
                                           std::move(tipState));
          } else {
            // No branch on this surface
            nBranchesOnSurface = 0;
          }
        }
        // Update state and stepper with material effects
        materialInteractor(surface, state, stepper, fullUpdate);
      } else {
        // Neither measurement nor material on surface, this branch is still
        // valid. Count the branch on current surface
        nBranchesOnSurface = 1;
      }

      // Reset current tip if there is no branch on current surface
      if (nBranchesOnSurface == 0) {
        ACTS_DEBUG("Branch on surface " << surface->geoID() << " is stopped");
        if (not result.activeTips.empty()) {
          ACTS_VERBOSE("Propagation jumps to branch with tip = "
                       << result.activeTips.back().first);
          reset(state, stepper, result);
        } else {
          ACTS_VERBOSE("Stop forward Kalman filtering with "
                       << result.trackTips.size() << " found tracks");
          result.forwardFiltered = true;
        }
      }

      return Result<void>::success();
    }

    /// @brief CombinatorialKalmanFilter actor operation : add track state with
    /// source link: measurement or outlier
    ///
    /// @param stateMask The bitmask that instructs which components to allocate
    /// and which to leave invalid
    /// @param boundState The bound state on current surface
    /// @param sourcelink The source link to be stored
    /// @param isOutlier Indicator for outlier or not
    /// @param result is the mutable result state object
    /// @param geoContext The geometry context (needed for Kalman update)
    /// @param neighborTip The neighbor state tip on this surface (the predicted
    /// parameters could be shared between neighbors)
    /// @param sharedTip The tip of state with shared source link
    ///
    /// @return The tip of added state and its state
    Result<std::pair<size_t, TipState>> addSourcelinkState(
        const TrackStatePropMask& stateMask, const BoundState& boundState,
        const source_link_t& sourcelink, bool isOutlier, result_type& result,
        std::reference_wrapper<const GeometryContext> geoContext,
        const size_t& prevTip, const TipState& prevTipState,
        size_t neighborTip = SIZE_MAX, size_t sharedTip = SIZE_MAX) const {
      // Inherit the tip state from the previous and will be updated later
      TipState tipState = prevTipState;

      // Add a track state
      auto currentTip = result.fittedStates.addTrackState(stateMask, prevTip);

      // Get the track state proxy
      auto trackStateProxy = result.fittedStates.getTrackState(currentTip);

      auto [boundParams, jacobian, pathLength] = boundState;

      // Fill the parametric part of the track state proxy
      if ((not ACTS_CHECK_BIT(stateMask, TrackStatePropMask::Predicted)) and
          neighborTip != SIZE_MAX) {
        // The predicted parameter is already stored, just set the index
        auto neighborState = result.fittedStates.getTrackState(neighborTip);
        trackStateProxy.data().ipredicted = neighborState.data().ipredicted;
      } else {
        trackStateProxy.predicted() = boundParams.parameters();
        trackStateProxy.predictedCovariance() = *boundParams.covariance();
      }
      trackStateProxy.jacobian() = jacobian;
      trackStateProxy.pathLength() = pathLength;

      // Assign the uncalibrated&calibrated measurement to the track
      // state (the uncalibrated could be already stored in other states)
      if ((not ACTS_CHECK_BIT(stateMask, TrackStatePropMask::Uncalibrated)) and
          sharedTip != SIZE_MAX) {
        // The uncalibrated are already stored, just set the
        // index
        auto shared = result.fittedStates.getTrackState(sharedTip);
        trackStateProxy.data().iuncalibrated = shared.data().iuncalibrated;
      } else {
        trackStateProxy.uncalibrated() = sourcelink;
      }
      std::visit(
          [&](const auto& calibrated) {
            trackStateProxy.setCalibrated(calibrated);
          },
          m_calibrator(trackStateProxy.uncalibrated(),
                       trackStateProxy.predicted()));

      // Get and set the type flags
      auto& typeFlags = trackStateProxy.typeFlags();
      typeFlags.set(TrackStateFlag::MaterialFlag);
      typeFlags.set(TrackStateFlag::ParameterFlag);

      // Increment of number of processedState and passed sensitive surfaces
      tipState.nSensitiveSurfaces++;
      tipState.nStates++;

      if (isOutlier) {
        ACTS_VERBOSE("Creating outlier track state with tip = " << currentTip);
        // Set the outlier flag
        typeFlags.set(TrackStateFlag::OutlierFlag);
        // Increment number of outliers
        tipState.nOutliers++;
        // No Kalman update for outlier
        // Set the filtered parameter index to be the same with predicted
        // parameter
        trackStateProxy.data().ifiltered = trackStateProxy.data().ipredicted;
      } else {
        // Kalman update
        auto updateRes = m_updater(geoContext, trackStateProxy);
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
      return std::make_pair(std::move(currentTip), std::move(tipState));
    }

    /// @brief CombinatorialKalmanFilter actor operation : add hole track state
    ///
    /// @param stateMask The bitmask that instructs which components to allocate
    /// @param boundState The bound state on current surface
    /// @param result is the mutable result state object
    /// and which to leave invalid
    /// @param prevTip The index of the previous state
    ///
    /// @return The tip of added state
    size_t addHoleState(const TrackStatePropMask& stateMask,
                        const BoundState& boundState, result_type& result,
                        size_t prevTip = SIZE_MAX) const {
      // Add a track state
      auto currentTip = result.fittedStates.addTrackState(stateMask, prevTip);
      ACTS_VERBOSE("Creating Hole track state with tip = " << currentTip);

      // now get track state proxy back
      auto trackStateProxy = result.fittedStates.getTrackState(currentTip);

      // Set the track state flags
      auto& typeFlags = trackStateProxy.typeFlags();
      typeFlags.set(TrackStateFlag::MaterialFlag);
      typeFlags.set(TrackStateFlag::ParameterFlag);
      typeFlags.set(TrackStateFlag::HoleFlag);

      auto [boundParams, jacobian, pathLength] = boundState;
      // Fill the track state
      trackStateProxy.predicted() = boundParams.parameters();
      trackStateProxy.predictedCovariance() = *boundParams.covariance();
      trackStateProxy.jacobian() = jacobian;
      trackStateProxy.pathLength() = pathLength;
      // Set the surface
      trackStateProxy.setReferenceSurface(
          boundParams.referenceSurface().getSharedPtr());
      // Set the filtered parameter index to be the same with predicted
      // parameter
      trackStateProxy.data().ifiltered = trackStateProxy.data().ipredicted;

      return currentTip;
    }

    /// @brief CombinatorialKalmanFilter actor operation : add passive track
    /// state
    ///
    /// @param stateMask The bitmask that instructs which components to allocate
    /// @param curvilinearState The curvilinear state on in-sensive material
    /// surface
    /// @param result is the mutable result state object
    /// and which to leave invalid
    /// @param prevTip The index of the previous state
    ///
    /// @return The tip of added state
    size_t addPassiveState(const TrackStatePropMask& stateMask,
                           const CurvilinearState& curvilinearState,
                           result_type& result,
                           size_t prevTip = SIZE_MAX) const {
      // Add a track state
      auto currentTip = result.fittedStates.addTrackState(stateMask, prevTip);
      ACTS_VERBOSE(
          "Creating track state on in-sensitive material surface with tip = "
          << currentTip);

      // now get track state proxy back
      auto trackStateProxy = result.fittedStates.getTrackState(currentTip);

      // Set the track state flags
      auto& typeFlags = trackStateProxy.typeFlags();
      typeFlags.set(TrackStateFlag::MaterialFlag);
      typeFlags.set(TrackStateFlag::ParameterFlag);

      auto [curvilinearParams, jacobian, pathLength] = curvilinearState;
      // Fill the track state
      trackStateProxy.predicted() = curvilinearParams.parameters();
      trackStateProxy.predictedCovariance() = *curvilinearParams.covariance();
      trackStateProxy.jacobian() = jacobian;
      trackStateProxy.pathLength() = pathLength;
      // Set the surface
      trackStateProxy.setReferenceSurface(Surface::makeShared<PlaneSurface>(
          curvilinearParams.position(), curvilinearParams.momentum()));
      // Set the filtered parameter index to be the same with predicted
      // parameter
      trackStateProxy.data().ifiltered = trackStateProxy.data().ipredicted;

      return currentTip;
    }

    /// @brief CombinatorialKalmanFilter actor operation : material interaction
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param surface The surface where the material interaction happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param updateStage The materal update stage
    ///
    template <typename propagator_state_t, typename stepper_t>
    void materialInteractor(
        const Surface* surface, propagator_state_t& state, stepper_t& stepper,
        const MaterialUpdateStage& updateStage = fullUpdate) const {
      // Indicator if having material
      bool hasMaterial = false;

      if (surface and surface->surfaceMaterial()) {
        // Prepare relevant input particle properties
        detail::PointwiseMaterialInteraction interaction(surface, state,
                                                         stepper);
        // Evaluate the material properties
        if (interaction.evaluateMaterialProperties(state, updateStage)) {
          // Surface has material at this stage
          hasMaterial = true;

          // Evaluate the material effects
          interaction.evaluatePointwiseMaterialInteraction(multipleScattering,
                                                           energyLoss);

          // Screen out material effects info
          ACTS_VERBOSE("Material effects on surface: "
                       << surface->geoID()
                       << " at update stage: " << updateStage << " are :");
          ACTS_VERBOSE("eLoss = "
                       << interaction.Eloss << ", "
                       << "variancePhi = " << interaction.variancePhi << ", "
                       << "varianceTheta = " << interaction.varianceTheta
                       << ", "
                       << "varianceQoverP = " << interaction.varianceQoverP);

          // Update the state and stepper with material effects
          interaction.updateState(state, stepper);
        }
      }

      if (not hasMaterial) {
        // Screen out message
        ACTS_VERBOSE("No material effects on surface: " << surface->geoID()
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
      // The tip of the track being smoothed
      const auto& currentTip = result.trackTips.at(result.iSmoothed);

      // Get the indices of measurement states;
      std::vector<size_t> measurementIndices;
      // Count track states to be smoothed
      size_t nStates = 0;
      result.fittedStates.applyBackwards(currentTip, [&](auto st) {
        bool isMeasurement =
            st.typeFlags().test(TrackStateFlag::MeasurementFlag);
        if (isMeasurement) {
          measurementIndices.emplace_back(st.index());
        } else if (measurementIndices.empty()) {
          // No smoothed parameter if the last measurment state has not been
          // found yet
          st.data().ismoothed = detail_lt::IndexData::kInvalid;
        }
        // Start count when the last measurement state is found
        if (not measurementIndices.empty()) {
          nStates++;
        }
      });
      // Return error if the track has no measurement states (but this should
      // not happen)
      if (measurementIndices.empty()) {
        ACTS_ERROR("Smoothing for a track without measurements.");
        return CombinatorialKalmanFilterError::SmoothFailed;
      }
      // Screen output for debugging
      ACTS_VERBOSE("Apply smoothing on " << nStates
                                         << " filtered track states.");
      // Smooth the track states
      auto smoothRes = m_smoother(state.geoContext, result.fittedStates,
                                  measurementIndices.front());
      if (!smoothRes.ok()) {
        ACTS_ERROR("Smoothing step failed: " << smoothRes.error());
        return smoothRes.error();
      }
      // Obtain the smoothed parameters at first measurement state
      auto firstMeasurement =
          result.fittedStates.getTrackState(measurementIndices.back());

      // Update the stepping parameters - in order to progress to destination
      ACTS_VERBOSE(
          "Smoothing successful, updating stepping state, "
          "set target surface.");
      stepper.update(state.stepping,
                     MultiTrajectoryHelpers::freeSmoothed(
                         state.options.geoContext, firstMeasurement),
                     firstMeasurement.smoothedCovariance());
      // Reverse the propagation direction
      state.stepping.stepSize =
          ConstrainedStep(-1. * state.options.maxStepSize);
      state.stepping.navDir = backward;
      // Set accumulatd path to zero before targeting surface
      state.stepping.pathAccumulated = 0.;
      // Not sure if the following line helps anything
      state.options.direction = backward;

      return Result<void>::success();
    }

    /// Pointer to a logger that is owned by the parent,
    /// CombinatorialKalmanFilter
    const Logger* m_logger;

    /// Getter for the logger, to support logging macros
    const Logger& logger() const { return *m_logger; }

    /// The CombinatorialKalmanFilter updater
    updater_t m_updater;

    /// The CombinatorialKalmanFilter smoother
    smoother_t m_smoother;

    /// The source link selector
    source_link_selector_t m_sourcelinkSelector;

    /// The branch propagation stopper
    branch_stopper_t m_branchStopper;

    /// The Measuremetn calibrator
    calibrator_t m_calibrator;

    /// The Surface beeing
    SurfaceReached targetReached;
  };

  template <typename source_link_t, typename parameters_t>
  class Aborter {
   public:
    /// Broadcast the result_type
    using action_type = Actor<source_link_t, parameters_t>;

    template <typename propagator_state_t, typename stepper_t,
              typename result_t>
    bool operator()(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const result_t& result) const {
      if (!result.result.ok()) {
        return true;
      }
      return false;
    }
  };

 public:
  /// Fit implementation of the foward filter, calls the
  /// the forward filter and backward smoother
  ///
  /// @tparam source_link_container_t Source link container type
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param tfOptions CombinatorialKalmanFilterOptions steering the track
  /// finding
  /// @note The input measurements are given in the form of @c SourceLinks.
  /// It's
  /// @c calibrator_t's job to turn them into calibrated measurements used in
  /// the track finding.
  ///
  /// @return the output as an output track
  template <typename source_link_container_t, typename start_parameters_t,
            typename parameters_t = BoundParameters>
  Result<CombinatorialKalmanFilterResult<
      typename source_link_container_t::value_type>>
  findTracks(const source_link_container_t& sourcelinks,
             const start_parameters_t& sParameters,
             const CombinatorialKalmanFilterOptions<source_link_selector_t>&
                 tfOptions) const {
    using SourceLink = typename source_link_container_t::value_type;
    static_assert(SourceLinkConcept<SourceLink>,
                  "Source link does not fulfill SourceLinkConcept");

    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_VERBOSE("Preparing " << sourcelinks.size() << " input measurements");
    std::unordered_map<const Surface*, std::vector<SourceLink>>
        inputMeasurements;
    for (const auto& sl : sourcelinks) {
      const Surface* srf = &sl.referenceSurface();
      inputMeasurements[srf].emplace_back(sl);
    }

    // Create the ActionList and AbortList
    using CombinatorialKalmanFilterAborter = Aborter<SourceLink, parameters_t>;
    using CombinatorialKalmanFilterActor = Actor<SourceLink, parameters_t>;
    using CombinatorialKalmanFilterResult =
        typename CombinatorialKalmanFilterActor::result_type;
    using Actors = ActionList<CombinatorialKalmanFilterActor>;
    using Aborters = AbortList<CombinatorialKalmanFilterAborter>;

    // Create relevant options for the propagation options
    PropagatorOptions<Actors, Aborters> propOptions(tfOptions.geoContext,
                                                    tfOptions.magFieldContext);
    // Set max steps
    propOptions.maxSteps = 10000;

    // Catch the actor and set the measurements
    auto& combKalmanActor =
        propOptions.actionList.template get<CombinatorialKalmanFilterActor>();
    combKalmanActor.m_logger = m_logger.get();
    combKalmanActor.inputMeasurements = std::move(inputMeasurements);
    combKalmanActor.targetSurface = tfOptions.referenceSurface;
    combKalmanActor.multipleScattering = tfOptions.multipleScattering;
    combKalmanActor.energyLoss = tfOptions.energyLoss;
    combKalmanActor.smoothing = tfOptions.smoothing;

    // Set config and logger for source link selector
    combKalmanActor.m_sourcelinkSelector.m_config =
        tfOptions.sourcelinkSelectorConfig;
    combKalmanActor.m_sourcelinkSelector.m_logger = m_logger;

    // also set logger on updater and smoother
    combKalmanActor.m_updater.m_logger = m_logger;
    combKalmanActor.m_smoother.m_logger = m_logger;

    // Run the CombinatorialKalmanFilter
    auto result = m_propagator.template propagate(sParameters, propOptions);

    if (!result.ok()) {
      ACTS_ERROR("Propapation failed: " << result.error());
      return result.error();
    }

    const auto& propRes = *result;

    /// Get the result of the CombinatorialKalmanFilter
    auto combKalmanResult =
        propRes.template get<CombinatorialKalmanFilterResult>();

    /// The propagation could already reach max step size
    /// before the track finding is finished during two phases:
    // -> forward filtering for track finding
    // -> surface targeting to get fitted parameters at target surface
    // @TODO: Implement distinguishment between the above two cases if necessary
    if (combKalmanResult.result.ok() and not combKalmanResult.finished) {
      combKalmanResult.result = Result<void>(
          CombinatorialKalmanFilterError::PropagationReachesMaxSteps);
    }

    if (!combKalmanResult.result.ok()) {
      ACTS_ERROR("CombinatorialKalmanFilter failed: "
                 << combKalmanResult.result.error());
      return combKalmanResult.result.error();
    }

    // Return the converted Track
    return combKalmanResult;
  }

};  // namespace Acts

}  // namespace Acts