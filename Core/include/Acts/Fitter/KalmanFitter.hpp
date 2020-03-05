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
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/EventData/TrackStateSorters.hpp"
#include "Acts/Fitter/KalmanFitterError.hpp"
#include "Acts/Fitter/detail/VoidKalmanComponents.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <functional>
#include <map>
#include <memory>

namespace Acts {

/// @brief Options struct how the Fitter is called
///
/// It contains the context of the fitter call, the outlier finder, the optional
/// surface where to express the fit result and configurations for material
/// effects and smoothing options
///
///
/// @note the context objects must be provided
template <typename outlier_finder_t = VoidOutlierFinder>
struct KalmanFitterOptions {
  // Broadcast the outlier finder type
  using OutlierFinder = outlier_finder_t;

  /// Deleted default constructor
  KalmanFitterOptions() = delete;

  /// PropagatorOptions with context
  ///
  /// @param gctx The goemetry context for this fit
  /// @param mctx The magnetic context for this fit
  /// @param cctx The calibration context for this fit
  /// @param olCfg The config for the outlier finder
  /// @param rSurface The reference surface for the fit to be expressed at
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  /// @param bwdFiltering Whether to run backward filtering as smoothing
  KalmanFitterOptions(std::reference_wrapper<const GeometryContext> gctx,
                      std::reference_wrapper<const MagneticFieldContext> mctx,
                      std::reference_wrapper<const CalibrationContext> cctx,
                      const OutlierFinder& outlierFinder,
                      const Surface* rSurface = nullptr,
                      bool mScattering = true, bool eLoss = true,
                      bool bwdFiltering = false)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        outlierFinder(outlierFinder),
        referenceSurface(rSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        backwardFiltering(bwdFiltering) {}

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  /// The config for the outlier finder
  OutlierFinder outlierFinder;

  /// The reference Surface
  const Surface* referenceSurface = nullptr;

  /// Whether to consider multiple scattering
  bool multipleScattering = true;

  /// Whether to consider energy loss
  bool energyLoss = true;

  /// Whether to run backward filtering.
  bool backwardFiltering = false;
};

template <typename source_link_t>
struct KalmanFitterResult {
  // Fitted states that the actor has handled.
  MultiTrajectory<source_link_t> fittedStates;

  // This is the index of the 'tip' of the track stored in multitrajectory.
  // Since this KF only stores one trajectory, it is unambiguous.
  // SIZE_MAX is the start of a trajectory.
  size_t trackTip = SIZE_MAX;

  // The optional Parameters at the provided surface
  std::optional<BoundParameters> fittedParameters;

  // Counter for states with measurements
  size_t measurementStates = 0;

  // Counter for handled states
  size_t processedStates = 0;

  // Indicator if smoothing has been done.
  bool smoothed = false;

  // Indicator if initialization has been performed.
  bool initialized = false;

  // Measurement surfaces without hits
  std::vector<const Surface*> missedActiveSurfaces;

  // Indicator if forward filtering has been done
  bool forwardFiltered = false;

  // Measurement surfaces handled in both forward and backward filtering
  std::vector<const Surface*> passedAgainSurfaces;

  Result<void> result{Result<void>::success()};
};

/// @brief Kalman fitter implementation of Acts as a plugin
///
/// to the Propgator
///
/// @tparam propagator_t Type of the propagation class
/// @tparam updater_t Type of the kalman updater class
/// @tparam smoother_t Type of the kalman smoother class
/// @tparam outlier_finder_t Type of the outlier finder class
/// @tparam calibrator_t Type of the calibrator class
/// @tparam input_converter_t Type of the input converter class
/// @tparam output_converter_t Type of the output converter class
///
/// The Kalman filter contains an Actor and a Sequencer sub-class.
/// The Sequencer has to be part of the Navigator of the Propagator
/// in order to initialize and provide the measurement surfaces.
///
/// The Actor is part of the Propagation call and does the Kalman update
/// and eventually the smoothing.  Updater, Smoother and Calibrator are
/// given to the Actor for further use:
/// - The Updater is the implemented kalman updater formalism, it
///   runs via a visitor pattern through the measurements.
/// - The Smoother is called at the end of the forward fit by the Actor.
/// - The outlier finder is called during the filtering by the Actor.
///   It determines if the measurement is an outlier
/// - The Calibrator is a dedicated calibration algorithm that allows
///   to calibrate measurements using track information, this could be
///    e.g. sagging for wires, module deformations, etc.
///
/// Measurements are not required to be ordered for the KalmanFilter,
/// measurement ordering needs to be figured out by the navigation of
/// the propagator.
///
/// The Input converter is a converter that transforms the input
/// measurement/track/segments into a set of FittableMeasurements
///
/// The Output converter is a converter that transforms the
/// set of track states into a given track/track particle class
///
/// The void components are provided mainly for unit testing.
template <typename propagator_t, typename updater_t = VoidKalmanUpdater,
          typename smoother_t = VoidKalmanSmoother,
          typename outlier_finder_t = VoidOutlierFinder,
          typename calibrator_t = VoidMeasurementCalibrator,
          typename input_converter_t = VoidKalmanComponents,
          typename output_converter_t = VoidKalmanComponents>
class KalmanFitter {
 public:
  /// Shorthand definition
  using MeasurementSurfaces = std::multimap<const Layer*, const Surface*>;

  /// Default constructor is deleted
  KalmanFitter() = delete;

  /// Constructor from arguments
  KalmanFitter(propagator_t pPropagator,
               std::unique_ptr<const Logger> logger =
                   getDefaultLogger("KalmanFilter", Logging::INFO),
               input_converter_t pInputCnv = input_converter_t(),
               output_converter_t pOutputCnv = output_converter_t())
      : m_propagator(std::move(pPropagator)),
        m_inputConverter(std::move(pInputCnv)),
        m_outputConverter(std::move(pOutputCnv)),
        m_logger(logger.release()) {}

 private:
  /// The propgator for the transport and material update
  propagator_t m_propagator;

  /// The input converter to Fittable measurements
  input_converter_t m_inputConverter;

  /// The output converter into a given format
  output_converter_t m_outputConverter;

  /// Logger getter to support macros
  const Logger& logger() const { return *m_logger; }

  /// Owned logging instance
  std::shared_ptr<const Logger> m_logger;

  /// The navigator type
  using KalmanNavigator = typename decltype(m_propagator)::Navigator;

  /// The navigator has DirectNavigator type or not
  static constexpr bool isDirectNavigator =
      std::is_same<KalmanNavigator, DirectNavigator>::value;

  /// @brief Propagator Actor plugin for the KalmanFilter
  ///
  /// @tparam source_link_t is an type fulfilling the @c SourceLinkConcept
  /// @tparam parameters_t The type of parameters used for "local" paremeters.
  ///
  /// The KalmanActor does not rely on the measurements to be
  /// sorted along the track.
  template <typename source_link_t, typename parameters_t>
  class Actor {
   public:
    using TrackStateType = TrackState<source_link_t, parameters_t>;

    /// Explicit constructor with updater and calibrator
    Actor(updater_t pUpdater = updater_t(), smoother_t pSmoother = smoother_t(),
          outlier_finder_t pOutlierFinder = outlier_finder_t(),
          calibrator_t pCalibrator = calibrator_t())
        : m_updater(std::move(pUpdater)),
          m_smoother(std::move(pSmoother)),
          m_outlierFinder(std::move(pOutlierFinder)),
          m_calibrator(std::move(pCalibrator)) {}

    /// Broadcast the result_type
    using result_type = KalmanFitterResult<source_link_t>;

    /// The target surface
    const Surface* targetSurface = nullptr;

    /// Allows retrieving measurements for a surface
    std::map<const Surface*, source_link_t> inputMeasurements;

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// Whether to consider energy loss.
    bool energyLoss = true;

    /// Whether run smoothing as backward filtering
    bool backwardFiltering = false;

    /// @brief Kalman actor operation
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
      ACTS_VERBOSE("KalmanFitter step");

      // Initialization:
      // - Only when track states are not set
      if (!result.initialized) {
        // -> Move the TrackState vector
        // -> Feed the KalmanSequencer with the measurements to be fitted
        ACTS_VERBOSE("Initializing");
        initialize(state, stepper, result);
        result.initialized = true;
      }

      // Update:
      // - Waiting for a current surface that has material
      // -> a trackState will be created on surface with material
      auto surface = state.navigation.currentSurface;
      if (surface and surface->surfaceMaterial()) {
        // Check if the surface is in the measurement map
        // -> Get the measurement / calibrate
        // -> Create the predicted state
        // -> Perform the kalman update
        // -> Check outlier behavior (@todo)
        // -> Fill strack state information & update stepper information
        if (state.stepping.navDir == forward and not result.smoothed and
            not result.forwardFiltered) {
          ACTS_VERBOSE("Perform forward filter step");
          auto res = filter(surface, state, stepper, result);
          if (!res.ok()) {
            ACTS_ERROR("Error in forward filter: " << res.error());
            result.result = res.error();
          }
        } else if (state.stepping.navDir == backward) {
          ACTS_VERBOSE("Perform backward filter step");
          auto res = backwardFilter(surface, state, stepper, result);
          if (!res.ok()) {
            ACTS_ERROR("Error in backward filter: " << res.error());
            result.result = res.error();
          }
        }
      }

      // Finalization:
      // when all track states have been handled or the navigation is breaked,
      // reset navigation&stepping before run backward filtering or
      // proceed to run smoothing
      if (state.stepping.navDir == forward) {
        if (result.measurementStates == inputMeasurements.size() or
            (result.measurementStates > 0 and
             state.navigation.navigationBreak)) {
          if (backwardFiltering and not result.forwardFiltered) {
            ACTS_VERBOSE("Forward filtering done");
            result.forwardFiltered = true;
            // Run backward filtering
            // Reverse navigation direction and reset navigation and stepping
            // state to last measurement
            ACTS_VERBOSE("Reverse navigation direction.");
            reverse(state, stepper, result);
          } else if (not result.smoothed) {
            // -> Sort the track states (as now the path length is set)
            // -> Call the smoothing
            // -> Set a stop condition when all track states have been handled
            ACTS_VERBOSE("Finalize/run smoothing");
            auto res = finalize(state, stepper, result);
            if (!res.ok()) {
              ACTS_ERROR("Error in finalize: " << res.error());
              result.result = res.error();
            }
          }
        }
      }

      // Post-finalization:
      // - Progress to target/reference surface and built the final track
      // parameters
      if ((result.smoothed or state.stepping.navDir == backward) and
          targetReached(state, stepper, *targetSurface)) {
        ACTS_VERBOSE("Completing");
        // Transport & bind the parameter to the final surface
        auto fittedState =
            stepper.boundState(state.stepping, *targetSurface, true);
        // Assign the fitted parameters
        result.fittedParameters = std::get<BoundParameters>(fittedState);
        // Break the navigation for stopping the Propagation
        state.navigation.navigationBreak = true;

        // Reset smoothed status of states missed in backward filtering
        if (backwardFiltering) {
          result.fittedStates.applyBackwards(result.trackTip, [&](auto state) {
            auto fSurface = &state.referenceSurface();
            auto surface_it = std::find_if(
                result.passedAgainSurfaces.begin(),
                result.passedAgainSurfaces.end(),
                [=](const Surface* surface) { return surface == fSurface; });
            if (surface_it == result.passedAgainSurfaces.end()) {
              // If backward filtering missed this surface, then there is no
              // smoothed parameter
              state.data().ismoothed = detail_lt::IndexData::kInvalid;
            }
          });
        }
      }
    }

    /// @brief Kalman actor operation : initialize
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    void initialize(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    result_type& /*result*/) const {}

    /// @brief Kalman actor operation : reverse direction
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    void reverse(propagator_state_t& state, stepper_t& stepper,
                 result_type& result) const {
      // Reset propagator options
      state.options.direction = backward;
      state.options.maxStepSize = -1.0 * state.options.maxStepSize;
      // Not sure if reset of pathLimit during propagation makes any sense
      state.options.pathLimit = -1.0 * state.options.pathLimit;

      // Reset stepping&navigation state using last measurement track state on
      // sensitive surface
      state.navigation = typename propagator_t::NavigatorState();
      result.fittedStates.applyBackwards(result.trackTip, [&](auto st) {
        if (st.hasUncalibrated()) {
          // Set the navigation state
          state.navigation.startSurface = &st.referenceSurface();
          state.navigation.startLayer =
              state.navigation.startSurface->associatedLayer();
          state.navigation.startVolume =
              state.navigation.startLayer->trackingVolume();
          state.navigation.targetSurface = targetSurface;
          state.navigation.currentSurface = state.navigation.startSurface;
          state.navigation.currentVolume = state.navigation.startVolume;

          // Update the stepping state
          stepper.update(state.stepping,
                         st.filteredParameters(state.options.geoContext));
          // Reverse stepping direction
          state.stepping.navDir = backward;
          state.stepping.stepSize = ConstrainedStep(state.options.maxStepSize);
          state.stepping.pathAccumulated = 0.;
          // Reinitialize the stepping jacobian
          st.referenceSurface().initJacobianToGlobal(
              state.options.geoContext, state.stepping.jacToGlobal,
              state.stepping.pos, state.stepping.dir,
              st.filteredParameters(state.options.geoContext).parameters());
          state.stepping.jacobian = BoundMatrix::Identity();
          state.stepping.jacTransport = FreeMatrix::Identity();
          state.stepping.derivative = FreeVector::Zero();

          // For the last measurement state, smoothed is filtered
          st.smoothed() = st.filtered();
          st.smoothedCovariance() = st.filteredCovariance();
          result.passedAgainSurfaces.push_back(&st.referenceSurface());

          // Update material effects for last measurement state in backward
          // direction
          materialInteractor(state.navigation.currentSurface, state, stepper);

          return false;  // abort execution
        }
        return true;  // continue execution
      });
    }

    /// @brief Kalman actor operation : update
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
      // Try to find the surface in the measurement surfaces
      auto sourcelink_it = inputMeasurements.find(surface);
      if (sourcelink_it != inputMeasurements.end()) {
        // Screen output message
        ACTS_VERBOSE("Measurement surface " << surface->geoID()
                                            << " detected.");

        // Update state and stepper with pre material effects
        materialInteractor(surface, state, stepper, preUpdate);

        // Transport & bind the state to the current surface
        auto [boundParams, jacobian, pathLength] =
            stepper.boundState(state.stepping, *surface, true);

        // add a full TrackState entry multi trajectory
        // (this allocates storage for all components, we will set them later)
        result.trackTip = result.fittedStates.addTrackState(
            TrackStatePropMask::All, result.trackTip);

        // now get track state proxy back
        auto trackStateProxy =
            result.fittedStates.getTrackState(result.trackTip);

        // assign the source link to the track state
        trackStateProxy.uncalibrated() = sourcelink_it->second;

        // Fill the track state
        trackStateProxy.predicted() = boundParams.parameters();
        trackStateProxy.predictedCovariance() = *boundParams.covariance();
        trackStateProxy.jacobian() = jacobian;
        trackStateProxy.pathLength() = pathLength;

        // We have predicted parameters, so calibrate the uncalibrated input
        // measuerement
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

        // If the update is successful, set covariance and
        auto updateRes = m_updater(state.geoContext, trackStateProxy, forward);
        if (!updateRes.ok()) {
          ACTS_ERROR("Update step failed: " << updateRes.error());
          return updateRes.error();
        } else {
          if (not m_outlierFinder(trackStateProxy)) {
            // Set the measurement type flag
            typeFlags.set(TrackStateFlag::MeasurementFlag);
            // Update the stepping state with filtered parameters
            ACTS_VERBOSE(
                "Filtering step successful, updated parameters are : \n"
                << trackStateProxy.filtered().transpose());
            // update stepping state using filtered parameters after kalman
            // update We need to (re-)construct a BoundParameters instance here,
            // which is a bit awkward.
            stepper.update(state.stepping, trackStateProxy.filteredParameters(
                                               state.options.geoContext));
            // We count the state with measurement
            ++result.measurementStates;
          } else {
            ACTS_VERBOSE(
                "Filtering step successful. But measurement is deterimined to "
                "be an outlier. Stepping state is not updated.")
            // Set the outlier type flag
            typeFlags.set(TrackStateFlag::OutlierFlag);
          }

          // Update state and stepper with post material effects
          materialInteractor(surface, state, stepper, postUpdate);
        }
        // We count the processed state
        ++result.processedStates;
      } else {
        if (result.measurementStates > 0) {
          // No source links on surface, add either hole or passive material
          // TrackState entry multi trajectory. No storage allocation for
          // uncalibrated/calibrated measurement and filtered parameter
          result.trackTip = result.fittedStates.addTrackState(
              ~(TrackStatePropMask::Uncalibrated |
                TrackStatePropMask::Calibrated | TrackStatePropMask::Filtered),
              result.trackTip);

          // now get track state proxy back
          auto trackStateProxy =
              result.fittedStates.getTrackState(result.trackTip);

          // Set the surface
          trackStateProxy.setReferenceSurface(surface->getSharedPtr());

          // Set the track state flags
          auto& typeFlags = trackStateProxy.typeFlags();
          typeFlags.set(TrackStateFlag::MaterialFlag);
          typeFlags.set(TrackStateFlag::ParameterFlag);

          if (surface->associatedDetectorElement() != nullptr) {
            ACTS_VERBOSE("Detected hole on " << surface->geoID());
            // If the surface is sensitive, set the hole type flag
            typeFlags.set(TrackStateFlag::HoleFlag);

            // Count the missed surface
            result.missedActiveSurfaces.push_back(surface);

            // Transport & bind the state to the current surface
            auto [boundParams, jacobian, pathLength] =
                stepper.boundState(state.stepping, *surface, true);

            // Fill the track state
            trackStateProxy.predicted() = boundParams.parameters();
            trackStateProxy.predictedCovariance() = *boundParams.covariance();
            trackStateProxy.jacobian() = jacobian;
            trackStateProxy.pathLength() = pathLength;
          } else {
            ACTS_VERBOSE("Detected in-sensitive surface " << surface->geoID());

            // Transport & get curvilinear state instead of bound state
            auto [curvilinearParams, jacobian, pathLength] =
                stepper.curvilinearState(state.stepping, true);

            // Fill the track state
            trackStateProxy.predicted() = curvilinearParams.parameters();
            trackStateProxy.predictedCovariance() =
                *curvilinearParams.covariance();
            trackStateProxy.jacobian() = jacobian;
            trackStateProxy.pathLength() = pathLength;
          }

          // Set the filtered parameter index to be the same with predicted
          // parameter
          trackStateProxy.data().ifiltered = trackStateProxy.data().ipredicted;

          // We count the processed state
          ++result.processedStates;
        }

        // Update state and stepper with material effects
        materialInteractor(surface, state, stepper, fullUpdate);
      }
      return Result<void>::success();
    }

    /// @brief Kalman actor operation : backward update
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result The mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    Result<void> backwardFilter(const Surface* surface,
                                propagator_state_t& state,
                                const stepper_t& stepper,
                                result_type& result) const {
      // Try to find the surface in the measurement surfaces
      auto sourcelink_it = inputMeasurements.find(surface);
      if (sourcelink_it != inputMeasurements.end()) {
        // Screen output message
        ACTS_VERBOSE("Measurement surface "
                     << surface->geoID()
                     << " detected in backward propagation.");

        // No backward filtering for last measurement state, but still update
        // with material effects
        if (state.stepping.navDir == backward and
            surface == state.navigation.startSurface) {
          materialInteractor(surface, state, stepper);
          return Result<void>::success();
        }

        // Update state and stepper with pre material effects
        materialInteractor(surface, state, stepper, preUpdate);

        // Transport & bind the state to the current surface
        auto [boundParams, jacobian, pathLength] =
            stepper.boundState(state.stepping, *surface, true);

        // Create a detached track state proxy
        auto tempTrackTip =
            result.fittedStates.addTrackState(TrackStatePropMask::All);

        // Get the detached track state proxy back
        auto trackStateProxy = result.fittedStates.getTrackState(tempTrackTip);

        // Assign the source link to the detached track state
        trackStateProxy.uncalibrated() = sourcelink_it->second;

        // Fill the track state
        trackStateProxy.predicted() = boundParams.parameters();
        trackStateProxy.predictedCovariance() = *boundParams.covariance();
        trackStateProxy.jacobian() = jacobian;
        trackStateProxy.pathLength() = pathLength;

        // We have predicted parameters, so calibrate the uncalibrated input
        // measuerement
        std::visit(
            [&](const auto& calibrated) {
              trackStateProxy.setCalibrated(calibrated);
            },
            m_calibrator(trackStateProxy.uncalibrated(),
                         trackStateProxy.predicted()));

        // If the update is successful, set covariance and
        auto updateRes = m_updater(state.geoContext, trackStateProxy, backward);
        if (!updateRes.ok()) {
          ACTS_ERROR("Backward update step failed: " << updateRes.error());
          return updateRes.error();
        } else {
          // Update the stepping state with filtered parameters
          ACTS_VERBOSE(
              "Backward Filtering step successful, updated parameters are : \n"
              << trackStateProxy.filtered().transpose());

          // Fill the smoothed parameter for the existing track state
          result.fittedStates.applyBackwards(result.trackTip, [&](auto state) {
            auto fSurface = &state.referenceSurface();
            if (fSurface == surface) {
              result.passedAgainSurfaces.push_back(surface);
              state.smoothed() = trackStateProxy.filtered();
              state.smoothedCovariance() = trackStateProxy.filteredCovariance();
              return false;
            }
            return true;
          });

          // update stepping state using filtered parameters after kalman update
          // We need to (re-)construct a BoundParameters instance here, which is
          // a bit awkward.
          stepper.update(state.stepping, trackStateProxy.filteredParameters(
                                             state.options.geoContext));

          // Update state and stepper with post material effects
          materialInteractor(surface, state, stepper, postUpdate);
        }
      } else {
        // Transport covariance
        if (surface->associatedDetectorElement() != nullptr) {
          ACTS_VERBOSE("Detected hole on " << surface->geoID()
                                           << " in backward filtering");
          if (state.stepping.covTransport) {
            stepper.covarianceTransport(state.stepping, *surface, true);
          }
        } else {
          ACTS_VERBOSE("Detected in-sensitive surface "
                       << surface->geoID() << " in backward filtering");
          if (state.stepping.covTransport) {
            stepper.covarianceTransport(state.stepping, true);
          }
        }

        // Not creating bound state here, so need manually reinitialize jacobian
        state.stepping.jacobian = BoundMatrix::Identity();

        // Update state and stepper with material effects
        materialInteractor(surface, state, stepper);
      }

      return Result<void>::success();
    }

    /// @brief Kalman actor operation : material interaction
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
      // Prepare relevant input particle properties
      detail::PointwiseMaterialInteraction interaction(surface, state, stepper);

      // Evaluate the material properties
      if (interaction.evaluateMaterialProperties(state, updateStage)) {
        // Evaluate the material effects
        interaction.evaluatePointwiseMaterialInteraction(multipleScattering,
                                                         energyLoss);

        ACTS_VERBOSE("Material effects on surface: "
                     << surface->geoID() << " at update stage: " << updateStage
                     << " are :");
        ACTS_VERBOSE("eLoss = "
                     << interaction.Eloss << ", "
                     << "variancePhi = " << interaction.variancePhi << ", "
                     << "varianceTheta = " << interaction.varianceTheta << ", "
                     << "varianceQoverP = " << interaction.varianceQoverP);

        // Update the state and stepper with material effects
        interaction.updateState(state, stepper);
      } else {
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
      // Remember you smoothed the track states
      result.smoothed = true;

      // Get the index of measurement states;
      std::vector<size_t> measurementIndices;
      auto lastState = result.fittedStates.getTrackState(result.trackTip);
      if (lastState.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
        measurementIndices.push_back(result.trackTip);
      }
      // Count track states to be smoothed
      size_t nStates = 0;
      result.fittedStates.applyBackwards(result.trackTip, [&](auto st) {
        // Smoothing will start from the last measurement state
        if (measurementIndices.empty()) {
          // No smoothed parameter for the last few non-measurment states
          st.data().ismoothed = detail_lt::IndexData::kInvalid;
        } else {
          nStates++;
        }
        size_t iprevious = st.previous();
        if (iprevious != Acts::detail_lt::IndexData::kInvalid) {
          auto previousState = result.fittedStates.getTrackState(iprevious);
          if (previousState.typeFlags().test(
                  Acts::TrackStateFlag::MeasurementFlag)) {
            measurementIndices.push_back(iprevious);
          }
        }
      });
      // Return error if the track has no measurement states (but this should
      // not happen)
      if (measurementIndices.empty()) {
        return KalmanFitterError::SmoothFailed;
      }
      // Screen output for debugging
      if (logger().doPrint(Logging::VERBOSE)) {
        ACTS_VERBOSE("Apply smoothing on " << nStates
                                           << " filtered track states.");
      }
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
      parameters_t smoothedPars =
          firstMeasurement.smoothedParameters(state.options.geoContext);

      // Update the stepping parameters - in order to progress to destination
      ACTS_VERBOSE(
          "Smoothing successful, updating stepping state, "
          "set target surface.");
      stepper.update(state.stepping, smoothedPars);
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

    /// Pointer to a logger that is owned by the parent, KalmanFilter
    const Logger* m_logger;

    /// Getter for the logger, to support logging macros
    const Logger& logger() const { return *m_logger; }

    /// The Kalman updater
    updater_t m_updater;

    /// The Kalman smoother
    smoother_t m_smoother;

    /// The outlier finder
    outlier_finder_t m_outlierFinder;

    /// The Measuremetn calibrator
    calibrator_t m_calibrator;

    /// The Surface beeing
    detail::SurfaceReached targetReached;
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
  /// @tparam source_link_t Source link type identifying uncalibrated input
  /// measurements.
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam kalman_fitter_options_t Type of the kalman fitter options
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param kfOptions KalmanOptions steering the fit
  /// @note The input measurements are given in the form of @c SourceLinks. It's
  /// @c calibrator_t's job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_t, typename start_parameters_t,
            typename kalman_fitter_options_t,
            typename parameters_t = BoundParameters,
            typename result_t = Result<KalmanFitterResult<source_link_t>>>
  auto fit(const std::vector<source_link_t>& sourcelinks,
           const start_parameters_t& sParameters,
           const kalman_fitter_options_t& kfOptions) const
      -> std::enable_if_t<!isDirectNavigator, result_t> {
    static_assert(SourceLinkConcept<source_link_t>,
                  "Source link does not fulfill SourceLinkConcept");

    static_assert(
        std::is_same<outlier_finder_t,
                     typename kalman_fitter_options_t::OutlierFinder>::value,
        "Inconsistent type of outlier finder between kalman fitter and "
        "kalman fitter options");

    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_VERBOSE("Preparing " << sourcelinks.size() << " input measurements");
    std::map<const Surface*, source_link_t> inputMeasurements;
    for (const auto& sl : sourcelinks) {
      const Surface* srf = &sl.referenceSurface();
      inputMeasurements.emplace(srf, sl);
    }

    // Create the ActionList and AbortList
    using KalmanAborter = Aborter<source_link_t, parameters_t>;
    using KalmanActor = Actor<source_link_t, parameters_t>;
    using KalmanResult = typename KalmanActor::result_type;
    using Actors = ActionList<KalmanActor>;
    using Aborters = AbortList<KalmanAborter>;

    // Create relevant options for the propagation options
    PropagatorOptions<Actors, Aborters> kalmanOptions(
        kfOptions.geoContext, kfOptions.magFieldContext);

    // Catch the actor and set the measurements
    auto& kalmanActor = kalmanOptions.actionList.template get<KalmanActor>();
    kalmanActor.m_logger = m_logger.get();
    kalmanActor.inputMeasurements = std::move(inputMeasurements);
    kalmanActor.targetSurface = kfOptions.referenceSurface;
    kalmanActor.multipleScattering = kfOptions.multipleScattering;
    kalmanActor.energyLoss = kfOptions.energyLoss;
    kalmanActor.backwardFiltering = kfOptions.backwardFiltering;

    // Set config for outlier finder
    kalmanActor.m_outlierFinder = kfOptions.outlierFinder;

    // also set logger on updater and smoother
    kalmanActor.m_updater.m_logger = m_logger;
    kalmanActor.m_smoother.m_logger = m_logger;

    // Run the fitter
    auto result = m_propagator.template propagate(sParameters, kalmanOptions);

    if (!result.ok()) {
      return result.error();
    }

    const auto& propRes = *result;

    /// Get the result of the fit
    auto kalmanResult = propRes.template get<KalmanResult>();

    /// It could happen that the fit ends in zero processed states.
    /// The result gets meaningless so such case is regarded as fit failure.
    //@TODO: should we require the number of measurments >0 ?
    if (kalmanResult.result.ok() and not kalmanResult.processedStates) {
      kalmanResult.result = Result<void>(KalmanFitterError::PropagationInVain);
    }

    if (!kalmanResult.result.ok()) {
      return kalmanResult.result.error();
    }

    // Return the converted Track
    return m_outputConverter(std::move(kalmanResult));
  }

  /// Fit implementation of the foward filter, calls the
  /// the forward filter and backward smoother
  ///
  /// @tparam source_link_t Source link type identifying uncalibrated input
  /// measurements.
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam kalman_fitter_options_t Type of the kalman fitter options
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param kfOptions KalmanOptions steering the fit
  /// @param sSequence surface sequence used to initialize a DirectNavigator
  /// @note The input measurements are given in the form of @c SourceLinks. It's
  /// @c calibrator_t's job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_t, typename start_parameters_t,
            typename kalman_fitter_options_t,
            typename parameters_t = BoundParameters,
            typename result_t = Result<KalmanFitterResult<source_link_t>>>
  auto fit(const std::vector<source_link_t>& sourcelinks,
           const start_parameters_t& sParameters,
           const kalman_fitter_options_t& kfOptions,
           const std::vector<const Surface*>& sSequence) const
      -> std::enable_if_t<isDirectNavigator, result_t> {
    static_assert(SourceLinkConcept<source_link_t>,
                  "Source link does not fulfill SourceLinkConcept");

    static_assert(
        std::is_same<outlier_finder_t,
                     typename kalman_fitter_options_t::OutlierFinder>::value,
        "Inconsistent type of outlier finder between kalman fitter and "
        "kalman fitter options");

    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_VERBOSE("Preparing " << sourcelinks.size() << " input measurements");
    std::map<const Surface*, source_link_t> inputMeasurements;
    for (const auto& sl : sourcelinks) {
      const Surface* srf = &sl.referenceSurface();
      inputMeasurements.emplace(srf, sl);
    }

    // Create the ActionList and AbortList
    using KalmanAborter = Aborter<source_link_t, parameters_t>;
    using KalmanActor = Actor<source_link_t, parameters_t>;
    using KalmanResult = typename KalmanActor::result_type;
    using Actors = ActionList<DirectNavigator::Initializer, KalmanActor>;
    using Aborters = AbortList<KalmanAborter>;

    // Create relevant options for the propagation options
    PropagatorOptions<Actors, Aborters> kalmanOptions(
        kfOptions.geoContext, kfOptions.magFieldContext);

    // Catch the actor and set the measurements
    auto& kalmanActor = kalmanOptions.actionList.template get<KalmanActor>();
    kalmanActor.m_logger = m_logger.get();
    kalmanActor.inputMeasurements = std::move(inputMeasurements);
    kalmanActor.targetSurface = kfOptions.referenceSurface;
    kalmanActor.multipleScattering = kfOptions.multipleScattering;
    kalmanActor.energyLoss = kfOptions.energyLoss;
    kalmanActor.backwardFiltering = kfOptions.backwardFiltering;

    // Set config for outlier finder
    kalmanActor.m_outlierFinder.m_config = kfOptions.outlierFinderConfig;

    // also set logger on updater and smoother
    kalmanActor.m_updater.m_logger = m_logger;
    kalmanActor.m_smoother.m_logger = m_logger;

    // Set the surface sequence
    auto& dInitializer =
        kalmanOptions.actionList.template get<DirectNavigator::Initializer>();
    dInitializer.surfaceSequence = sSequence;

    // Run the fitter
    auto result = m_propagator.template propagate(sParameters, kalmanOptions);

    if (!result.ok()) {
      return result.error();
    }

    const auto& propRes = *result;

    /// Get the result of the fit
    auto kalmanResult = propRes.template get<KalmanResult>();

    /// It could happen that the fit ends in zero processed states.
    /// The result gets meaningless so such case is regarded as fit failure.
    if (kalmanResult.result.ok() and not kalmanResult.processedStates) {
      kalmanResult.result = Result<void>(KalmanFitterError::PropagationInVain);
    }

    if (!kalmanResult.result.ok()) {
      return kalmanResult.result.error();
    }

    // Return the converted Track
    return m_outputConverter(std::move(kalmanResult));
  }
};

}  // namespace Acts
