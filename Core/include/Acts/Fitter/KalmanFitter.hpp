// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/variant.hpp>
#include <memory>
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackStateSorters.hpp"
#include "Acts/Fitter/detail/VoidKalmanComponents.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @brief Kalman fitter implementation of Acts as a plugin
/// to the Propgator
///
/// @tparam propagator_t Type of the propagation class
/// @tparam updator_t Type of the kalman updator class
/// @tparam smoother_t Type of the kalman smoother class
/// @tparam calibrator_t Type of the calibrator class
/// @tparam input_converter_t Type of the input converter class
/// @tparam output_converter_t Type of the output converter class
///
/// The Kalman filter contains an Actor and a Sequencer sub-class.
/// The Sequencer has to be part of the Navigator of the Propagator
/// in order to initialize and provide the measurement surfaces.
///
/// The Actor is part of the Propagation call and does the Kalman update
/// and eventually the smoothing.  Updator, Smoother and Calibrator are
/// given to the Actor for further use:
/// - The Updator is the implemented kalman updator formalism, it
///   runs via a visitor pattern through the measurements.
/// - The Smoother is called at the end of the forward fit by the Actor.
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
template <typename propagator_t,
          typename updator_t          = VoidKalmanUpdator,
          typename smoother_t         = VoidKalmanSmoother,
          typename calibrator_t       = VoidKalmanComponents,
          typename input_converter_t  = VoidKalmanComponents,
          typename output_converter_t = VoidKalmanComponents>
class KalmanFitter
{
public:
  /// Shorthand definition
  using MeasurementSurfaces = std::multimap<const Layer*, const Surface*>;

  /// Default constructor is deleted
  KalmanFitter() = delete;

  /// Constructor from arguments
  KalmanFitter(propagator_t       pPropagator,
               input_converter_t  pInputCnv  = input_converter_t(),
               output_converter_t pOutputCnv = output_converter_t())
    : m_propagator(std::move(pPropagator))
    , m_inputConverter(std::move(pInputCnv))
    , m_outputConverter(std::move(pOutputCnv))
  {
  }

  /// Fit implementation of the foward filter, calls the
  /// the forward filter and backward smoother
  ///
  /// @tparam input_measurements_t Type of the fittable measurements
  /// @tparam parameters_t Type of the initial parameters
  /// @tparam surface_t Type of the reference surface
  ///
  /// @param measurements The fittable measurements
  /// @param sParameters The initial track parameters
  /// @param rSurface The reference surface
  ///
  /// @return the output as an output track
  template <typename input_measurements_t,
            typename parameters_t,
            typename surface_t>
  auto
  fit(input_measurements_t measurements,
      const parameters_t&  sParameters,
      const surface_t*     rSurface = nullptr) const
  {
    // Bring the measurements into Acts style
    auto trackStates = m_inputConverter(measurements);

    // Create the ActionList and AbortList
    using KalmanActor  = Actor<decltype(trackStates)>;
    using KalmanResult = typename KalmanActor::result_type;
    using Actors       = ActionList<KalmanActor>;
    using Aborters     = AbortList<>;

    // Create relevant options for the propagation options
    PropagatorOptions<Actors, Aborters> kalmanOptions;
    // Catch the actor and set the measurements
    auto& kalmanActor = kalmanOptions.actionList.template get<KalmanActor>();
    kalmanActor.trackStates   = std::move(trackStates);
    kalmanActor.targetSurface = rSurface;

    // Run the fitter
    const auto& result
        = m_propagator.template propagate(sParameters, kalmanOptions);

    /// Get the result of the fit
    auto kalmanResult = result.template get<KalmanResult>();

    // Return the converted Track
    return m_outputConverter(std::move(kalmanResult));
  }

private:
  /// The propgator for the transport and material update
  propagator_t m_propagator;

  /// The input converter to Fittable measurements
  input_converter_t m_inputConverter;

  /// The output converter into a given format
  output_converter_t m_outputConverter;

  /// @brief Propagator Actor plugin for the KalmanFilter
  ///
  /// @tparam track_states_t is any iterable std::container of
  /// boost::variant TrackState objects.
  ///
  /// @tparam updator_t The Kalman updator used for this fitter
  ///
  /// @tparam calibrator_t The Measurement calibrator for Fittable
  /// measurements to be calibrated
  ///
  /// The KalmanActor does not rely on the measurements to be
  /// sorted along the track.
  template <typename track_states_t>
  class Actor
  {
  public:
    using TrackState = typename track_states_t::value_type;

    /// Explicit constructor with updator and calibrator
    Actor(updator_t    pUpdator    = updator_t(),
          smoother_t   pSmoother   = smoother_t(),
          calibrator_t pCalibrator = calibrator_t())
      : m_updator(std::move(pUpdator))
      , m_smoother(std::move(pSmoother))
      , m_calibrator(std::move(pCalibrator))
    {
    }

    /// Simple result struct to be returned
    /// It mainly acts as an internal state which is
    /// created for every propagation/extrapolation step
    struct this_result
    {
      // Move the result into the fitted states
      track_states_t fittedStates = {};

      // The optional Parameters at the provided surface
      boost::optional<BoundParameters> fittedParameters;

      // Counter for handled states
      size_t processedStates = 0;

      // Indicator if you smoothed
      bool smoothed = false;

      // Measurement surfaces without hits
      std::vector<const Surface*> missedActiveSurfaces = {};

      // The index map for accessing the track state in order
      std::map<const Surface*, size_t> accessIndices = {};
    };

    /// Broadcast the result_type
    using result_type = this_result;

    /// The target surface
    const Surface* targetSurface = nullptr;

    /// The Track states with which the Actor is initialized
    track_states_t trackStates = {};

    /// @brief Kalman actor operation
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    void
    operator()(propagator_state_t& state,
               const stepper_t&    stepper,
               result_type&        result) const
    {
      // Initialization:
      // - Only when track states are not set
      if (result.fittedStates.empty()) {
        // -> Move the TrackState vector
        // -> Feed the KalmanSequencer with the measurements to be fitted
        initialize(state, stepper, result);
      }

      // Update:
      // - Waiting for a current surface that appears in the measurement list
      auto surface = state.navigation.currentSurface;
      if (surface and not result.smoothed) {
        // Check if the surface is in the measurement map
        // -> Get the measurement / calibrate
        // -> Create the predicted state
        // -> Perform the kalman update
        // -> Check outlier behavior (@todo)
        // -> Fill strack state information & update stepper information
        filter(surface, state, stepper, result);
      }

      // Finalization:
      // - When all track states have been handled
      if (result.processedStates == trackStates.size()
          and not result.smoothed) {
        // -> Sort the track states (as now the path length is set)
        // -> Call the smoothing
        // -> Set a stop condition when all track states have been handled
        finalize(state, stepper, result);
      }
      // Post-finalization:
      // - Progress to target/reference surface and built the final track
      // parameters
      if (result.smoothed and targetReached(state, stepper, *targetSurface)) {
        // Transport & bind the parameter to the final surface
        auto fittedState
            = stepper.boundState(state.stepping, *targetSurface, true);
        // Assign the fitted parameters
        result.fittedParameters = std::get<BoundParameters>(fittedState);
        // Break the navigation for stopping the Propagation
        state.navigation.navigationBreak = true;
      }
    }

  private:
    /// @brief Kalman actor operation : initialize
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    void
    initialize(propagator_state_t& state,
               const stepper_t&    stepper,
               result_type&        result) const
    {
      // Screen output message
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Initializing KalmanFitter with ";
        dstream << trackStates.size();
        dstream << " measurements to fit.";
        return dstream.str();
      });
      // Create the multimap
      MeasurementSurfaces measurementSurfaces;
      // Move the track states
      result.fittedStates = std::move(trackStates);
      // Memorize the index to access the state
      size_t stateIndex = 0;
      for (auto& tState : result.fittedStates) {
        // Get the Surface
        const Surface& surface = tState.referenceSurface();
        // Get the associated Layer to this Surface
        auto layer = surface.associatedLayer();
        if (layer == nullptr) {
          // Find the intersection to allocate the layer
          auto surfaceIntersection
              = surface.intersectionEstimate(stepper.position(state.stepping),
                                             stepper.direction(state.stepping),
                                             state.stepping.navDir,
                                             false);
          // Allocate the layer via the tracking geometry search
          if (surfaceIntersection and state.navigation.worldVolume) {
            auto intersection = surfaceIntersection.position;
            auto layerVolume
                = state.navigation.worldVolume->trackingVolume(intersection);
            layer = layerVolume ? layerVolume->associatedLayer(intersection)
                                : nullptr;
          }
        }
        // Insert the surface into the measurementsurfaces multimap
        if (layer) {
          measurementSurfaces.insert(
              std::pair<const Layer*, const Surface*>(layer, &surface));
          // Insert the fitted state into the fittedStates map
          result.accessIndices[&surface] = stateIndex;
        }
        ++stateIndex;
      }
      // Screen output message
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Set ";
        dstream << measurementSurfaces.size();
        dstream << " measurements surfaces to the navigation.";
        return dstream.str();
      });
      // Feed the KalmanSequencer with the measurement surfaces
      state.navigation.externalSurfaces = std::move(measurementSurfaces);
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
    void
    filter(const Surface*      surface,
           propagator_state_t& state,
           const stepper_t&    stepper,
           result_type&        result) const
    {
      // Try to find the surface in the measurement surfaces
      auto cindexItr = result.accessIndices.find(surface);
      if (cindexItr != result.accessIndices.end()) {
        // Screen output message
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << "Measurement surface ";
          dstream << surface->geoID().toString();
          dstream << " detected.";
          return dstream.str();
        });

        // Get the current TrackState
        TrackState& trackState = result.fittedStates[cindexItr->second];

        // Transport & bind the state to the current surface
        std::tuple<BoundParameters,
                   typename TrackState::Parameters::CovMatrix_t,
                   double>
            boundState = stepper.boundState(state.stepping, *surface, true);

        trackState.parameter.predicted  = std::get<0>(boundState);
        trackState.parameter.jacobian   = std::get<1>(boundState);
        trackState.parameter.pathLength = std::get<2>(boundState);

        // If the update is successful, set covariance and
        if (m_updator(trackState)) {
          // Update the stepping state
          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << "Filtering step successful, updated parameters are : ";
            dstream << *trackState.parameter.filtered;
            return dstream.str();
          });
          // update stepping state using filtered parameters
          // after kalman update
          stepper.update(state.stepping, *trackState.parameter.filtered);
        }
        // We count the processed state
        ++result.processedStates;
      } else if (surface->associatedDetectorElement()) {
        // Count the missed surface
        result.missedActiveSurfaces.push_back(surface);
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
    void
    finalize(propagator_state_t& state,
             const stepper_t&    stepper,
             result_type&        result) const
    {
      // Remember you smoothed the track states
      result.smoothed = true;

      // Sort the TrackStates according to the path length
      TrackStatePathLengthSorter plSorter;
      std::sort(
          result.fittedStates.begin(), result.fittedStates.end(), plSorter);
      // Screen output for debugging
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Apply smoothing on ";
        dstream << result.fittedStates.size();
        dstream << " filtered track states.";
        return dstream.str();
      });
      // Smooth the track states and obtain the last smoothed track parameters
      const auto& smoothedPars = m_smoother(result.fittedStates);
      // Update the stepping parameters - in order to progress to destination
      if (smoothedPars) {
        // Update the stepping state
        debugLog(state, [&] {
          return std::string("Smoothing successful, updating stepping state, "
                             "set target surface.");
        });
        stepper.update(state.stepping, smoothedPars.get());
        // Reverse the propagation direction
        state.stepping.stepSize
            = detail::ConstrainedStep(-1. * state.options.maxStepSize);
        state.options.direction = backward;
      }
    }

    /// The private KalmanActor debug logging
    ///
    /// It needs to be fed by a lambda function that returns a string,
    /// that guarantees that the lambda is only called in the
    /// options.debug == true case in order not to spend time when not needed.
    ///
    /// @tparam propagator_state_t Type of the nested propagator state object
    ///
    /// @param state the propagator state for the debug flag, prefix/length
    /// @param logAction is a callable function that returns a stremable object
    template <typename propagator_state_t>
    void
    debugLog(propagator_state_t&                 state,
             const std::function<std::string()>& logAction) const
    {
      if (state.options.debug) {
        std::stringstream dstream;
        dstream << "K->" << std::setw(state.options.debugPfxWidth);
        dstream << "KalmanActor"
                << " | ";
        dstream << std::setw(state.options.debugMsgWidth) << logAction()
                << '\n';
        state.options.debugString += dstream.str();
      }
    }

    /// The Kalman updator
    updator_t m_updator;

    /// The Kalman smoother
    smoother_t m_smoother;

    /// The Measuremetn calibrator
    calibrator_t m_calibrator;

    /// The Surface beeing
    detail::SurfaceReached targetReached;
  };
};

}  // namespace Acts
