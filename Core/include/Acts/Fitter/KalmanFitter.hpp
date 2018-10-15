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
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/surface_getter.hpp"
#include "Acts/EventData/detail/trackstate_manipulation.hpp"
#include "Acts/EventData/detail/trackstate_sorters.hpp"
#include "Acts/Fitter/detail/VoidKalmanComponents.hpp"
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
    /// Explicit constructor with updagor and calibrator
    Actor(updator_t    pUpdator    = updator_t(),
          smoother_t   pSmoother   = smoother_t(),
          calibrator_t pCalibrator = calibrator_t())
      : m_updator(std::move(pUpdator))
      , m_smoother(std::move(pSmoother))
      , m_calibrator(std::move(pCalibrator))
    {
    }

    /// Simple result struct to be returned
    /// It mainly acts as an internal state state which is
    /// created for every propagation/extrapolation step
    struct this_result
    {
      // Move the result into the fitted states
      track_states_t fittedStates;

      // Measurement surfaces without hits
      std::vector<const Surface*> missedActiveSurfaces = {};

      // Counter for handled states
      size_t processedStates = 0;

      // The index map for accessing the track state in order
      std::map<const Surface*, size_t> accessIndex;
    };
    using result_type = this_result;

    /// The Track states with which the Actor is initialized
    track_states_t trackStates = {};

    /// @brief Kalman actor operation
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    ///
    /// @param state is the mutable propagator state object
    /// @param result is the mutable result state object
    template <typename propagator_state_t>
    void
    operator()(propagator_state_t& state, result_type& result) const
    {
      // Initialization:
      // - Only when track states are not set
      if (result.fittedStates.empty()) {
        // -> Move the TrackState vector
        // -> Feed the KalmanSequencer with the measurements to be fitted
        initialize(state, result);
      }

      // Update:
      // - Waiting for a current surface that appears in the measurement list
      auto surface = state.navigation.currentSurface;
      if (surface) {
        // Check if the surface is in the measurement map
        // -> Get the measurement / calibrate
        // -> Create the predicted state
        // -> Perform the kalman update
        // -> Check outlier behavior (@todo)
        // -> Fill strack state information & update stepper information
        filter(surface, state, result);
      }

      // Finalization:
      // - When all track states have been handled
      if (result.processedStates == trackStates.size()) {
        // -> Sort the track states (as now the path length is set)
        // -> Call the smoothing
        // -> Set a stop condition when all track states have been handled
        finalize(state, result);
      }
    }

  private:
    /// @brief Kalman actor operation : initialize
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    ///
    /// @param state is the mutable propagator state object
    /// @param result is the mutable result state object
    template <typename propagator_state_t>
    void
    initialize(propagator_state_t& state, result_type& result) const
    {
      // Create the multimap
      MeasurementSurfaces measurementSurfaces;
      // Move the track states
      result.fittedStates = std::move(trackStates);
      // @todo -----> outsource this to a assign to layer method

      // Memorize the index to access the state
      size_t stateIndex = 0;
      for (auto& tState : result.fittedStates) {
        // Get the Surface
        auto surface = &(detail::getSurface(tState));
        // Get the associated Layer to this Surface
        auto layer = surface->associatedLayer();
        if (layer == nullptr) {
          // Find the intersection to allocate the layer
          auto surfaceIntersection
              = surface->intersectionEstimate(state.stepping.position(),
                                              state.stepping.direction(),
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
              std::pair<const Layer*, const Surface*>(layer, surface));
          // Insert the fitted state into the fittedStates map
          result.accessIndex[surface] = stateIndex;
        }
        ++stateIndex;
      }
      // Feed the KalmanSequencer with the measurement surfaces
      state.navigation.sequence.externalSurfaces
          = std::move(measurementSurfaces);
    }

    /// @brief Kalman actor operation : update
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param result The mutable result state object
    template <typename propagator_state_t>
    void
    filter(const Surface*      surface,
           propagator_state_t& state,
           result_type&        result) const
    {
      // Try to find the surface in the measurement surfaces
      auto cindexItr = result.accessIndex.find(surface);
      if (cindexItr != result.accessIndex.end()) {
        // Transport & bind the stateto the current surface
        auto boundState = state.stepping.bind(*surface, true);
        // Get the current VariantTrackState
        auto trackState = result.fittedStates[cindexItr->second];
        // Perform the update and obtain the filtered parameters
        auto filteredPars = m_updator(trackState, std::move(boundState));
        // If the update is successful, set covariance and
        if (filteredPars) {
          state.stepping.update(filteredPars.get());
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
    ///
    /// @param state is the mutable propagator state object
    /// @param result is the mutable result state object
    template <typename propagator_state_t>
    void
    finalize(propagator_state_t& state, result_type& result) const
    {
      // Sort the TrackStates according to the path length
      detail::path_length_sorter plSorter;
      std::sort(
          result.fittedStates.begin(), result.fittedStates.end(), plSorter);
      // Smooth the track states and obtain the last smoothed track parameters
      auto smoothedPars = m_smoother(result.fittedStates);
      // Update the stepping parameters - in order to progress to destination
      if (smoothedPars) {
        state.stepping.update(smoothedPars.get());
      }
      // Set the destination surface - we should re-do the navigation
    }

    /// The Kalman updator
    updator_t m_updator;

    /// The Kalman smoother
    smoother_t m_smoother;

    /// The Measuremetn calibrator
    calibrator_t m_calibrator;
  };

  /// Default constructor is deleted
  KalmanFitter() = delete;

  /// Constructor from arguments
  KalmanFitter(propagator_t       pPropagator,
               calibrator_t       pCalibrator = calibrator_t(),
               input_converter_t  pInputCnv   = input_converter_t(),
               output_converter_t pOutputCnv  = output_converter_t())
    : m_propagator(std::move(pPropagator))
    , m_calibrator(std::move(pCalibrator))
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
  /// @param measurements are the fittable measurements
  /// @param initalParameters is the initial track parameters
  ///
  /// @return the output as an output track
  template <typename input_measurements_t,
            typename parameters_t,
            typename surface_t>
  auto
  fit(const input_measurements_t& measurements,
      const parameters_t& /*initalParameters*/,
      const surface_t* /*pReferenceSurface = nullptr*/) const
  {
    // Bring the measurements into Acts style
    auto trackStates = m_inputConverter(measurements);

    // Return the converted Track
    return m_outputConverter(trackStates);
  }

private:
  /// The propgator for the transport and material update
  propagator_t m_propagator;

  /// The measurement calibrator
  calibrator_t m_calibrator;

  /// The input converter to Fittable measurements
  input_converter_t m_inputConverter;

  /// The output converter into a given format
  output_converter_t m_outputConverter;
};

}  // namespace Acts
