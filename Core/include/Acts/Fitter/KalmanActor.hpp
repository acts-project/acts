// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/EventData/detail/surface_getter.hpp"
#include "Acts/Layers/Layer.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

/// @brief Propagator Actor plugin for the KalmanFilter
///
template <typename track_states_t>
class KalmanActor
{
public:
  using MeasurementSurfaces = std::multimap<const Layer*, const Surface*>;

  /// Simple result struct to be returned
  /// It mainly acts as an internal state state which is
  /// created for every propagation/extrapolation step
  struct this_result
  {
    // Move the result into the fitted states
    track_states_t fittedStates;
    // The index map for accessing the track state in order
    std::map<const Surface*, size_t> accessIndex;
  };

  using result_type = this_result;

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
    // Initialization of the Actor:
    // Move the TrackState vector and feed the Sequencer
    // with the measurements to be fitted
    if (result.fittedStates.empty()) {
      initialize(state, result);
    }

    // Waiting for a current surface that appears in the measurement list
    auto surface = state.navigation.currentSurface;
    if (surface) {
      auto cindexItr = result.accessIndex.find(surface);
      if (cindexItr != result.accessIndex.end()) {
        // The current index & thus the current measurement
        auto cmeasurement = result.fittedStates[cindexItr->second];
        // Create the predicted state
        state.stepping.covarianceTransport(*surface, true);
        auto jacobian = state.stepping.jacobian;
        // Prediction
        auto covPtr
            = std::make_unique<const ActsMatrixD<5, 5>>(state.stepping.cov);
        BoundParameters predicted(std::move(covPtr),
                                  state.stepping.pos,
                                  state.stepping.p * state.stepping.dir,
                                  state.stepping.q,
                                  *surface);
      }
    }
  }

private:
  /// @brief Kalman actor operation
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
    state.navigation.sequence.externalSurfaces = std::move(measurementSurfaces);
  }

  /// The private navigation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the
  /// state.options.debug == true
  /// case in order not to spend time when not needed.
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param[in,out] state the propagator state for the debug flag,
  ///      prefix and length
  /// @param logAction is a callable function that returns a stremable object
  template <typename propagator_state_t>
  void
  debugLog(propagator_state_t&                 state,
           const std::function<std::string()>& logAction) const
  {
    if (state.options.debug) {
      std::stringstream dstream;
      std::string       ds = (state.stepping.navDir == forward) ? "K->" : "<-K";
      dstream << ds << std::setw(state.options.debugPfxWidth);
      dstream << "KalmanActor"
              << " | ";
      dstream << std::setw(state.options.debugMsgWidth) << logAction() << '\n';
      state.options.debugString += dstream.str();
    }
  }
};

}  // namespace Acts
