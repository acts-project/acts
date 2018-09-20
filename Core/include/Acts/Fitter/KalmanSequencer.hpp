// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <sstream>

namespace Acts {

/// @brief Extension to Acts::Navigator for Kalman fitting/filtering
///
/// This is implemented as an actor and should be provided
/// to the Propagator::Options<ActionList,AbortList> as one
/// of the abort list otions (ideally the first)
///
/// It will then set the navigator surface memory to "ON" and
/// enables the bounceback for the smoother
///
/// @note This works only with the Acts::Navigator
/// and obviously not with the Acts::VoidNavigator.
struct KalmanSequencer
{
  /// Simple result struct to be returned
  /// It has all the SurfaceHit objects that
  /// are collected (and thus have been selected)
  struct this_state
  {
    // Initialization
    bool initialized = false;
    // If Smoothing is on - remember you bounced
    bool bounced = false;
  };
  using state_type = this_state;

  /// Run the smoothing
  bool smoothing = true;

  /// Surface provider action for the ActionList of the Propagator
  ///
  /// @tparam propagator_state_t is the type of Propagator state
  ///
  /// @param[in,out] state is the mutable stepper state object
  /// @param[in,out] result is the mutable result object
  template <typename propagator_state_t>
  bool
  operator()(propagator_state_t& state) const
  {
    if (!state.navigation.sequence.initialized) {
      /// Now set the initialized flag to true
      state.navigation.sequence.initialized = true;
      /// Remember all processed surfaces (in smoothing mode) surfaces
      state.navigation.rememberStates = smoothing;
      debugLog(state, [&] {
        return std::string("Initialization - switching navigation memory on.");
      });
      /// Sequencer went into action
      return true;
    }

    // Check if we have to switch to smoothing mode - only do once though
    if ((state.navigation.targetReached || state.navigation.navigationBreak)
        && smoothing
        && !state.navigation.sequence.bounced) {

      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "End of forward filter - switching to smoother mode: "
                << state.navigation.processedSurfaces.size()
                << " surfaces to reverse-navigate.";
        return dstream.str();
      });

      if (!state.navigation.processedSurfaces.empty()) {
        debugLog(state, [&] {
          return std::string(
              "Reversing navigation information & step direction.");
        });
        // Unset the targetReached and navigationBreak flags, set reverseMode to
        // true
        state.navigation.targetReached   = false;
        state.navigation.navigationBreak = false;
        state.navigation.reverseMode     = true;
        // Set the start parameters to be the new target parameters
        state.navigation.targetSurface = state.navigation.startSurface;
        state.navigation.targetLayer   = state.navigation.startLayer;
        state.navigation.targetVolume  = state.navigation.startVolume;
        // Reverse the direction & (last) step size
        state.stepping.navDir
            = state.stepping.navDir == forward ? backward : forward;
        // The total Accumulated path so far
        double tPath = state.stepping.pathAccumulated;
        // Generate a reverse list of surface intersections from
        // the remembered surfaces
        auto& rSurfaces = state.navigation.processedSurfaces;
        std::for_each(rSurfaces.begin(),
                      rSurfaces.end(),
                      [&state, tPath](SurfaceIntersection& si) {
                        si.intersection.pathLength -= tPath;
                        si.pDirection = state.stepping.navDir;
                      });
        // Reverse it
        std::reverse(rSurfaces.begin(), rSurfaces.end());
        state.navigation.navSurfaces    = std::move(rSurfaces);
        state.navigation.navSurfaceIter = state.navigation.navSurfaces.begin();
      }
      state.stepping.stepSize
          = state.navigation.navSurfaceIter->intersection.pathLength;
      // For the moment set the current volume to the start volume
      // @todo remember volumes as well for volume based material update
      state.navigation.currentVolume = state.navigation.targetVolume;
      // Screen output
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Navigation stepSize updated to ";
        dstream << state.stepping.stepSize.toString();
        return dstream.str();
      });
      state.navigation.sequence.bounced = true;
      // The sequencer was in action
      return true;
    }
    return false;
  }

private:
  /// The private propagation debug logging
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
      std::string       ds = (state.stepping.navDir == forward) ? "K->" : "<-K";
      dstream << ds << std::setw(state.options.debugPfxWidth);
      dstream << "KalmanSequencer"
              << " | ";
      dstream << std::setw(state.options.debugMsgWidth) << logAction() << '\n';
      state.options.debugString += dstream.str();
      std::cout << dstream.str();
    }
  }
};

}  // namespace Acts
