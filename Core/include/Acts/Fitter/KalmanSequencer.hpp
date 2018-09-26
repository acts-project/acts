// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <sstream>
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

/// @brief Extension to Acts::Navigator for Kalman fitting without
/// reference trajectory, i.e. it is designed for remembering the
/// forward navigation states and run them in reverse order
///
/// This is implemented as an actor and should be provided
/// to the Propagator::Options<ActionList,AbortList> as one
/// of the abort list otions (ideally the first)
///
/// It will then set the navigator surface memory to "ON" and
/// enables the bounceback for the smoother
///
/// @note This works only with the Acts::Navigator
/// and not with the Acts::VoidNavigator.
struct KalmanSequencer
{

  using NavigationSurfaces  = std::vector<SurfaceIntersection>;
  using MeasurementSurfaces = std::multimap<const Layer*, const Surface*>;

  /// Simple result struct to be returned
  /// It has all the SurfaceHit objects that
  /// are collected (and thus have been selected)
  struct this_state
  {

    /// If Smoothing is on - remember you bounced
    bool bounced = false;

    /// Measurement surfaces
    MeasurementSurfaces externalSurfaces = {};

    /// The memory of the sequence - these navigationSurfaces are ordered
    NavigationSurfaces navigationSurfaces = {};
  };
  using state_type = this_state;

  /// Run the smoothing
  bool smoothing = true;

  /// Surface provider action for the ActionList of the Propagator
  ///
  /// @tparam propagator_state_t is the type of Propagator state
  ///
  /// @param[in,out] state is the mutable stepper state object
  /// @param[in,out] memory Boolean flag to memorize surface
  template <typename propagator_state_t>
  bool
  operator()(propagator_state_t& state, bool memory = false) const
  {

    // Sequence: Reverse navigation
    // Check if we have to switch to smoothing mode - only do once though
    if ((state.navigation.targetReached or state.navigation.navigationBreak)
        and smoothing
        and not state.navigation.sequence.bounced) {

      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "End of forward filter - switching to smoother mode: "
                << state.navigation.sequence.navigationSurfaces.size()
                << " navigationSurfaces to reverse-navigate.";
        return dstream.str();
      });

      if (!state.navigation.sequence.navigationSurfaces.empty()) {
        debugLog(state, [&] {
          return std::string(
              "Reversing navigation information & step direction.");
        });
        // Unset the targetReached and navigationBreak flags,
        // set reverseMode to true
        state.navigation.targetReached   = false;
        state.navigation.navigationBreak = false;
        state.navigation.reverseMode     = true;
        // The total Accumulated path so far
        double tPath = state.stepping.pathAccumulated;
        // Reverse the path limit
        state.options.pathLimit = -tPath;
        // Set the Stepsize to the maximum path accumulated
        state.stepping.stepSize = detail::ConstrainedStep(-tPath);
        // Set the start parameters to be the new target parameters
        state.navigation.targetSurface = state.navigation.startSurface;
        state.navigation.targetLayer   = state.navigation.startLayer;
        state.navigation.targetVolume  = state.navigation.startVolume;
        // Reverse the direction & (last) step size
        state.stepping.navDir
            = state.stepping.navDir == forward ? backward : forward;
        // Generate a reverse list of surface intersections from
        // the remembered navigationSurfaces
        auto& rSurfaces = state.navigation.sequence.navigationSurfaces;
        std::reverse(rSurfaces.begin(), rSurfaces.end());
        double lStep = 0;
        // Set the path lengths
        std::for_each(rSurfaces.begin(),
                      rSurfaces.end(),
                      [&state, &lStep, tPath](SurfaceIntersection& si) {
                        si.intersection.pathLength -= (tPath);
                        double cStep = si.intersection.pathLength;
                        si.intersection.pathLength -= lStep;
                        lStep         = cStep;
                        si.pDirection = state.stepping.navDir;
                      });
        state.navigation.navSurfaces    = std::move(rSurfaces);
        state.navigation.navSurfaceIter = state.navigation.navSurfaces.begin();
      }
      double step = state.navigation.navSurfaceIter->intersection.pathLength;
      state.stepping.stepSize.update(step, detail::ConstrainedStep::actor);
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
      // Finally remember that you have bounced
      state.navigation.sequence.bounced = true;
      // The sequencer was in action
      return true;
    }

    // Sequence: surface memory - for material and sensitive
    if (memory and state.navigation.currentSurface
        and (state.navigation.currentSurface->associatedMaterial()
             or state.navigation.currentSurface->associatedDetectorElement())
        and not state.navigation.sequence.bounced) {
      // Screen output
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Remember the sequence surface ";
        dstream << state.navigation.currentSurface->geoID().toString();
        return dstream.str();
      });
      state.navigation.sequence.navigationSurfaces.push_back(
          SurfaceIntersection(Intersection(state.stepping.position(),
                                           state.stepping.pathAccumulated,
                                           true),
                              state.navigation.currentSurface,
                              state.stepping.navDir));
      /// Sequencer went into action
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
    }
  }
};

}  // namespace Acts
