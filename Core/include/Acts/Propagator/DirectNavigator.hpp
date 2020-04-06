// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/algorithm/string.hpp>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <string>
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

/// DirectNavigator class
///
/// This is a fully guided navigator that progresses through
/// a pre-given sequence of surfaces.
///
/// This can either be used as a validation tool, for truth
/// tracking, or track refitting
class DirectNavigator {
 public:
  /// The sequentially crossed surfaces
  using SurfaceSequence = std::vector<const Surface*>;
  using SurfaceIter = std::vector<const Surface*>::iterator;

  /// Defaulted Constructed
  DirectNavigator() = default;

  /// The tolerance used to define "surface reached"
  double tolerance = s_onSurfaceTolerance;

  /// Nested Actor struct, called Initializer
  ///
  /// This is needed for the initialization of the
  /// surface sequence
  struct Initializer {
    /// The Surface sequence
    SurfaceSequence surfaceSequence = {};

    /// Actor result / state
    struct this_result {
      bool initialized = false;
    };
    using result_type = this_result;

    /// Defaulting the constructor
    Initializer() = default;

    /// Actor operator call
    /// @tparam statet Type of the full propagator state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state the entire propagator state
    /// @param r the result of this Actor
    template <typename propagator_state_t, typename stepper_t>
    void operator()(propagator_state_t& state, const stepper_t& /*unused*/,
                    result_type& r) const {
      // Only act once
      if (not r.initialized) {
        // Initialize the surface sequence
        state.navigation.surfaceSequence = surfaceSequence;
        state.navigation.nextSurfaceIter =
            state.navigation.surfaceSequence.begin();
        r.initialized = true;
      }
    }

    /// Actor operator call - resultless, unused
    template <typename propagator_state_t, typename stepper_t>
    void operator()(propagator_state_t& /*unused*/,
                    const stepper_t& /*unused*/) const {}
  };

  /// Nested State struct
  ///
  /// It acts as an internal state which is
  /// created for every propagation/extrapolation step
  /// and keep thread-local navigation information
  struct State {
    /// Externally provided surfaces - expected to be ordered
    /// along the path
    SurfaceSequence surfaceSequence = {};

    /// Iterator the the next surface
    SurfaceIter nextSurfaceIter = surfaceSequence.begin();

    /// Navigation state - external interface: the start surface
    const Surface* startSurface = nullptr;
    /// Navigation state - external interface: the current surface
    const Surface* currentSurface = nullptr;
    /// Navigation state - external interface: the target surface
    const Surface* targetSurface = nullptr;
    /// Navigation state: the start volume
    const TrackingVolume* startVolume = nullptr;
    /// Navigation state: the current volume
    const TrackingVolume* currentVolume = nullptr;
    /// Navigation state: the target volume
    const TrackingVolume* targetVolume = nullptr;

    /// Navigation state - external interface: target is reached
    bool targetReached = false;
    /// Navigation state - external interface: a break has been detected
    bool navigationBreak = false;
  };

  /// @brief Navigator status call
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void status(propagator_state_t& state, const stepper_t& stepper) const {
    // Screen output
    debugLog(state, [&] { return std::string("Entering navigator::status."); });

    // Navigator status always resets the current surface
    state.navigation.currentSurface = nullptr;
    // Output the position in the sequence
    debugLog(state, [&] {
      std::stringstream dstream;
      dstream << std::distance(state.navigation.nextSurfaceIter,
                               state.navigation.surfaceSequence.end());
      dstream << " out of " << state.navigation.surfaceSequence.size();
      dstream << " surfaces remain to try.";
      return dstream.str();
    });
    // Check if we are on surface
    if (state.navigation.nextSurfaceIter !=
        state.navigation.surfaceSequence.end()) {
      // Establish the surface status
      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, **state.navigation.nextSurfaceIter, false);
      if (surfaceStatus == Intersection::Status::onSurface) {
        // Set the current surface
        state.navigation.currentSurface = *state.navigation.nextSurfaceIter;
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << "Current surface set to  "
                  << state.navigation.currentSurface->geoID();
          return dstream.str();
        });
        // Move the sequence to the next surface
        ++state.navigation.nextSurfaceIter;
        if (state.navigation.nextSurfaceIter !=
            state.navigation.surfaceSequence.end()) {
          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << "Next surface candidate is  "
                    << (*state.navigation.nextSurfaceIter)->geoID();
            return dstream.str();
          });
          stepper.releaseStepSize(state.stepping);
        }
      } else if (surfaceStatus == Intersection::Status::reachable) {
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << "Next surface reachable at distance  "
                  << stepper.outputStepSize(state.stepping);
          return dstream.str();
        });
      }
    }
  }

  /// @brief Navigator target call
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void target(propagator_state_t& state, const stepper_t& stepper) const {
    // Screen output
    debugLog(state, [&] { return std::string("Entering navigator::target."); });

    // Navigator target always resets the current surface
    state.navigation.currentSurface = nullptr;
    // Output the position in the sequence
    debugLog(state, [&] {
      std::stringstream dstream;
      dstream << std::distance(state.navigation.nextSurfaceIter,
                               state.navigation.surfaceSequence.end());
      dstream << " out of " << state.navigation.surfaceSequence.size();
      dstream << " surfaces remain to try.";
      return dstream.str();
    });
    if (state.navigation.nextSurfaceIter !=
        state.navigation.surfaceSequence.end()) {
      // Establish & update the surface status
      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, **state.navigation.nextSurfaceIter, false);
      if (surfaceStatus == Intersection::Status::unreachable) {
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << "Surface not reachable anymore, switching to next one in "
                     "sequence";
          return dstream.str();
        });
        // Move the sequence to the next surface
        ++state.navigation.nextSurfaceIter;
      } else {
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << "Navigation stepSize set to ";
          dstream << stepper.outputStepSize(state.stepping);
          return dstream.str();
        });
      }
    } else {
      // Set the navigation break
      state.navigation.navigationBreak = true;
      // If no externally provided target is given, the target is reached
      if (state.navigation.targetSurface == nullptr) {
        state.navigation.targetReached = true;
        // Announce it then
        debugLog(state,
                 [&] { return std::string("No target Surface, job done."); });
      }
    }
  }

 private:
  /// The private navigation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the
  /// state.options.debug == true case in order not to spend time
  /// when not needed.
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param[in,out] state the propagator state for the debug flag,
  ///      prefix and length
  /// @param logAction is a callable function that returns a streamable object
  template <typename propagator_state_t>
  void debugLog(propagator_state_t& state,
                const std::function<std::string()>& logAction) const {
    if (state.options.debug) {
      std::string vName = "Direct Navigator";
      std::vector<std::string> lines;
      std::string input = logAction();
      boost::split(lines, input, boost::is_any_of("\n"));
      for (const auto& line : lines) {
        std::stringstream dstream;
        dstream << ">>>" << std::setw(state.options.debugPfxWidth) << vName
                << " | ";
        dstream << std::setw(state.options.debugMsgWidth) << line << '\n';
        state.options.debugString += dstream.str();
      }
    }
  }
};

}  // namespace Acts
