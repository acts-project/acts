// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <algorithm>
#include <iterator>
#include <memory>
#include <vector>

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

  DirectNavigator(std::unique_ptr<const Logger> _logger =
                      getDefaultLogger("DirectNavigator", Logging::INFO))
      : m_logger{std::move(_logger)} {}

  /// Nested Actor struct, called Initializer
  ///
  /// This is needed for the initialization of the
  /// surface sequence
  struct Initializer {
    /// The Surface sequence
    SurfaceSequence navSurfaces = {};

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
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param state the entire propagator state
    /// @param r the result of this Actor
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void operator()(propagator_state_t& state, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, result_type& r,
                    const Logger& /*logger*/) const {
      // Only act once
      if (not r.initialized) {
        // Initialize the surface sequence
        state.navigation.navSurfaces = navSurfaces;
        state.navigation.navSurfaceIter = state.navigation.navSurfaces.begin();
        r.initialized = true;
      }
    }
  };

  /// Nested State struct
  ///
  /// It acts as an internal state which is
  /// created for every propagation/extrapolation step
  /// and keep thread-local navigation information
  struct State {
    /// Externally provided surfaces - expected to be ordered
    /// along the path
    SurfaceSequence navSurfaces = {};

    /// Iterator the next surface
    SurfaceIter navSurfaceIter = navSurfaces.begin();

    /// Navigation state - external interface: the start surface
    const Surface* startSurface = nullptr;
    /// Navigation state - external interface: the current surface
    const Surface* currentSurface = nullptr;
    /// Navigation state - external interface: the target surface
    const Surface* targetSurface = nullptr;
    /// Navigation state - starting layer
    const Layer* startLayer = nullptr;
    /// Navigation state - target layer
    const Layer* targetLayer = nullptr;
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

  State makeState(const Surface* startSurface,
                  const Surface* targetSurface) const {
    State result;
    result.startSurface = startSurface;
    result.targetSurface = targetSurface;
    return result;
  }

  /// Reset state
  ///
  /// @param state is the state to reset
  /// @param ssurface is the new starting surface
  /// @param tsurface is the target surface
  void resetState(State& state, const GeometryContext& /*geoContext*/,
                  const Vector3& /*pos*/, const Vector3& /*dir*/,
                  Direction /*navDir*/, const Surface* ssurface,
                  const Surface* tsurface) const {
    // Reset everything except the navSurfaces
    auto navSurfaces = state.navSurfaces;
    state = State();
    state.navSurfaces = navSurfaces;

    // Reset others
    state.navSurfaceIter =
        std::find(state.navSurfaces.begin(), state.navSurfaces.end(), ssurface);
    state.startSurface = ssurface;
    state.currentSurface = ssurface;
    state.targetSurface = tsurface;
  }

  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  const TrackingVolume* currentVolume(const State& state) const {
    return state.currentVolume;
  }

  const IVolumeMaterial* currentVolumeMaterial(const State& state) const {
    if (state.currentVolume == nullptr) {
      return nullptr;
    }
    return state.currentVolume->volumeMaterial();
  }

  const Surface* startSurface(const State& state) const {
    return state.startSurface;
  }

  const Surface* targetSurface(const State& state) const {
    return state.targetSurface;
  }

  bool targetReached(const State& state) const { return state.targetReached; }

  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

  void currentSurface(State& state, const Surface* surface) const {
    state.currentSurface = surface;
  }

  void targetReached(State& state, bool targetReached) const {
    state.targetReached = targetReached;
  }

  void navigationBreak(State& state, bool navigationBreak) const {
    state.navigationBreak = navigationBreak;
  }

  /// @brief Initialize call - start of propagation
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void initialize(propagator_state_t& state, const stepper_t& stepper) const {
    (void)state;
    (void)stepper;
  }

  /// @brief Navigator pre step call
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void preStep(propagator_state_t& state, const stepper_t& stepper) const {
    // Screen output
    ACTS_VERBOSE("Entering navigator::target.");

    // Navigator target always resets the current surface
    state.navigation.currentSurface = nullptr;
    // Output the position in the sequence
    ACTS_VERBOSE(std::distance(state.navigation.navSurfaceIter,
                               state.navigation.navSurfaces.end())
                 << " out of " << state.navigation.navSurfaces.size()
                 << " surfaces remain to try.");

    if (state.navigation.navSurfaceIter != state.navigation.navSurfaces.end()) {
      // Establish & update the surface status
      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, **state.navigation.navSurfaceIter, false);
      if (surfaceStatus == Intersection3D::Status::unreachable) {
        ACTS_VERBOSE(
            "Surface not reachable anymore, switching to next one in "
            "sequence");
        // Move the sequence to the next surface
        ++state.navigation.navSurfaceIter;
      } else {
        ACTS_VERBOSE("Navigation stepSize set to "
                     << stepper.outputStepSize(state.stepping));
      }
    } else {
      // Set the navigation break
      state.navigation.navigationBreak = true;
      // If no externally provided target is given, the target is reached
      if (state.navigation.targetSurface == nullptr) {
        state.navigation.targetReached = true;
        // Announce it then
        ACTS_VERBOSE("No target Surface, job done.");
      }
    }
  }

  /// @brief Navigator post step call
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void postStep(propagator_state_t& state, const stepper_t& stepper) const {
    // Screen output
    ACTS_VERBOSE("Entering navigator::postStep.");

    // Navigator post step always resets the current surface
    state.navigation.currentSurface = nullptr;
    // Output the position in the sequence
    ACTS_VERBOSE(std::distance(state.navigation.navSurfaceIter,
                               state.navigation.navSurfaces.end())
                 << " out of " << state.navigation.navSurfaces.size()
                 << " surfaces remain to try.");

    // Check if we are on surface
    if (state.navigation.navSurfaceIter != state.navigation.navSurfaces.end()) {
      // Establish the surface status
      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, **state.navigation.navSurfaceIter, false);
      if (surfaceStatus == Intersection3D::Status::onSurface) {
        // Set the current surface
        state.navigation.currentSurface = *state.navigation.navSurfaceIter;
        ACTS_VERBOSE("Current surface set to  "
                     << state.navigation.currentSurface->geometryId())
        // Move the sequence to the next surface
        ++state.navigation.navSurfaceIter;
        if (state.navigation.navSurfaceIter !=
            state.navigation.navSurfaces.end()) {
          ACTS_VERBOSE("Next surface candidate is  "
                       << (*state.navigation.navSurfaceIter)->geometryId());
          stepper.releaseStepSize(state.stepping);
        }
      } else if (surfaceStatus == Intersection3D::Status::reachable) {
        ACTS_VERBOSE("Next surface reachable at distance  "
                     << stepper.outputStepSize(state.stepping));
      }
    }
  }

 private:
  const Logger& logger() const { return *m_logger; }

  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts
