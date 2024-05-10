// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <memory>
#include <vector>

namespace Acts {

/// This is a fully guided navigator that progresses through a pre-given
/// sequence of surfaces.
///
/// This can either be used as a validation tool, for truth tracking, or track
/// refitting.
class DirectNavigator {
 public:
  /// The sequentially crossed surfaces
  using SurfaceSequence = std::vector<const Surface*>;
  using SurfaceIter = std::vector<const Surface*>::iterator;

  struct Config {};

  struct Options : public NavigatorPlainOptions {
    void setPlainOptions(const NavigatorPlainOptions& options) {
      static_cast<NavigatorPlainOptions&>(*this) = options;
    }
  };

  /// @brief Nested Actor struct, called Initializer
  ///
  /// This is needed for the initialization of the surface sequence.
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
      if (!r.initialized) {
        // Initialize the surface sequence
        state.navigation.navSurfaces = navSurfaces;
        state.navigation.navSurfaceIter = state.navigation.navSurfaces.begin();

        // In case the start surface is in the list of nav surfaces
        // we need to correct the iterator to point to the next surface
        // in the vector
        if (state.navigation.options.startSurface) {
          auto surfaceIter = std::find(state.navigation.navSurfaces.begin(),
                                       state.navigation.navSurfaces.end(),
                                       state.navigation.options.startSurface);
          // if current surface in the list, point to the next surface
          if (surfaceIter != state.navigation.navSurfaces.end()) {
            state.navigation.navSurfaceIter = ++surfaceIter;
          }
        }

        r.initialized = true;
      }
    }
  };

  /// @brief Nested State struct
  ///
  /// It acts as an internal state which is created for every
  /// propagation/extrapolation step and keep thread-local navigation
  /// information
  struct State {
    Options options;

    /// Externally provided surfaces - expected to be ordered along the path
    SurfaceSequence navSurfaces = {};

    /// Iterator the next surface
    SurfaceIter navSurfaceIter = navSurfaces.begin();

    /// Navigation state - external interface: the current surface
    const Surface* currentSurface = nullptr;

    /// Navigation state - external interface: target is reached
    bool targetReached = false;
    /// Navigation state - external interface: a break has been detected
    bool navigationBreak = false;
  };

  DirectNavigator(std::unique_ptr<const Logger> _logger =
                      getDefaultLogger("DirectNavigator", Logging::INFO))
      : m_logger{std::move(_logger)} {}

  State makeState(const Options& options) const {
    State state;
    state.options = options;
    return state;
  }

  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  const TrackingVolume* currentVolume(const State& /*state*/) const {
    return nullptr;
  }

  const IVolumeMaterial* currentVolumeMaterial(const State& /*state*/) const {
    return nullptr;
  }

  const Surface* startSurface(const State& state) const {
    return state.options.startSurface;
  }

  const Surface* targetSurface(const State& state) const {
    return state.options.targetSurface;
  }

  bool targetReached(const State& state) const { return state.targetReached; }

  bool endOfWorldReached(State& /*state*/) const { return false; }

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
  /// @tparam propagator_state_t The state type of the propagator
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  template <typename propagator_state_t, typename stepper_t>
  void initialize(propagator_state_t& state,
                  const stepper_t& /*stepper*/) const {
    ACTS_VERBOSE("initialize");

    // We set the current surface to the start surface
    state.navigation.currentSurface = state.navigation.options.startSurface;
    if (state.navigation.currentSurface) {
      ACTS_VERBOSE("Current surface set to start surface "
                   << state.navigation.currentSurface->geometryId());
    }
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
    ACTS_VERBOSE("pre step");

    // Navigator target always resets the current surface
    state.navigation.currentSurface = nullptr;
    // Output the position in the sequence
    ACTS_VERBOSE(std::distance(state.navigation.navSurfaceIter,
                               state.navigation.navSurfaces.end())
                 << " out of " << state.navigation.navSurfaces.size()
                 << " surfaces remain to try.");

    if (state.navigation.navSurfaceIter != state.navigation.navSurfaces.end()) {
      // Establish & update the surface status
      // TODO we do not know the intersection index - passing the closer one
      const auto& surface = **state.navigation.navSurfaceIter;
      const double farLimit = std::numeric_limits<double>::max();
      const auto index =
          chooseIntersection(
              state.geoContext, surface, stepper.position(state.stepping),
              state.options.direction * stepper.direction(state.stepping),
              BoundaryCheck(false), m_nearLimit, farLimit,
              state.options.surfaceTolerance)
              .index();
      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, surface, index, state.options.direction,
          BoundaryCheck(false), state.options.surfaceTolerance, *m_logger);
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
      if (state.navigation.options.targetSurface == nullptr) {
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
    ACTS_VERBOSE("post step");

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
      // TODO we do not know the intersection index - passing the closer one
      const auto& surface = **state.navigation.navSurfaceIter;
      const double farLimit = std::numeric_limits<double>::max();
      const auto index =
          chooseIntersection(
              state.geoContext, surface, stepper.position(state.stepping),
              state.options.direction * stepper.direction(state.stepping),
              BoundaryCheck(false), m_nearLimit, farLimit,
              state.options.surfaceTolerance)
              .index();
      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, surface, index, state.options.direction,
          BoundaryCheck(false), state.options.surfaceTolerance, *m_logger);
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
        }
      } else if (surfaceStatus == Intersection3D::Status::reachable) {
        ACTS_VERBOSE("Next surface reachable at distance  "
                     << stepper.outputStepSize(state.stepping));
      }
    }
  }

 private:
  ObjectIntersection<Surface> chooseIntersection(
      const GeometryContext& gctx, const Surface& surface,
      const Vector3& position, const Vector3& direction,
      const BoundaryCheck& bcheck, double nearLimit, double farLimit,
      double tolerance) const {
    auto intersections =
        surface.intersect(gctx, position, direction, bcheck, tolerance);

    for (auto& intersection : intersections.split()) {
      if (detail::checkIntersection(intersection, nearLimit, farLimit,
                                    logger())) {
        return intersection;
      }
    }

    return ObjectIntersection<Surface>::invalid();
  }

  const Logger& logger() const { return *m_logger; }

  std::unique_ptr<const Logger> m_logger;

  // TODO https://github.com/acts-project/acts/issues/2738
  /// Distance limit to discard intersections "behind us"
  /// @note this is only necessary because some surfaces have more than one
  ///       intersection
  double m_nearLimit = -100 * UnitConstants::um;
};

}  // namespace Acts
