// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Propagator/NavigatorStatistics.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <limits>
#include <memory>
#include <vector>

namespace Acts {

/// @brief A fully guided navigator
///
/// This is a fully guided navigator that progresses through a provided sequence
/// of surfaces.
///
/// This can either be used as a validation tool, for truth tracking, or track
/// refitting.
///
class DirectNavigator {
 public:
  /// The sequentially crossed surfaces
  using SurfaceSequence = std::vector<const Surface*>;

  /// @brief The nested configuration struct
  struct Config {};

  /// @brief The nested options struct
  struct Options : public NavigatorPlainOptions {
    /// Constructor from geometry context
    /// @param gctx The geometry context
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    /// The Surface sequence
    SurfaceSequence surfaces;

    /// The surface tolerance
    double surfaceTolerance = s_onSurfaceTolerance;

    // TODO https://github.com/acts-project/acts/issues/2738
    /// Distance limit to discard intersections "behind us"
    /// @note this is only necessary because some surfaces have more than one
    ///       intersection
    double nearLimit = -100 * UnitConstants::um;

    /// The far limit to resolve surfaces
    double farLimit = std::numeric_limits<double>::max();

    /// Set the plain navigator options
    /// @param options The plain navigator options to copy
    void setPlainOptions(const NavigatorPlainOptions& options) {
      static_cast<NavigatorPlainOptions&>(*this) = options;
    }
  };

  /// @brief Nested State struct
  ///
  /// It acts as an internal state which is created for every
  /// propagation/extrapolation step and keep thread-local navigation
  /// information
  struct State {
    /// Constructor from options
    /// @param options_ The navigator options
    explicit State(const Options& options_) : options(options_) {}

    /// Configuration options for the direct navigator
    Options options;

    /// Propagation direction (forward or backward)
    Direction direction = Direction::Forward();

    /// Index of the next surface to try
    /// @note -1 means before the first surface in the sequence and size()
    ///       means after the last surface in the sequence
    int surfaceIndex = -1;

    /// Navigation state - external interface: the current surface
    const Surface* currentSurface = nullptr;

    /// Navigation state - external interface: a break has been detected
    bool navigationBreak = false;

    /// Navigation statistics
    NavigatorStatistics statistics;

    /// Get the current navigation surface
    /// @return Reference to the surface at the current surface index
    const Surface& navSurface() const {
      return *options.surfaces.at(surfaceIndex);
    }

    /// Move to the next surface in the sequence
    /// Increments or decrements surface index based on propagation direction
    void nextSurface() {
      if (direction == Direction::Forward()) {
        ++surfaceIndex;
      } else {
        --surfaceIndex;
      }
    }

    /// Check if we have reached the end of the surface sequence
    /// @return True if no more surfaces remain in the propagation direction
    bool endOfSurfaces() const {
      if (direction == Direction::Forward()) {
        return surfaceIndex >= static_cast<int>(options.surfaces.size());
      }
      return surfaceIndex < 0;
    }

    /// Get the number of surfaces remaining in the sequence
    /// @return Number of surfaces left to process in the propagation direction
    int remainingSurfaces() const {
      if (direction == Direction::Forward()) {
        return options.surfaces.size() - surfaceIndex;
      }
      return surfaceIndex + 1;
    }

    /// Reset the surface index to the initial position
    /// Sets index to before first surface (forward) or after last surface
    /// (backward)
    void resetSurfaceIndex() {
      surfaceIndex = direction == Direction::Forward()
                         ? -1
                         : static_cast<int>(options.surfaces.size());
    }
  };

  /// Constructor with optional logger
  /// @param _logger Logger instance for navigation messages
  explicit DirectNavigator(std::unique_ptr<const Logger> _logger =
                               getDefaultLogger("DirectNavigator",
                                                Logging::INFO))
      : m_logger{std::move(_logger)} {}

  /// Create a new navigation state from options
  /// @param options The navigator options
  /// @return New state initialized with the provided options
  State makeState(const Options& options) const {
    State state(options);
    return state;
  }

  /// Get the current surface from the navigation state
  /// @param state The navigation state
  /// @return Pointer to current surface, nullptr if none
  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  /// Get the current tracking volume (not used by DirectNavigator)
  /// @param state The navigation state (unused)
  /// @return Always returns nullptr as DirectNavigator doesn't use volumes
  const TrackingVolume* currentVolume(const State& state) const {
    static_cast<void>(state);
    return nullptr;
  }

  /// Get the current volume material (not used by DirectNavigator)
  /// @param state The navigation state (unused)
  /// @return Always returns nullptr as DirectNavigator doesn't use volume material
  const IVolumeMaterial* currentVolumeMaterial(const State& state) const {
    static_cast<void>(state);
    return nullptr;
  }

  /// Get the start surface from the navigation state
  /// @param state The navigation state
  /// @return Pointer to start surface, nullptr if none
  const Surface* startSurface(const State& state) const {
    return state.options.startSurface;
  }

  /// Get the target surface from the navigation state
  /// @param state The navigation state
  /// @return Pointer to target surface, nullptr if none
  const Surface* targetSurface(const State& state) const {
    return state.options.targetSurface;
  }

  /// Check if the end of world has been reached (not applicable for
  /// DirectNavigator)
  /// @param state The navigation state (unused)
  /// @return Always returns false as DirectNavigator operates on a surface sequence
  bool endOfWorldReached(State& state) const {
    static_cast<void>(state);
    return false;
  }

  /// Check if navigation should break/stop
  /// @param state The navigation state
  /// @return True if navigation should break, false otherwise
  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

  /// @brief Initialize the navigator
  ///
  /// This function initializes the navigator for a new propagation.
  ///
  /// @param state The navigation state
  /// @param position The start position
  /// @param direction The start direction
  /// @param propagationDirection The propagation direction
  /// @return Always returns success result as DirectNavigator initialization cannot fail
  [[nodiscard]] Result<void> initialize(State& state, const Vector3& position,
                                        const Vector3& direction,
                                        Direction propagationDirection) const {
    static_cast<void>(position);
    static_cast<void>(direction);

    ACTS_VERBOSE("Initialize. Surface sequence for navigation:");
    for (const Surface* surface : state.options.surfaces) {
      ACTS_VERBOSE(surface->geometryId()
                   << " - "
                   << surface->center(state.options.geoContext).transpose());
    }

    state.direction = propagationDirection;
    ACTS_VERBOSE("Navigation direction is " << propagationDirection);

    // We set the current surface to the start surface
    state.currentSurface = state.options.startSurface;
    if (state.currentSurface != nullptr) {
      ACTS_VERBOSE("Current surface set to start surface "
                   << state.currentSurface->geometryId());
    } else {
      ACTS_VERBOSE("Current surface set to nullptr");
    }

    // Find initial index.
    auto found =
        std::ranges::find(state.options.surfaces, state.options.startSurface);

    if (found != state.options.surfaces.end()) {
      // The index should be the index before the start surface, depending on
      // the direction
      state.surfaceIndex = std::distance(state.options.surfaces.begin(), found);
      state.surfaceIndex += state.direction == Direction::Backward() ? 1 : -1;
    } else {
      ACTS_DEBUG(
          "Did not find the start surface in the sequence. Assuming it is not "
          "part of the sequence. Trusting the correctness of the input "
          "sequence. Resetting the surface index.");
      state.resetSurfaceIndex();
    }

    state.navigationBreak = false;

    return Result<void>::success();
  }

  /// @brief Get the next target surface
  ///
  /// This function gets the next target surface for the propagation. For
  /// the direct navigator this is always the next surface in the sequence.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return The next target surface
  NavigationTarget nextTarget(State& state, const Vector3& position,
                              const Vector3& direction) const {
    // Navigator target always resets the current surface
    state.currentSurface = nullptr;

    if (state.navigationBreak) {
      return NavigationTarget::None();
    }

    ACTS_VERBOSE("DirectNavigator::nextTarget");

    // Move the sequence to the next surface
    state.nextSurface();

    while (!state.endOfSurfaces()) {
      ACTS_VERBOSE("Next surface candidate is "
                   << state.navSurface().geometryId() << ". "
                   << state.remainingSurfaces() << " out of "
                   << state.options.surfaces.size()
                   << " surfaces remain to try.");

      // Establish & update the surface status
      // TODO we do not know the intersection index - passing the closer one
      const Surface& surface = state.navSurface();
      const double farLimit = std::numeric_limits<double>::max();
      const NavigationTarget target = chooseIntersection(
          state.options.geoContext, surface, position, direction,
          BoundaryTolerance::Infinite(), state.options.nearLimit, farLimit,
          state.options.surfaceTolerance);
      if (target.isValid()) {
        return target;
      }

      ACTS_VERBOSE("No valid intersection found with surface "
                   << surface.geometryId() << ", trying next surface.");
      state.nextSurface();
    }

    ACTS_VERBOSE("End of surfaces reached, navigation break.");
    state.navigationBreak = true;
    return NavigationTarget::None();
  }

  /// @brief Check if the current target is still valid
  ///
  /// This function checks if the target is valid. For the direct navigator this
  /// is always true.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return True if the target is valid
  bool checkTargetValid(const State& state, const Vector3& position,
                        const Vector3& direction) const {
    static_cast<void>(state);
    static_cast<void>(position);
    static_cast<void>(direction);

    return true;
  }

  /// @brief Handle the surface reached
  ///
  /// This function handles the surface reached. For the direct navigator this
  /// effectively sets the current surface to the reached surface.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  /// @param surface The surface reached
  void handleSurfaceReached(State& state, const Vector3& position,
                            const Vector3& direction,
                            const Surface& surface) const {
    static_cast<void>(position);
    static_cast<void>(direction);
    static_cast<void>(surface);

    if (state.navigationBreak) {
      return;
    }

    ACTS_VERBOSE("DirectNavigator::handleSurfaceReached");

    // Set the current surface
    state.currentSurface = &state.navSurface();
    ACTS_VERBOSE("Current surface set to "
                 << state.currentSurface->geometryId());
  }

 private:
  NavigationTarget chooseIntersection(
      const GeometryContext& gctx, const Surface& surface,
      const Vector3& position, const Vector3& direction,
      const BoundaryTolerance& boundaryTolerance, double nearLimit,
      double farLimit, double tolerance) const {
    auto intersections = surface.intersect(gctx, position, direction,
                                           boundaryTolerance, tolerance);

    for (auto [intersectionIndex, intersection] :
         Acts::enumerate(intersections)) {
      if (intersection.isValid() &&
          detail::checkPathLength(intersection.pathLength(), nearLimit,
                                  farLimit, logger())) {
        return NavigationTarget(intersection, intersectionIndex, surface,
                                boundaryTolerance);
      }
    }

    return NavigationTarget::None();
  }

  const Logger& logger() const { return *m_logger; }

  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts
