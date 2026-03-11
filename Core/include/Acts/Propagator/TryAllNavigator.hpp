// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/NavigatorError.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Propagator/NavigatorStatistics.hpp"
#include "Acts/Propagator/detail/NavigationHelpers.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

namespace Acts::Experimental {

/// @brief Alternative @c Navigator which tries all possible intersections
///
/// See @c Navigator for more general information about the Navigator concept.
///
/// This Navigator tries all possible intersections with all surfaces in the
/// current volume. It does not use any information about the geometry to
/// optimize the search. It is therefore very slow, but can be used as a
/// reference implementation.
///
/// Additionally, this implementation tries to discovers additional
/// intersections after stepping forward and then checking for intersections
/// based on the previous and current positions. This is slower, but more robust
/// against bent tracks.
class TryAllNavigator final {
 public:
  /// Configuration for this Navigator
  struct Config final {
    /// Tracking Geometry for this Navigator
    std::shared_ptr<const TrackingGeometry> trackingGeometry;

    /// stop at every sensitive surface (whether it has material or not)
    bool resolveSensitive = true;
    /// stop at every material surface (whether it is passive or not)
    bool resolveMaterial = true;
    /// stop at every surface regardless what it is
    bool resolvePassive = false;

    /// Which boundary checks to perform for surface approach
    BoundaryTolerance boundaryToleranceSurfaceApproach =
        BoundaryTolerance::None();
  };

  /// Options for this Navigator
  struct Options final : public NavigatorPlainOptions {
    /// @param gctx The geometry context for this navigator instance
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    /// The surface tolerance
    double surfaceTolerance = s_onSurfaceTolerance;

    /// The near limit to resolve surfaces
    double nearLimit = s_onSurfaceTolerance;

    /// The far limit to resolve surfaces
    double farLimit = std::numeric_limits<double>::max();

    /// @param options The plain options to copy
    void setPlainOptions(const NavigatorPlainOptions& options) {
      static_cast<NavigatorPlainOptions&>(*this) = options;
    }
  };

  /// Nested state struct
  struct State final {
    /// @param options_ Navigator options to initialize state with
    explicit State(const Options& options_) : options(options_) {}

    /// Navigation options containing configuration for this propagation
    Options options;

    // Starting geometry information of the navigation which should only be set
    // while initialization. NOTE: This information is mostly used by actors to
    // check if we are on the starting surface (e.g. MaterialInteraction).
    /// Surface where the propagation started
    const Surface* startSurface = nullptr;

    // Target geometry information of the navigation which should only be set
    // while initialization. NOTE: This information is mostly used by actors to
    // check if we are on the target surface (e.g. MaterialInteraction).
    /// Surface that is the target of the propagation
    const Surface* targetSurface = nullptr;

    // Current geometry information of the navigation which is set during
    // initialization and potentially updated after each step.
    /// Currently active surface during propagation
    const Surface* currentSurface = nullptr;
    /// Currently active tracking volume during propagation
    const TrackingVolume* currentVolume = nullptr;

    /// The vector of navigation candidates to work through
    std::vector<detail::NavigationObjectCandidate> navigationCandidates;

    /// If a break has been detected
    bool navigationBreak = false;

    /// Navigation statistics
    NavigatorStatistics statistics;

    /// The vector of active targets ahead of the current position
    std::vector<NavigationTarget> activeTargetsAhead;
    /// The vector of active targets behind the current position
    std::vector<NavigationTarget> activeTargetsBehind;
    /// Index to keep track of which target behind the current position is
    /// currently active
    std::int32_t activeTargetBehindIndex = -1;

    /// The position before the last step
    std::optional<Vector3> lastPosition;

    /// Provides easy access to the current active targets, prioritizing targets
    /// behind the current position.
    /// @return Reference to the vector of currently active intersection
    /// candidates, prioritizing targets behind the current position
    const std::vector<NavigationTarget>& currentTargets() const {
      if (hasTargetsBehind()) {
        return activeTargetsBehind;
      }
      return activeTargetsAhead;
    }

    /// Provides easy access to the active intersection target
    /// @return Reference to the currently active intersection candidate
    const NavigationTarget& activeTargetBehind() const {
      return activeTargetsBehind.at(activeTargetBehindIndex);
    }

    /// Checks if there are still active targets behind the current position
    /// that have not been tried yet
    /// @return True if there are still active targets behind the current position, false otherwise
    bool hasTargetsBehind() const {
      return !activeTargetsBehind.empty() &&
             activeTargetBehindIndex <
                 static_cast<std::int32_t>(activeTargetsBehind.size());
    }
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param logger a logger instance
  explicit TryAllNavigator(Config cfg, std::unique_ptr<const Logger> logger =
                                           getDefaultLogger("TryAllNavigator",
                                                            Logging::INFO))
      : m_cfg(std::move(cfg)), m_logger(std::move(logger)) {}

  /// Creates a new navigator state
  /// @param options Navigator options for state initialization
  /// @return Initialized navigator state with current candidates storage
  State makeState(const Options& options) const {
    State state(options);
    return state;
  }

  /// Get the current surface from the navigation state
  /// @param state The navigation state
  /// @return Pointer to the current surface, or nullptr if none
  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  /// Get the current tracking volume from the navigation state
  /// @param state The navigation state
  /// @return Pointer to the current tracking volume, or nullptr if none
  const TrackingVolume* currentVolume(const State& state) const {
    return state.currentVolume;
  }

  /// Get the material of the current tracking volume
  /// @param state The navigation state
  /// @return Pointer to the volume material, or nullptr if no volume or no material
  const IVolumeMaterial* currentVolumeMaterial(const State& state) const {
    if (state.currentVolume == nullptr) {
      return nullptr;
    }
    return state.currentVolume->volumeMaterial();
  }

  /// Get the start surface from the navigation state
  /// @param state The navigation state
  /// @return Pointer to the start surface, or nullptr if none
  const Surface* startSurface(const State& state) const {
    return state.startSurface;
  }

  /// Get the target surface from the navigation state
  /// @param state The navigation state
  /// @return Pointer to the target surface, or nullptr if none
  const Surface* targetSurface(const State& state) const {
    return state.targetSurface;
  }

  /// Check if the end of the world has been reached
  /// @param state The navigation state
  /// @return True if no current volume is set (end of world reached)
  bool endOfWorldReached(State& state) const {
    return state.currentVolume == nullptr;
  }

  /// Check if navigation has been interrupted
  /// @param state The navigation state
  /// @return True if navigation break flag is set
  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

  /// @brief Initialize the navigator
  ///
  /// This method initializes the navigator for a new propagation. It sets the
  /// current volume and surface to the start volume and surface, respectively.
  ///
  /// @param state The navigation state
  /// @param position The starting position
  /// @param direction The starting direction
  /// @param propagationDirection The propagation direction
  /// @return Result indicating success or failure of initialization
  [[nodiscard]] Result<void> initialize(State& state, const Vector3& position,
                                        const Vector3& direction,
                                        Direction propagationDirection) const {
    static_cast<void>(propagationDirection);

    ACTS_VERBOSE("initialize");

    state.startSurface = state.options.startSurface;
    state.targetSurface = state.options.targetSurface;

    const TrackingVolume* startVolume = nullptr;

    if (state.startSurface != nullptr &&
        state.startSurface->associatedLayer() != nullptr) {
      ACTS_VERBOSE(
          "Fast start initialization through association from Surface.");
      const auto* startLayer = state.startSurface->associatedLayer();
      startVolume = startLayer->trackingVolume();
    } else {
      ACTS_VERBOSE("Slow start initialization through search.");
      ACTS_VERBOSE("Starting from position " << toString(position)
                                             << " and direction "
                                             << toString(direction));
      startVolume = m_cfg.trackingGeometry->lowestTrackingVolume(
          state.options.geoContext, position);
    }

    // Initialize current volume, layer and surface
    {
      state.currentVolume = startVolume;
      if (state.currentVolume != nullptr) {
        ACTS_VERBOSE(volInfo(state) << "Start volume resolved.");
      } else {
        ACTS_DEBUG("Start volume not resolved.");
        state.navigationBreak = true;
        return NavigatorError::NoStartVolume;
      }

      state.currentSurface = state.startSurface;
      if (state.currentSurface != nullptr) {
        ACTS_VERBOSE(volInfo(state) << "Current surface set to start surface "
                                    << state.currentSurface->geometryId());
      } else {
        ACTS_VERBOSE(volInfo(state) << "No start surface set.");
      }
    }

    // Initialize navigation candidates for the start volume
    reinitializeCandidates(state);

    state.lastPosition.reset();

    return Result<void>::success();
  }

  /// @brief Get the next target surface
  ///
  /// This method gets the next target surface based on the current
  /// position and direction. It returns a none target if no target can be
  /// found.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return The next target surface
  NavigationTarget nextTarget(State& state, const Vector3& position,
                              const Vector3& direction) const {
    // Navigator preStep always resets the current surface
    state.currentSurface = nullptr;

    // Check if the navigator is inactive
    if (state.navigationBreak) {
      return NavigationTarget::None();
    }

    ACTS_VERBOSE(volInfo(state) << "nextTarget");

    if (state.lastPosition.has_value() && !state.hasTargetsBehind()) {
      ACTS_VERBOSE(volInfo(state) << "Evaluate blind step");

      const Vector3 stepStart = state.lastPosition.value();
      state.lastPosition.reset();
      const Vector3 stepEnd = position;
      const Vector3 step = stepEnd - stepStart;
      const double stepDistance = step.norm();

      ACTS_VERBOSE("- from: " << stepStart.transpose());
      ACTS_VERBOSE("- to: " << stepEnd.transpose());
      ACTS_VERBOSE("- distance: " << stepDistance);

      if (stepDistance < std::numeric_limits<double>::epsilon()) {
        ACTS_DEBUG(volInfo(state)
                   << "Step distance is zero: " << stepDistance
                   << ". Try another to resolve the target again.");
        return nextTarget(state, position, direction);
      }

      const Vector3 stepDirection = step.normalized();

      const double nearLimit = -stepDistance + state.options.surfaceTolerance;
      const double farLimit = 0;

      state.activeTargetsBehind =
          resolveTargets(state, stepEnd, stepDirection, nearLimit, farLimit);
      state.activeTargetBehindIndex = -1;

      ACTS_VERBOSE(volInfo(state)
                   << "Found " << state.activeTargetsBehind.size()
                   << " intersections behind");

      for (const auto& target : state.activeTargetsBehind) {
        ACTS_VERBOSE("Found target behind " << target.surface().geometryId());
      }
    }

    // Prioritize targets behind the current position
    ++state.activeTargetBehindIndex;
    if (state.hasTargetsBehind()) {
      ACTS_VERBOSE(volInfo(state) << "Handle active candidates behind");

      ACTS_VERBOSE(volInfo(state)
                   << (state.activeTargetsBehind.size() -
                       state.activeTargetBehindIndex)
                   << " out of " << state.activeTargetsBehind.size()
                   << " surfaces remain to try.");

      const NavigationTarget& nextTarget = state.activeTargetBehind();

      ACTS_VERBOSE(volInfo(state) << "Next target behind selected: "
                                  << nextTarget.surface().geometryId());

      return nextTarget;
    }

    // No more targets behind, now try to find targets ahead as usual

    state.lastPosition = position;
    state.activeTargetsBehind.clear();
    state.activeTargetBehindIndex = -1;

    ACTS_VERBOSE(volInfo(state)
                 << "No targets behind, try to find targets ahead");

    const double nearLimit = state.options.nearLimit;
    const double farLimit = state.options.farLimit;

    state.activeTargetsAhead =
        resolveTargets(state, position, direction, nearLimit, farLimit);

    NavigationTarget nextTarget = NavigationTarget::None();

    for (const auto& target : state.activeTargetsAhead) {
      const Intersection3D& intersection = target.intersection();

      if (intersection.status() == IntersectionStatus::onSurface) {
        ACTS_ERROR(volInfo(state)
                   << "We are on surface " << target.surface().geometryId()
                   << " before trying to reach it. This should not happen. "
                      "Good luck.");
        continue;
      }

      if (intersection.status() == IntersectionStatus::reachable) {
        nextTarget = target;
        break;
      }
    }

    if (nextTarget.isNone()) {
      ACTS_VERBOSE(volInfo(state)
                   << "No target ahead found. Step blindly forward.");
    } else {
      ACTS_VERBOSE(volInfo(state) << "Next target ahead selected: "
                                  << nextTarget.surface().geometryId());
    }

    return nextTarget;
  }

  /// @brief Check if the target is still valid
  ///
  /// This method checks if the target is valid based on the current position
  /// and direction. It returns true if the target is still valid.
  ///
  /// For the TryAllNavigator, the target is always invalid since we do not want
  /// to assume any specific surface sequence over multiple steps.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return True if the target is still valid
  bool checkTargetValid(const State& state, const Vector3& position,
                        const Vector3& direction) const {
    static_cast<void>(state);
    static_cast<void>(position);
    static_cast<void>(direction);

    return false;
  }

  /// @brief Handle the surface reached
  ///
  /// This method is called when a surface is reached. It sets the current
  /// surface in the navigation state and updates the navigation candidates.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void handleSurfaceReached(State& state, const Vector3& position,
                            const Vector3& direction,
                            const Surface& /*surface*/) const {
    // Check if the navigator is inactive
    if (state.navigationBreak) {
      return;
    }

    ACTS_VERBOSE(volInfo(state) << "handleSurfaceReached");

    const std::vector<NavigationTarget>& currentTargets =
        state.currentTargets();

    if (currentTargets.empty()) {
      ACTS_VERBOSE(volInfo(state) << "No current target set.");
      return;
    }

    assert(state.currentSurface == nullptr && "Current surface must be reset.");

    // handle multiple surface intersections due to increased bounds

    std::vector<NavigationTarget> hitTargets;

    for (const auto& target : currentTargets) {
      const std::uint8_t index = target.intersectionIndex();
      const Surface& surface = target.surface();
      const BoundaryTolerance boundaryTolerance = BoundaryTolerance::None();

      const Intersection3D intersection =
          surface
              .intersect(state.options.geoContext, position, direction,
                         boundaryTolerance, state.options.surfaceTolerance)
              .at(index);

      if (intersection.status() == IntersectionStatus::onSurface) {
        hitTargets.emplace_back(target);
      }
    }

    ACTS_VERBOSE(volInfo(state)
                 << "Found " << hitTargets.size()
                 << " intersections on surface with bounds check.");

    // reset stored targets
    state.lastPosition.reset();
    state.activeTargetsAhead.clear();
    state.activeTargetsBehind.clear();
    state.activeTargetBehindIndex = -1;

    if (hitTargets.empty()) {
      ACTS_VERBOSE(volInfo(state) << "No hit targets found.");
      return;
    }

    if (hitTargets.size() > 1) {
      ACTS_VERBOSE(volInfo(state)
                   << "Only using first intersection within bounds.");
    }

    // we can only handle a single surface hit so we pick the first one
    const NavigationTarget& target = hitTargets.front();
    const Surface& surface = target.surface();

    ACTS_VERBOSE(volInfo(state) << "Surface " << surface.geometryId()
                                << " successfully hit, storing it.");
    state.currentSurface = &surface;

    if (target.isSurfaceTarget()) {
      ACTS_VERBOSE(volInfo(state) << "This is a surface");
    } else if (target.isLayerTarget()) {
      ACTS_VERBOSE(volInfo(state) << "This is a layer");
    } else if (target.isPortalTarget()) {
      ACTS_VERBOSE(volInfo(state)
                   << "This is a boundary. Reinitialize navigation");

      const BoundarySurface& boundary = target.boundarySurface();

      state.currentVolume = boundary.attachedVolume(state.options.geoContext,
                                                    position, direction);

      ACTS_VERBOSE(volInfo(state) << "Switched volume");

      reinitializeCandidates(state);
    } else {
      ACTS_ERROR(volInfo(state) << "Unknown intersection type");
    }
  }

 private:
  /// Configuration object for this navigator
  Config m_cfg;

  /// Logger instance for this navigator
  std::unique_ptr<const Logger> m_logger;

  /// @brief Get the logger instance
  /// @return Reference to the logger instance
  const Logger& logger() const { return *m_logger; }

  /// Helper method to reset and reinitialize the navigation candidates.
  void reinitializeCandidates(State& state) const {
    state.navigationCandidates.clear();
    state.activeTargetsAhead.clear();
    state.activeTargetsBehind.clear();
    state.activeTargetBehindIndex = -1;

    initializeVolumeCandidates(state);
  }

  /// Helper method to initialize navigation candidates for the current volume.
  /// @param state Navigation state to initialize candidates for
  void initializeVolumeCandidates(State& state) const {
    const TrackingVolume* volume = state.currentVolume;
    ACTS_VERBOSE(volInfo(state) << "Initialize volume");

    if (volume == nullptr) {
      state.navigationBreak = true;
      ACTS_VERBOSE(volInfo(state) << "No volume set. Good luck.");
      return;
    }

    emplaceAllVolumeCandidates(
        state.navigationCandidates, *volume, m_cfg.resolveSensitive,
        m_cfg.resolveMaterial, m_cfg.resolvePassive,
        m_cfg.boundaryToleranceSurfaceApproach, logger());
  }

  std::vector<NavigationTarget> resolveTargets(State& state,
                                               const Vector3& position,
                                               const Vector3& direction,
                                               double nearLimit,
                                               double farLimit) const {
    std::vector<NavigationTarget> targets;

    // Find intersections with all candidates
    for (const auto& candidate : state.navigationCandidates) {
      auto intersections =
          candidate.intersect(state.options.geoContext, position, direction,
                              state.options.surfaceTolerance);
      for (auto [intersectionIndex, intersection] :
           Acts::enumerate(intersections)) {
        // exclude invalid intersections
        if (!intersection.isValid() ||
            !detail::checkPathLength(intersection.pathLength(), nearLimit,
                                     farLimit)) {
          continue;
        }
        // store candidate
        targets.emplace_back(candidate.target(intersection, intersectionIndex));
      }
    }

    std::ranges::sort(targets, NavigationTarget::pathLengthOrder);

    return targets;
  }

  /// @brief Get volume information string for logging
  /// @param state The navigation state
  /// @return String containing volume name or "No Volume" followed by separator
  std::string volInfo(const State& state) const {
    return (state.currentVolume != nullptr ? state.currentVolume->volumeName()
                                           : "No Volume") +
           " | ";
  }
};

}  // namespace Acts::Experimental
