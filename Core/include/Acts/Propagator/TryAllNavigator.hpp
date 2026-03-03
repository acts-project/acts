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

namespace Acts {

/// @brief Captures the common functionality of the try-all navigators
///
/// This class is not meant to be used directly, but to be inherited by the
/// actual navigator implementations.
///
class TryAllNavigatorBase {
 public:
  /// @brief Configuration for this Navigator
  struct Config {
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

  /// @brief Options for this Navigator
  struct Options : public NavigatorPlainOptions {
    /// @brief Constructor with geometry context
    /// @param gctx The geometry context for this navigator instance
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    /// The surface tolerance
    double surfaceTolerance = s_onSurfaceTolerance;

    /// The near limit to resolve surfaces
    double nearLimit = s_onSurfaceTolerance;

    /// The far limit to resolve surfaces
    double farLimit = std::numeric_limits<double>::max();

    /// @brief Set plain options from NavigatorPlainOptions
    /// @param options The plain options to copy
    void setPlainOptions(const NavigatorPlainOptions& options) {
      static_cast<NavigatorPlainOptions&>(*this) = options;
    }
  };

  /// @brief Nested State struct
  ///
  /// It acts as an internal state which is created for every propagation and
  /// meant to keep thread-local navigation information.
  struct State {
    /// @brief Constructor with options
    /// @param options_ The navigator options for this state
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
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param logger a logger instance
  TryAllNavigatorBase(Config cfg, std::unique_ptr<const Logger> logger)
      : m_cfg(std::move(cfg)), m_logger{std::move(logger)} {}

  /// @brief Get the current surface from the navigation state
  /// @param state The navigation state
  /// @return Pointer to the current surface, or nullptr if none
  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  /// @brief Get the current tracking volume from the navigation state
  /// @param state The navigation state
  /// @return Pointer to the current tracking volume, or nullptr if none
  const TrackingVolume* currentVolume(const State& state) const {
    return state.currentVolume;
  }

  /// @brief Get the material of the current tracking volume
  /// @param state The navigation state
  /// @return Pointer to the volume material, or nullptr if no volume or no material
  const IVolumeMaterial* currentVolumeMaterial(const State& state) const {
    if (state.currentVolume == nullptr) {
      return nullptr;
    }
    return state.currentVolume->volumeMaterial();
  }

  /// @brief Get the start surface from the navigation state
  /// @param state The navigation state
  /// @return Pointer to the start surface, or nullptr if none
  const Surface* startSurface(const State& state) const {
    return state.startSurface;
  }

  /// @brief Get the target surface from the navigation state
  /// @param state The navigation state
  /// @return Pointer to the target surface, or nullptr if none
  const Surface* targetSurface(const State& state) const {
    return state.targetSurface;
  }

  /// @brief Check if the end of the world has been reached
  /// @param state The navigation state
  /// @return True if no current volume is set (end of world reached)
  bool endOfWorldReached(State& state) const {
    return state.currentVolume == nullptr;
  }

  /// @brief Check if navigation has been interrupted
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

    return Result<void>::success();
  }

 protected:
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

  /// @brief Get volume information string for logging
  /// @param state The navigation state
  /// @return String containing volume name or "No Volume" followed by separator
  std::string volInfo(const State& state) const {
    return (state.currentVolume != nullptr ? state.currentVolume->volumeName()
                                           : "No Volume") +
           " | ";
  }

  /// @brief Get the logger instance
  /// @return Reference to the logger instance
  const Logger& logger() const { return *m_logger; }

  /// Configuration object for this navigator
  Config m_cfg;

  /// Logger instance for this navigator
  std::unique_ptr<const Logger> m_logger;
};

/// @brief Alternative @c Navigator which tries all possible intersections
///
/// See @c Navigator for more general information about the Navigator concept.
///
/// This Navigator tries all possible intersections with all surfaces in the
/// current volume. It does not use any information about the geometry to
/// optimize the search. It is therefore very slow, but can be used as a
/// reference implementation.
///
class TryAllNavigator final : public TryAllNavigatorBase {
 public:
  /// Type alias for navigator configuration
  using Config = TryAllNavigatorBase::Config;
  /// Type alias for navigator options
  using Options = TryAllNavigatorBase::Options;

  /// @brief Nested State struct
  struct State : public TryAllNavigatorBase::State {
    /// @brief Constructor for navigator state
    /// @param options_ Navigator options to initialize state with
    explicit State(const Options& options_)
        : TryAllNavigatorBase::State(options_) {}

    /// Current navigation candidates with intersection information
    std::vector<NavigationTarget> currentTargets;
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param logger a logger instance
  explicit TryAllNavigator(Config cfg, std::unique_ptr<const Logger> logger =
                                           getDefaultLogger("TryAllNavigator",
                                                            Logging::INFO))
      : TryAllNavigatorBase(std::move(cfg), std::move(logger)) {}

  /// @brief Creates a new navigator state
  /// @param options Navigator options for state initialization
  /// @return Initialized navigator state with current candidates storage
  State makeState(const Options& options) const {
    State state(options);
    return state;
  }

  using TryAllNavigatorBase::currentSurface;
  using TryAllNavigatorBase::currentVolume;
  using TryAllNavigatorBase::currentVolumeMaterial;
  using TryAllNavigatorBase::endOfWorldReached;
  using TryAllNavigatorBase::navigationBreak;
  using TryAllNavigatorBase::startSurface;
  using TryAllNavigatorBase::targetSurface;

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
    auto baseRes = TryAllNavigatorBase::initialize(state, position, direction,
                                                   propagationDirection);
    if (!baseRes.ok()) {
      return baseRes.error();
    }

    // Initialize navigation candidates for the start volume
    reinitializeCandidates(state);

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

    double nearLimit = state.options.nearLimit;
    double farLimit = state.options.farLimit;

    // handle overstepping
    if (!state.currentTargets.empty()) {
      const NavigationTarget& previousTarget = state.currentTargets.front();

      const Surface& surface = previousTarget.surface();
      IntersectionIndex index = previousTarget.intersectionIndex();
      BoundaryTolerance boundaryTolerance = previousTarget.boundaryTolerance();

      auto intersection =
          surface
              .intersect(state.options.geoContext, position, direction,
                         boundaryTolerance, state.options.surfaceTolerance)
              .at(index);

      if (intersection.pathLength() < 0) {
        nearLimit = std::min(nearLimit, intersection.pathLength() -
                                            state.options.surfaceTolerance);
        farLimit = -state.options.surfaceTolerance;

        ACTS_VERBOSE(volInfo(state)
                     << "handle overstepping with nearLimit " << nearLimit
                     << " and farLimit " << farLimit);
      }
    }

    std::vector<NavigationTarget> nextTargets;

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
        nextTargets.emplace_back(
            candidate.target(intersection, intersectionIndex));
      }
    }

    std::ranges::sort(nextTargets, NavigationTarget::pathLengthOrder);

    ACTS_VERBOSE(volInfo(state)
                 << "found " << nextTargets.size() << " intersections");

    NavigationTarget nextTarget = NavigationTarget::None();
    state.currentTargets.clear();

    for (const auto& target : nextTargets) {
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

    state.currentTargets = std::move(nextTargets);

    if (nextTarget.isNone()) {
      ACTS_VERBOSE(volInfo(state) << "no target found");
    } else {
      ACTS_VERBOSE(volInfo(state)
                   << "next target is " << nextTarget.surface().geometryId());
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

    if (state.currentTargets.empty()) {
      ACTS_VERBOSE(volInfo(state) << "No current target set.");
      return;
    }

    assert(state.currentSurface == nullptr && "Current surface must be reset.");

    // handle multiple surface intersections due to increased bounds

    std::vector<NavigationTarget> hitTargets;

    for (const auto& target : state.currentTargets) {
      std::uint8_t index = target.intersectionIndex();
      const Surface& surface = target.surface();
      BoundaryTolerance boundaryTolerance = BoundaryTolerance::None();

      Intersection3D intersection =
          surface
              .intersect(state.options.geoContext, position, direction,
                         boundaryTolerance, state.options.surfaceTolerance)
              .at(index);

      if (intersection.status() == IntersectionStatus::onSurface) {
        hitTargets.emplace_back(target);
      }
    }

    state.currentTargets.clear();

    ACTS_VERBOSE(volInfo(state)
                 << "Found " << hitTargets.size()
                 << " intersections on surface with bounds check.");

    if (hitTargets.empty()) {
      ACTS_VERBOSE(volInfo(state) << "No hit targets found.");
      return;
    }

    if (hitTargets.size() > 1) {
      ACTS_VERBOSE(volInfo(state)
                   << "Only using first intersection within bounds.");
    }

    // we can only handle a single surface hit so we pick the first one
    const auto target = hitTargets.front();
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
  /// Helper method to reset and reinitialize the navigation candidates.
  void reinitializeCandidates(State& state) const {
    state.navigationCandidates.clear();
    state.currentTargets.clear();

    initializeVolumeCandidates(state);
  }
};

/// @brief Alternative @c Navigator which tries all possible intersections
///
/// See @c Navigator for more general information about the Navigator concept.
///
/// This Navigator tries all possible intersections with all surfaces in the
/// current volume. It does not use any information about the geometry to
/// optimize the search. It is therefore very slow, but can be used as a
/// reference implementation.
///
/// Different to @c TryAllNavigator, this Navigator discovers intersections by
/// stepping forward blindly and then checking for intersections with all
/// surfaces in the current volume. This is slower, but more robust against
/// bent tracks.
///
class TryAllOverstepNavigator final : public TryAllNavigatorBase {
 public:
  /// Type alias for navigator configuration
  using Config = TryAllNavigatorBase::Config;

  /// Type alias for navigator options
  using Options = TryAllNavigatorBase::Options;

  /// @brief Nested State struct
  ///
  /// It acts as an internal state which is created for every propagation and
  /// meant to keep thread-local navigation information.
  struct State : public TryAllNavigatorBase::State {
    /// @brief Constructor for overstep navigator state
    /// @param options_ Navigator options to initialize state with
    explicit State(const Options& options_)
        : TryAllNavigatorBase::State(options_) {}

    /// The vector of active targets to work through
    std::vector<NavigationTarget> activeTargets;
    /// The current active target index of the navigation state
    int activeTargetIndex = -1;

    /// The position before the last step
    std::optional<Vector3> lastPosition;

    /// Provides easy access to the active intersection target
    /// @return Reference to the currently active intersection candidate
    const NavigationTarget& activeTarget() const {
      return activeTargets.at(activeTargetIndex);
    }

    /// @brief Check if all navigation candidates have been processed
    /// @return True if no more candidates are available for navigation
    bool endOfTargets() const {
      return activeTargetIndex >= static_cast<int>(activeTargets.size());
    }
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param logger a logger instance
  explicit TryAllOverstepNavigator(
      Config cfg, std::unique_ptr<const Logger> logger = getDefaultLogger(
                      "TryAllOverstepNavigator", Logging::INFO))
      : TryAllNavigatorBase(std::move(cfg), std::move(logger)) {}

  /// @brief Creates a new overstep navigator state
  /// @param options Navigator options for state initialization
  /// @return Initialized navigator state with active candidates and position tracking
  State makeState(const Options& options) const {
    State state(options);
    return state;
  }

  using TryAllNavigatorBase::currentSurface;
  using TryAllNavigatorBase::currentVolume;
  using TryAllNavigatorBase::currentVolumeMaterial;
  using TryAllNavigatorBase::endOfWorldReached;
  using TryAllNavigatorBase::navigationBreak;
  using TryAllNavigatorBase::startSurface;
  using TryAllNavigatorBase::targetSurface;

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
    auto baseRes = TryAllNavigatorBase::initialize(state, position, direction,
                                                   propagationDirection);
    if (!baseRes.ok()) {
      return baseRes.error();
    }

    // Initialize navigation candidates for the start volume
    reinitializeCandidates(state);

    state.lastPosition.reset();

    return Result<void>::success();
  }

  /// @brief Get the next target surface
  ///
  /// This method gets the next target surface based on the current
  /// position and direction. It returns an invalid target if no target can be
  /// found.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return The next target surface
  NavigationTarget nextTarget(State& state, const Vector3& position,
                              const Vector3& direction) const {
    static_cast<void>(direction);

    // Navigator preStep always resets the current surface
    state.currentSurface = nullptr;

    // Check if the navigator is inactive
    if (state.navigationBreak) {
      return NavigationTarget::None();
    }

    ACTS_VERBOSE(volInfo(state) << "nextTarget");

    // We cannot do anything without a last position
    if (!state.lastPosition.has_value() && state.endOfTargets()) {
      ACTS_VERBOSE(
          volInfo(state)
          << "Initial position, nothing to do, blindly stepping forward.");
      state.lastPosition = position;
      return NavigationTarget::None();
    }

    if (state.endOfTargets()) {
      ACTS_VERBOSE(volInfo(state) << "evaluate blind step");

      Vector3 stepStart = state.lastPosition.value();
      Vector3 stepEnd = position;
      Vector3 step = stepEnd - stepStart;
      double stepDistance = step.norm();
      if (stepDistance < std::numeric_limits<double>::epsilon()) {
        ACTS_ERROR(volInfo(state) << "Step distance is zero. " << stepDistance);
      }
      Vector3 stepDirection = step.normalized();

      double nearLimit = -stepDistance + state.options.surfaceTolerance;
      double farLimit = 0;

      state.lastPosition.reset();
      state.activeTargets.clear();
      state.activeTargetIndex = -1;

      // Find intersections with all candidates
      for (const auto& candidate : state.navigationCandidates) {
        auto intersections =
            candidate.intersect(state.options.geoContext, stepEnd,
                                stepDirection, state.options.surfaceTolerance);
        for (auto [intersectionIndex, intersection] :
             Acts::enumerate(intersections)) {
          // exclude invalid intersections
          if (!intersection.isValid() ||
              !detail::checkPathLength(intersection.pathLength(), nearLimit,
                                       farLimit)) {
            continue;
          }
          // store candidate
          state.activeTargets.emplace_back(
              candidate.target(intersection, intersectionIndex));
        }
      }

      std::ranges::sort(state.activeTargets, NavigationTarget::pathLengthOrder);

      ACTS_VERBOSE(volInfo(state) << "Found " << state.activeTargets.size()
                                  << " intersections");

      for (const auto& target : state.activeTargets) {
        ACTS_VERBOSE("found target " << target.surface().geometryId());
      }
    }

    ++state.activeTargetIndex;

    if (state.endOfTargets()) {
      ACTS_VERBOSE(volInfo(state)
                   << "No target found, blindly stepping forward.");
      state.lastPosition = position;
      return NavigationTarget::None();
    }

    ACTS_VERBOSE(volInfo(state) << "handle active candidates");

    ACTS_VERBOSE(volInfo(state)
                 << (state.activeTargets.size() - state.activeTargetIndex)
                 << " out of " << state.activeTargets.size()
                 << " surfaces remain to try.");

    const auto& target = state.activeTarget();

    ACTS_VERBOSE(volInfo(state) << "Next surface candidate will be "
                                << target.surface().geometryId());

    return target;
  }

  /// @brief Check if the target is still valid
  ///
  /// This method checks if the target is valid based on the current position
  /// and direction. It returns true if the target is still valid.
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

    return true;
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
    if (state.navigationBreak) {
      return;
    }

    ACTS_VERBOSE(volInfo(state) << "handleSurfaceReached");

    assert(state.currentSurface == nullptr && "Current surface must be reset.");

    if (state.endOfTargets()) {
      ACTS_VERBOSE(volInfo(state) << "No active candidate set.");
      return;
    }

    std::vector<NavigationTarget> hitTargets;

    while (!state.endOfTargets()) {
      const auto& candidate = state.activeTarget();
      IntersectionIndex index = candidate.intersectionIndex();
      const Surface& surface = candidate.surface();
      BoundaryTolerance boundaryTolerance = candidate.boundaryTolerance();

      // first with boundary tolerance
      IntersectionStatus surfaceStatus =
          surface
              .intersect(state.options.geoContext, position, direction,
                         boundaryTolerance, state.options.surfaceTolerance)
              .at(index)
              .status();

      if (surfaceStatus != IntersectionStatus::onSurface) {
        break;
      }

      // now without boundary tolerance
      boundaryTolerance = BoundaryTolerance::None();
      surfaceStatus =
          surface
              .intersect(state.options.geoContext, position, direction,
                         boundaryTolerance, state.options.surfaceTolerance)
              .at(index)
              .status();

      if (surfaceStatus == IntersectionStatus::onSurface) {
        hitTargets.emplace_back(candidate);
      }

      ++state.activeTargetIndex;
      ACTS_VERBOSE("skip target " << surface.geometryId());
    }

    // we increased the target index one too many times
    --state.activeTargetIndex;

    ACTS_VERBOSE(volInfo(state)
                 << "Found " << hitTargets.size()
                 << " intersections on surface with bounds check.");

    if (hitTargets.empty()) {
      ACTS_VERBOSE(volInfo(state)
                   << "Surface successfully hit, but outside bounds.");
      return;
    }

    if (hitTargets.size() > 1) {
      ACTS_VERBOSE(volInfo(state)
                   << "Only using first intersection within bounds.");
    }

    // we can only handle a single surface hit so we pick the first one
    const auto& candidate = hitTargets.front();
    const Surface& surface = candidate.surface();

    ACTS_VERBOSE(volInfo(state) << "Surface successfully hit, storing it.");
    state.currentSurface = &surface;

    if (state.currentSurface != nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Current surface set to surface "
                                  << surface.geometryId());
    }

    if (candidate.isSurfaceTarget()) {
      ACTS_VERBOSE(volInfo(state) << "This is a surface");
    } else if (candidate.isLayerTarget()) {
      ACTS_VERBOSE(volInfo(state) << "This is a layer");
    } else if (candidate.isPortalTarget()) {
      ACTS_VERBOSE(volInfo(state)
                   << "This is a portal. Reinitialize navigation");

      const BoundarySurface& boundary = candidate.boundarySurface();

      state.currentVolume = boundary.attachedVolume(state.options.geoContext,
                                                    position, direction);

      ACTS_VERBOSE(volInfo(state) << "Switched volume");

      reinitializeCandidates(state);
    } else {
      ACTS_ERROR(volInfo(state) << "Unknown intersection type");
    }
  }

 private:
  /// Helper method to reset and reinitialize the navigation candidates.
  void reinitializeCandidates(State& state) const {
    state.navigationCandidates.clear();
    state.activeTargets.clear();
    state.activeTargetIndex = -1;

    initializeVolumeCandidates(state);
  }
};

}  // namespace Acts
