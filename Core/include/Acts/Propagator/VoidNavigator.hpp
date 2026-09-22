// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/NavigatorInitializeArguments.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Propagator/NavigatorStatistics.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"

#include <optional>
#include <stdexcept>
#include <vector>

namespace Acts {

class TrackingVolume;
class IVolumeMaterial;

/// @brief A navigator without a tracking geometry
///
/// It offers the additional surfaces of the options and the target surface,
/// which is one of them. Without either it does nothing. The navigators with a
/// tracking geometry hold one to offer the additional surfaces on top of the
/// geometry.
class VoidNavigator {
 public:
  /// Nested Config struct
  struct Config {};

  /// Nested Options struct
  struct Options : public NavigatorPlainOptions {
    /// Constructor for void navigator options
    /// @param gctx Geometry context
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    /// Constructor from the plain options of another navigator
    /// @param options The plain navigator options to copy
    explicit Options(const NavigatorPlainOptions& options)
        : NavigatorPlainOptions(options) {}

    /// Sets the plain navigator options
    /// @param options The plain navigator options to copy
    void setPlainOptions(const NavigatorPlainOptions& options) {
      static_cast<NavigatorPlainOptions&>(*this) = options;
    }
  };

  /// Nested State struct
  struct State {
    /// Constructor for void navigator state
    /// @param options_ The navigator options to store in state
    explicit State(const Options& options_) : options(options_) {}

    /// Configuration options for the navigator
    Options options;

    /// Surface that is the target of the navigation
    const Surface* targetSurface = nullptr;

    /// The additional surface the propagation is on
    const Surface* currentSurface = nullptr;

    /// An additional surface with its bookkeeping
    struct AdditionalSurfaceState {
      /// The additional surface
      AdditionalSurface surface;
      /// Whether the propagation reached the surface
      bool reached = false;
    };

    /// Additional surfaces of the options and the target surface
    std::vector<AdditionalSurfaceState> additionalSurfaces;

    /// The candidate of another navigator an additional surface took
    /// precedence over. It is handed out once the additional surface is no
    /// longer the closer one.
    std::optional<NavigationTarget> heldBack;

    /// Navigation statistics
    NavigatorStatistics statistics;
  };

  /// Creates a new navigator state for void navigation
  /// @param options The navigator options
  /// @return Initialized void navigator state
  State makeState(const Options& options) const {
    State state(options);
    return state;
  }

  /// Returns the additional surface the propagation is on
  /// @param state The navigation state
  /// @return The current surface, nullptr if none
  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  /// Returns the current tracking volume (always nullptr for void navigator)
  /// @return Always nullptr since void navigator has no volumes
  const TrackingVolume* currentVolume(const State& /*state*/) const {
    return nullptr;
  }

  /// Returns the current volume material (always nullptr for void navigator)
  /// @return Always nullptr since void navigator has no material
  const IVolumeMaterial* currentVolumeMaterial(const State& /*state*/) const {
    return nullptr;
  }

  /// Returns the start surface (always nullptr for void navigator)
  /// @return Always nullptr
  const Surface* startSurface(const State& /*state*/) const { return nullptr; }

  /// Returns the target surface
  /// @param state The navigation state
  /// @return The target surface, nullptr if none
  const Surface* targetSurface(const State& state) const {
    return state.targetSurface;
  }

  /// Checks if the end of the world has been reached (always false for void
  /// navigator)
  /// @return Always false since void navigator has no world boundaries
  bool endOfWorldReached(const State& /*state*/) const { return false; }

  /// Checks if navigation should break (always true for void navigator)
  /// @return Always true since void navigator has no geometry to navigate
  bool navigationBreak(const State& /*state*/) const { return true; }

  /// Initializes the void navigator for a new run
  /// @param state The navigation state
  /// @param args The initialization arguments of this navigation run
  /// @return Always successful result
  [[nodiscard]] Result<void> initialize(
      State& state, const NavigatorInitializeArguments& args) const {
    state.targetSurface = args.targetSurface;
    state.currentSurface = nullptr;

    state.additionalSurfaces.clear();
    state.additionalSurfaces.reserve(state.options.additionalSurfaces.size() +
                                     1);
    for (const AdditionalSurface& surface : state.options.additionalSurfaces) {
      if (surface.surface == nullptr) {
        throw std::invalid_argument("Navigator: additional surface is nullptr");
      }
      state.additionalSurfaces.push_back({surface});
    }
    // The target is an additional surface with the bounds check the target
    // aborter applies
    if (args.targetSurface != nullptr) {
      state.additionalSurfaces.push_back(
          {{args.targetSurface, BoundaryTolerance::None()}});
    }

    state.heldBack.reset();

    return Result<void>::success();
  }

  /// Returns the closest additional surface
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  /// @return The next target, or none
  NavigationTarget nextTarget(State& state, const Vector3& position,
                              const Vector3& direction) const {
    state.currentSurface = nullptr;
    return closest(state, position, direction, nullptr);
  }

  /// Returns the closer of the next candidate of another navigator and the
  /// closest additional surface. The candidate is held back, not skipped,
  /// while an additional surface is closer.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  /// @param volume The current volume of the other navigator, read after
  ///        @p nextCandidate ran
  /// @param nextCandidate Callable returning the next candidate of the other
  ///        navigator
  /// @return The next target
  template <typename next_candidate_t>
  NavigationTarget nextTarget(State& state, const Vector3& position,
                              const Vector3& direction,
                              const TrackingVolume* const& volume,
                              next_candidate_t&& nextCandidate) const {
    state.currentSurface = nullptr;

    if (state.additionalSurfaces.empty()) {
      return nextCandidate();
    }

    // The other navigator can have no candidate for one step only
    if (!state.heldBack.has_value() || state.heldBack->isNone()) {
      state.heldBack = nextCandidate();
    }

    if (const NavigationTarget additional =
            closest(state, position, direction, volume);
        !additional.isNone()) {
      NavigationTarget& candidate = state.heldBack.value();
      if (!candidate.isNone()) {
        // The stored path length is stale by the distance travelled since
        const Intersection3D refreshed =
            candidate.surface()
                .intersect(state.options.geoContext, position, direction,
                           candidate.boundaryTolerance(),
                           state.options.surfaceTolerance)
                .at(candidate.intersectionIndex());
        // Keep the stored one if the straight-line estimate lost the surface
        if (refreshed.isValid()) {
          candidate.intersection() = refreshed;
        }
      }
      // A tie goes to the candidate of the other navigator
      if (candidate.isNone() ||
          additional.pathLength() < candidate.intersection().pathLength()) {
        return additional;
      }
    }

    const NavigationTarget candidate = state.heldBack.value();
    state.heldBack.reset();
    return candidate;
  }

  /// Checks if the current target is valid (always true for void navigator)
  /// @return Always true since the stepper re-intersects the target
  bool checkTargetValid(const State& /*state*/, const Vector3& /*position*/,
                        const Vector3& /*direction*/) const {
    return true;
  }

  /// Handles reaching an additional surface
  /// @param state The navigation state
  /// @param surface The surface reached
  void handleSurfaceReached(State& state, const Vector3& /*position*/,
                            const Vector3& /*direction*/,
                            const Surface& surface) const {
    handleAdditionalSurfaceReached(state, surface);
  }

  /// Handles reaching a surface for another navigator
  /// @param state The navigation state
  /// @param surface The surface reached
  /// @return True if the surface is an additional surface the other navigator
  ///         did not target, so the other navigator has nothing to do
  bool handleAdditionalSurfaceReached(State& state,
                                      const Surface& surface) const {
    bool found = false;
    for (State::AdditionalSurfaceState& additional : state.additionalSurfaces) {
      if (additional.surface.surface == &surface) {
        additional.reached = true;
        found = true;
      }
    }
    if (!found) {
      return false;
    }
    state.currentSurface = &surface;
    return state.heldBack.has_value();
  }

 private:
  /// Closest additional surface ahead of the position
  NavigationTarget closest(const State& state, const Vector3& position,
                           const Vector3& direction,
                           const TrackingVolume* volume) const {
    NavigationTarget closest = NavigationTarget::None();

    for (const State::AdditionalSurfaceState& additional :
         state.additionalSurfaces) {
      const AdditionalSurface& surface = additional.surface;
      if (additional.reached && surface.dropAfterReached) {
        continue;
      }
      if (surface.volume != nullptr && surface.volume != volume) {
        continue;
      }

      // The closer solution can be behind the propagation, so check all
      const MultiIntersection3D multiIntersection = surface.surface->intersect(
          state.options.geoContext, position, direction,
          surface.boundaryTolerance, state.options.surfaceTolerance);

      // The stepper decides whether the surface is reached, so everything
      // ahead and on the position is offered
      for (const auto [intersectionIndex, intersection] :
           enumerate(multiIntersection)) {
        if (!intersection.isValid() ||
            (intersection.status() != IntersectionStatus::onSurface &&
             intersection.pathLength() <= 0) ||
            intersection.pathLength() > state.options.farLimit) {
          continue;
        }
        if (closest.isNone() ||
            intersection.pathLength() < closest.pathLength()) {
          closest = NavigationTarget(
              intersection, static_cast<IntersectionIndex>(intersectionIndex),
              *surface.surface, surface.boundaryTolerance);
        }
      }
    }

    return closest;
  }
};

}  // namespace Acts
