// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
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

class IVolumeMaterial;

/// @brief A navigator without a tracking geometry
///
/// It offers the additional surfaces of the options and the target surface,
/// which is one of them. Without either it does nothing. The navigators with a
/// tracking geometry hold one to offer the additional surfaces on top of the
/// geometry.
///
/// The additional surfaces and the candidates of the other navigator are two
/// streams. Both are intersected again before every step, and the closer one
/// is the target. A surface the last step crossed, ahead before it and behind
/// after it, is offered behind the position, so the propagation steps back
/// onto it.
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
      /// Whether the surface is ahead of the last position the surfaces were
      /// intersected from
      bool ahead = false;
      /// The same for the position before
      bool aheadBefore = false;
      /// Whether the last position the surfaces were intersected from is on
      /// the surface
      bool on = false;
      /// The same for the position before
      bool onBefore = false;
    };

    /// Additional surfaces of the options and the target surface
    std::vector<AdditionalSurfaceState> additionalSurfaces;

    /// The candidate of another navigator an additional surface took
    /// precedence over. It is handed out once the additional surface is no
    /// longer the closer one.
    std::optional<NavigationTarget> heldBack;

    /// The target handed out to the propagation
    std::optional<NavigationTarget> handedOut;
    /// Whether the target handed out is an additional surface
    bool handedOutIsAdditional = false;

    /// Last position the additional surfaces were intersected from
    std::optional<Vector3> lastPosition;
    /// The position before
    std::optional<Vector3> positionBefore;

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
    state.handedOut.reset();
    state.handedOutIsAdditional = false;
    state.lastPosition.reset();
    state.positionBefore.reset();

    return Result<void>::success();
  }

  /// Returns the closest additional surface
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  /// @return The next target, or none
  NavigationTarget nextTarget(
      State& state, const Vector3& position, const Vector3& direction,
      const Vector3& /*curvature*/ = Vector3::Zero()) const {
    state.currentSurface = nullptr;
    return handOut(state, closest(state, position, direction, nullptr), true);
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

    // Ask the other navigator again if it had no candidate
    if (!state.heldBack.has_value() || state.heldBack->isNone()) {
      state.heldBack = nextCandidate();
    }

    if (const NavigationTarget additional =
            closest(state, position, direction, volume);
        !additional.isNone()) {
      NavigationTarget& candidate = state.heldBack.value();
      refresh(state, candidate, position, direction);
      // A tie goes to the candidate of the other navigator
      if (candidate.isNone() ||
          additional.pathLength() < candidate.intersection().pathLength()) {
        return handOut(state, additional, true);
      }
    }

    const NavigationTarget candidate = state.heldBack.value();
    state.heldBack.reset();
    return handOut(state, candidate, false);
  }

  /// Checks if the target is still the closest of the two streams
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  /// @return True if the target handed out is still the closest
  bool checkTargetValid(State& state, const Vector3& position,
                        const Vector3& direction,
                        const Vector3& /*curvature*/ = Vector3::Zero()) const {
    return checkAdditionalTargetValid(state, position, direction, nullptr,
                                      true);
  }

  /// Checks the target handed out against the other stream, for another
  /// navigator holding this one. Intersects the additional surfaces and the
  /// held back candidate again, so the target is swapped as soon as the other
  /// stream is closer, and a crossed additional surface is found.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  /// @param volume The current volume of the other navigator
  /// @param candidateValid Whether the other navigator holds its candidates
  ///        valid. If not, it is asked for a new one.
  /// @return True if the target handed out is still the closest
  bool checkAdditionalTargetValid(State& state, const Vector3& position,
                                  const Vector3& direction,
                                  const TrackingVolume* volume,
                                  bool candidateValid) const {
    if (!candidateValid) {
      state.heldBack.reset();
      state.handedOut.reset();
      return false;
    }
    if (state.additionalSurfaces.empty() || !state.handedOut.has_value()) {
      return true;
    }

    const NavigationTarget additional =
        closest(state, position, direction, volume);
    NavigationTarget& handedOut = state.handedOut.value();

    if (state.handedOutIsAdditional) {
      if (additional.isNone() ||
          &additional.surface() != &handedOut.surface() ||
          additional.intersectionIndex() != handedOut.intersectionIndex()) {
        return false;
      }
      if (state.heldBack.has_value() && !state.heldBack->isNone()) {
        NavigationTarget& candidate = state.heldBack.value();
        refresh(state, candidate, position, direction);
        // A tie goes to the candidate of the other navigator
        if (candidate.intersection().pathLength() <= additional.pathLength()) {
          return false;
        }
      }
      return true;
    }

    if (additional.isNone()) {
      return true;
    }
    refresh(state, handedOut, position, direction);
    if (additional.pathLength() < handedOut.intersection().pathLength()) {
      // Hold the candidate back instead of asking the other navigator again,
      // which would move on to its next candidate
      state.heldBack = handedOut;
      state.handedOut.reset();
      return false;
    }
    return true;
  }

  /// Handles reaching an additional surface
  /// @param state The navigation state
  /// @param surface The surface reached
  void handleSurfaceReached(
      State& state, const Vector3& /*position*/, const Vector3& /*direction*/,
      const Surface& surface,
      const Vector3& /*curvature*/ = Vector3::Zero()) const {
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
    return state.handedOutIsAdditional;
  }

 private:
  /// Slack of the crossing check, relative to the step. Covers the sagitta of
  /// a curved step.
  static constexpr double s_crossingSlack = 0.05;

  /// Records the target handed out
  static NavigationTarget handOut(State& state, const NavigationTarget& target,
                                  bool isAdditional) {
    if (target.isNone()) {
      state.handedOut.reset();
    } else {
      state.handedOut = target;
    }
    state.handedOutIsAdditional = isAdditional && !target.isNone();
    return target;
  }

  /// Intersects a candidate of the other navigator again. Keeps the stored
  /// intersection if the straight-line estimate lost the surface.
  static void refresh(const State& state, NavigationTarget& candidate,
                      const Vector3& position, const Vector3& direction) {
    if (candidate.isNone()) {
      return;
    }
    const Intersection3D refreshed =
        candidate.surface()
            .intersect(state.options.geoContext, position, direction,
                       candidate.boundaryTolerance(),
                       state.options.surfaceTolerance)
            .at(candidate.intersectionIndex());
    if (refreshed.isValid()) {
      candidate.intersection() = refreshed;
    }
  }

  /// Whether the last step crossed the surface at @p intersection, which is
  /// behind the position: the surface was ahead before the step, the step did
  /// not start on it, and the intersection lies on the step
  static bool crossed(const State& state,
                      const State::AdditionalSurfaceState& additional,
                      const Intersection3D& intersection) {
    if (!additional.aheadBefore || additional.onBefore ||
        !state.positionBefore.has_value()) {
      return false;
    }
    const Vector3& before = state.positionBefore.value();
    const Vector3& after = state.lastPosition.value();
    const double step = (after - before).norm();
    const double along = (intersection.position() - before).norm() +
                         (intersection.position() - after).norm();
    return along <=
           (1 + s_crossingSlack) * step + state.options.surfaceTolerance;
  }

  /// Closest additional surface ahead of the position, or crossed by the last
  /// step, in @p volume for a surface restricted to a volume
  NavigationTarget closest(State& state, const Vector3& position,
                           const Vector3& direction,
                           const TrackingVolume* volume) const {
    // Intersected from a new position: the last one becomes the one before
    const bool moved =
        !state.lastPosition.has_value() || state.lastPosition != position;
    if (moved) {
      state.positionBefore = state.lastPosition;
      state.lastPosition = position;
    }

    NavigationTarget closest = NavigationTarget::None();

    for (State::AdditionalSurfaceState& additional : state.additionalSurfaces) {
      if (moved) {
        additional.aheadBefore = additional.ahead;
        additional.ahead = false;
        additional.onBefore = additional.on;
        additional.on = false;
      }

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

      for (const auto [intersectionIndex, intersection] :
           enumerate(multiIntersection)) {
        if (!intersection.isValid() ||
            intersection.pathLength() > state.options.farLimit) {
          continue;
        }
        const double pathLength = intersection.pathLength();
        const bool onSurface =
            intersection.status() == IntersectionStatus::onSurface;
        if (onSurface) {
          additional.on = true;
        } else if (pathLength > 0) {
          additional.ahead = true;
        }

        // Offer what is ahead, what is on the position unless reached, so a
        // kept surface is not reached again without a step, and what the last
        // step crossed
        const bool offer =
            onSurface
                ? !additional.reached
                : pathLength > 0 || crossed(state, additional, intersection);
        if (!offer) {
          continue;
        }
        if (closest.isNone() || pathLength < closest.pathLength()) {
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
