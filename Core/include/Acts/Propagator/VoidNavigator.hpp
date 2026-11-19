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
#include "Acts/Propagator/detail/AdditionalSurfaceNavigation.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

class IVolumeMaterial;

/// A navigator without a tracking geometry. It visits the additional surfaces
/// in its options and the target surface, if one is provided.
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

    /// Additional surfaces and candidate history
    detail::AdditionalSurfaceNavigation additional;

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

    state.additional.initialize(state.options, args.targetSurface);
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
    return state.additional.nextTarget(state.options, position, direction);
  }

  /// Checks whether the selected additional surface is still the closest
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  /// @return True if the selected target remains valid
  bool checkTargetValid(State& state, const Vector3& position,
                        const Vector3& direction,
                        const Vector3& /*curvature*/ = Vector3::Zero()) const {
    return state.additional.checkTargetValid(state.options, position, direction,
                                             nullptr, true);
  }

  /// Records reaching an additional surface
  /// @param state The navigation state
  /// @param surface The surface reached
  void handleSurfaceReached(
      State& state, const Vector3& /*position*/, const Vector3& /*direction*/,
      const Surface& surface,
      const Vector3& /*curvature*/ = Vector3::Zero()) const {
    if (state.additional.handleSurfaceReached(surface)) {
      state.currentSurface = &surface;
    }
  }
};

}  // namespace Acts
