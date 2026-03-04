// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Propagator/NavigatorStatistics.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

class TrackingVolume;
class IVolumeMaterial;
class Surface;

/// @brief A navigator that does nothing
///
/// It does not provide any navigation action
///
class VoidNavigator {
 public:
  /// @brief Nested Config struct
  struct Config {};

  /// @brief Nested Options struct
  struct Options : public NavigatorPlainOptions {
    /// @brief Constructor for void navigator options
    /// @param gctx Geometry context (required but unused by void navigator)
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    /// @brief Sets the plain navigator options
    /// @param options The plain navigator options to copy (unused by void navigator)
    void setPlainOptions(const NavigatorPlainOptions& options) {
      static_cast<NavigatorPlainOptions&>(*this) = options;
    }
  };

  /// @brief Nested State struct
  struct State {
    /// @brief Constructor for void navigator state
    /// @param options_ The navigator options to store in state
    explicit State(const Options& options_) : options(options_) {}

    /// Configuration options for the navigator
    Options options;

    /// Navigation statistics
    NavigatorStatistics statistics;
  };

  /// @brief Creates a new navigator state for void navigation
  /// @param options The navigator options
  /// @return Initialized void navigator state
  State makeState(const Options& options) const {
    State state(options);
    return state;
  }

  /// @brief Returns the current surface (always nullptr for void navigator)
  /// @return Always nullptr since void navigator has no surfaces
  const Surface* currentSurface(const State& /*state*/) const {
    return nullptr;
  }

  /// @brief Returns the current tracking volume (always nullptr for void navigator)
  /// @return Always nullptr since void navigator has no volumes
  const TrackingVolume* currentVolume(const State& /*state*/) const {
    return nullptr;
  }

  /// @brief Returns the current volume material (always nullptr for void navigator)
  /// @return Always nullptr since void navigator has no material
  const IVolumeMaterial* currentVolumeMaterial(const State& /*state*/) const {
    return nullptr;
  }

  /// @brief Returns the start surface (always nullptr for void navigator)
  /// @return Always nullptr since void navigator has no surfaces
  const Surface* startSurface(const State& /*state*/) const { return nullptr; }

  /// @brief Returns the target surface (always nullptr for void navigator)
  /// @return Always nullptr since void navigator has no surfaces
  const Surface* targetSurface(const State& /*state*/) const { return nullptr; }

  /// @brief Checks if navigation should break (always true for void navigator)
  /// @return Always true to immediately stop navigation
  bool navigationBreak(const State& /*state*/) const { return true; }

  /// @brief Initializes the void navigator (always succeeds and does nothing)
  /// @return Always successful result since no initialization is needed
  [[nodiscard]] Result<void> initialize(
      State& /*state*/, const Vector3& /*position*/,
      const Vector3& /*direction*/, Direction /*propagationDirection*/) const {
    return Result<void>::success();
  }

  /// @brief Returns the next navigation target (always None for void navigator)
  /// @return NavigationTarget::None() since there are no targets in void space
  NavigationTarget nextTarget(State& /*state*/, const Vector3& /*position*/,
                              const Vector3& /*direction*/) const {
    return NavigationTarget::None();
  }

  /// @brief Checks if the current target is valid (always true for void navigator)
  /// @return Always true since there are no targets to invalidate
  bool checkTargetValid(const State& /*state*/, const Vector3& /*position*/,
                        const Vector3& /*direction*/) const {
    return true;
  }

  /// @brief Handles reaching a surface (does nothing for void navigator)
  void handleSurfaceReached(State& /*state*/, const Vector3& /*position*/,
                            const Vector3& /*direction*/,
                            const Surface& /*surface*/) const {
    return;
  }
};

}  // namespace Acts
