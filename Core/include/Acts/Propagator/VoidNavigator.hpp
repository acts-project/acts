// This file is part of the Acts project.
//
// Copyright (C) 2018-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

class Surface;

/// @brief The void navigator struct as a default navigator
///
/// It does not provide any navigation action, the compiler
/// should eventually optimise that the function call is not done
///
struct VoidNavigator {
  struct Config {};

  struct Options : public NavigatorPlainOptions {
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    void setPlainOptions(const NavigatorPlainOptions& options) {
      static_cast<NavigatorPlainOptions&>(*this) = options;
    }
  };

  /// @brief Nested State struct, minimal requirement
  struct State {
    explicit State(const Options& options_) : options(options_) {}

    Options options;
  };

  State makeState(const Options& options) const {
    State state(options);
    return state;
  }

  const Surface* currentSurface(const State& /*state*/) const {
    return nullptr;
  }

  const Surface* startSurface(const State& /*state*/) const { return nullptr; }

  const Surface* targetSurface(const State& /*state*/) const { return nullptr; }

  bool targetReached(const State& /*state*/) const { return false; }

  bool navigationBreak(const State& /*state*/) const { return false; }

  void currentSurface(State& /*state*/, const Surface* /*surface*/) const {}

  void targetReached(State& /*state*/, bool /*targetReached*/) const {}

  void navigationBreak(State& /*state*/, bool /*navigationBreak*/) const {}

  void initialize(State& /*state*/, const Vector3& /*position*/,
                  const Vector3& /*direction*/) const {}

  SurfaceIntersection estimateNextTarget(State& /*state*/,
                                         const Vector3& /*position*/,
                                         const Vector3& /*direction*/) const {
    return SurfaceIntersection::invalid();
  }

  void registerSurfaceStatus(State& /*state*/, const Vector3& /*position*/,
                             const Vector3& /*direction*/,
                             const Surface& /*surface*/,
                             IntersectionStatus /*surfaceStatus*/) const {}
};

}  // namespace Acts
