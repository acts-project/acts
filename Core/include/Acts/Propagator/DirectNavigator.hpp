// This file is part of the Acts project.
//
// Copyright (C) 2019-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <memory>
#include <vector>

namespace Acts {
class TrackingVolume;

/// This is a fully guided navigator that progresses through a pre-given
/// sequence of surfaces.
///
/// This can either be used as a validation tool, for truth tracking, or track
/// refitting.
class DirectNavigator {
 public:
  /// The sequentially crossed surfaces
  using SurfaceSequence = std::vector<const Surface*>;
  using SurfaceIter = SurfaceSequence::const_iterator;

  struct Config {};

  struct Options : public NavigatorPlainOptions {
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    /// The Surface sequence
    SurfaceSequence surfaces;

    // TODO https://github.com/acts-project/acts/issues/2738
    /// Distance limit to discard intersections "behind us"
    /// @note this is only necessary because some surfaces have more than one
    ///       intersection
    double nearLimit = -100 * UnitConstants::um;

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
    explicit State(const Options& options_) : options(options_) {}

    Options options;

    /// Index of the next surface to try
    std::size_t surfaceIndex = 0;

    const Surface& surface() {
      const Surface* surface = options.surfaces.at(surfaceIndex);
      assert(surface != nullptr && "Surface is nullptr");
      return *surface;
    }

    /// Navigation state - external interface: the current surface
    const Surface* currentSurface = nullptr;

    /// Navigation state - external interface: target is reached
    bool targetReached = false;
    /// Navigation state - external interface: a break has been detected
    bool navigationBreak = false;
  };

  explicit DirectNavigator(std::unique_ptr<const Logger> _logger =
                               getDefaultLogger("DirectNavigator",
                                                Logging::INFO))
      : m_logger{std::move(_logger)} {}

  State makeState(const Options& options) const {
    State state(options);
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
  /// @param [in,out] state is the propagation state object
  void initialize(State& state, const Vector3& /*position*/,
                  const Vector3& /*direction*/) const {
    ACTS_VERBOSE("Initialize. Surface sequence for navigation:");
    for (auto surface : state.options.surfaces) {
      ACTS_VERBOSE(surface->geometryId()
                   << " - "
                   << surface->center(state.options.geoContext).transpose());
    }

    // We set the current surface to the start surface
    state.currentSurface = state.options.startSurface;

    if (state.currentSurface != nullptr) {
      ACTS_VERBOSE("Current surface set to start surface "
                   << state.currentSurface->geometryId());
    }
  }

  /// @brief Navigator pre step call
  ///
  /// @param [in,out] state is the mutable propagator state object
  SurfaceIntersection estimateNextTarget(State& state, const Vector3& position,
                                         const Vector3& direction) const {
    if (state.navigationBreak) {
      return SurfaceIntersection::invalid();
    }

    ACTS_VERBOSE("pre step");

    // Navigator target always resets the current surface
    state.currentSurface = nullptr;

    // Output the position in the sequence
    ACTS_VERBOSE((state.options.surfaces.size() - state.surfaceIndex)
                 << " out of " << state.options.surfaces.size()
                 << " surfaces remain to try.");

    if (state.surfaceIndex >= state.options.surfaces.size()) {
      // Set the navigation break
      state.navigationBreak = true;
      // If no externally provided target is given, the target is reached
      if (state.options.targetSurface == nullptr) {
        state.targetReached = true;
        // Announce it then
        ACTS_VERBOSE("No target Surface, job done.");
      }
      return SurfaceIntersection::invalid();
    }

    while (state.surfaceIndex < state.options.surfaces.size()) {
      const Surface& surface = state.surface();
      const double farLimit = std::numeric_limits<double>::max();
      const SurfaceIntersection intersection = chooseIntersection(
          state.options.geoContext, surface, position, direction,
          BoundaryTolerance::Infinite(), state.options.nearLimit, farLimit,
          s_onSurfaceTolerance);
      if (intersection.status() >= Intersection3D::Status::reachable) {
        ACTS_VERBOSE("Next surface reachable at distance  "
                     << intersection.pathLength());
        return intersection;
      }

      ACTS_VERBOSE(
          "Surface not reachable anymore, switching to next one in "
          "sequence");
      // Move the sequence to the next surface
      ++state.surfaceIndex;
    }

    // Set the navigation break
    state.navigationBreak = true;
    return SurfaceIntersection::invalid();
  }

  void registerSurfaceStatus(State& state, const Vector3& /*position*/,
                             const Vector3& /*direction*/,
                             const Surface& surface,
                             IntersectionStatus surfaceStatus) const {
    if (state.navigationBreak) {
      return;
    }

    ACTS_VERBOSE("post step");

    // Navigator post step always resets the current surface
    state.currentSurface = nullptr;

    // Output the position in the sequence
    ACTS_VERBOSE((state.options.surfaces.size() - state.surfaceIndex)
                 << " out of " << state.options.surfaces.size()
                 << " surfaces remain to try.");

    if (state.surfaceIndex >= state.options.surfaces.size()) {
      return;
    }

    assert(&surface == &state.surface() &&
           "Surfaces do not match, this is a bug");

    if (surfaceStatus == Intersection3D::Status::onSurface) {
      // Set the current surface
      state.currentSurface = &surface;
      ACTS_VERBOSE("Current surface set to  "
                   << state.currentSurface->geometryId());
      // Move the sequence to the next surface
      ++state.surfaceIndex;
      if (state.surfaceIndex < state.options.surfaces.size()) {
        ACTS_VERBOSE("Next surface candidate is  "
                     << state.surface().geometryId());
      }
    }
  }

 private:
  ObjectIntersection<Surface> chooseIntersection(
      const GeometryContext& gctx, const Surface& surface,
      const Vector3& position, const Vector3& direction,
      const BoundaryTolerance& boundaryTolerance, double nearLimit,
      double farLimit, double tolerance) const {
    auto intersections = surface.intersect(gctx, position, direction,
                                           boundaryTolerance, tolerance);

    for (auto& intersection : intersections.split()) {
      if (detail::checkPathLength(intersection.pathLength(), nearLimit,
                                  farLimit, logger())) {
        return intersection;
      }
    }

    return ObjectIntersection<Surface>::invalid();
  }

  const Logger& logger() const { return *m_logger; }

  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts
