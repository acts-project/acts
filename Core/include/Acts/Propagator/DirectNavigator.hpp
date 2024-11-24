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
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <limits>
#include <memory>
#include <vector>

namespace Acts {

class GeometryContext;

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

    /// The surface tolerance
    double surfaceTolerance = s_onSurfaceTolerance;

    // TODO https://github.com/acts-project/acts/issues/2738
    /// Distance limit to discard intersections "behind us"
    /// @note this is only necessary because some surfaces have more than one
    ///       intersection
    double nearLimit = -100 * UnitConstants::um;

    /// The far limit to resolve surfaces
    double farLimit = std::numeric_limits<double>::max();

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

    Direction direction = Direction::Forward;

    /// Index of the next surface to try
    int surfaceIndex = -1;

    /// Navigation state - external interface: the current surface
    const Surface* currentSurface = nullptr;

    /// Navigation state - external interface: a break has been detected
    bool navigationBreak = false;

    /// Navigation statistics
    NavigatorStatistics statistics;

    const Surface& navSurface() const {
      return *options.surfaces.at(surfaceIndex);
    }

    void nextSurface() {
      if (direction == Direction::Forward) {
        ++surfaceIndex;
      } else {
        --surfaceIndex;
      }
    }

    bool endOfSurfaces() const {
      return surfaceIndex < 0 ||
             surfaceIndex >= static_cast<int>(options.surfaces.size());
    }

    int remainingSurfaces() const {
      if (direction == Direction::Forward) {
        return options.surfaces.size() - surfaceIndex;
      }
      return surfaceIndex + 1;
    }

    void resetSurfaceIndex() {
      surfaceIndex =
          direction == Direction::Forward ? -1 : options.surfaces.size();
    }
  };

  explicit DirectNavigator(std::unique_ptr<const Logger> _logger =
                               getDefaultLogger("DirectNavigator",
                                                Logging::INFO))
      : m_logger{std::move(_logger)} {}

  State makeState(const Options& options) const {
    State state(options);
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

  bool endOfWorldReached(State& /*state*/) const { return false; }

  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

  void initialize(State& state, const Vector3& /*position*/,
                  const Vector3& /*direction*/,
                  Direction propagationDirection) const {
    ACTS_VERBOSE("Initialize. Surface sequence for navigation:");
    for (const Surface* surface : state.options.surfaces) {
      ACTS_VERBOSE(surface->geometryId()
                   << " - "
                   << surface->center(state.options.geoContext).transpose());
    }

    state.direction = propagationDirection;

    // We set the current surface to the start surface
    state.currentSurface = state.options.startSurface;
    if (state.currentSurface != nullptr) {
      ACTS_VERBOSE("Current surface set to start surface "
                   << state.currentSurface->geometryId());
    } else {
      ACTS_VERBOSE("Current surface set to nullptr");
    }

    // Reset the surface index
    state.resetSurfaceIndex();
    bool foundStartSurface = false;
    for (const Surface* surface : state.options.surfaces) {
      if (surface == state.currentSurface) {
        foundStartSurface = true;
        break;
      }
      state.nextSurface();
    }
    ACTS_VERBOSE("Initial surface index set to " << state.surfaceIndex);
    if (!foundStartSurface) {
      ACTS_DEBUG(
          "Did not find the start surface in the sequence. Assuming it is not "
          "part of the sequence. Trusting the correctness of the input "
          "sequence. Resetting the surface index.");
      state.resetSurfaceIndex();
    }

    state.navigationBreak = false;
  }

  NavigationTarget estimateNextTarget(State& state, const Vector3& position,
                                      const Vector3& direction) const {
    if (state.navigationBreak) {
      return NavigationTarget::invalid();
    }

    ACTS_VERBOSE("DirectNavigator::estimateNextTarget");

    // Navigator target always resets the current surface
    state.currentSurface = nullptr;

    // Move the sequence to the next surface
    state.nextSurface();

    if (!state.endOfSurfaces()) {
      ACTS_VERBOSE("Next surface candidate is  "
                   << state.navSurface().geometryId() << ". "
                   << state.remainingSurfaces() << " out of "
                   << state.options.surfaces.size()
                   << " surfaces remain to try.");
    } else {
      ACTS_VERBOSE("End of surfaces reached, navigation break.");
      state.navigationBreak = true;
      return NavigationTarget::invalid();
    }

    // Establish & update the surface status
    // TODO we do not know the intersection index - passing the closer one
    const Surface& surface = state.navSurface();
    const double farLimit = std::numeric_limits<double>::max();
    const auto intersection = chooseIntersection(
        state.options.geoContext, surface, position, direction,
        BoundaryTolerance::Infinite(), state.options.nearLimit, farLimit,
        state.options.surfaceTolerance);
    return NavigationTarget(surface, intersection.index(),
                            BoundaryTolerance::Infinite());
  }

  bool checkTargetValid(const State& /*state*/, const Vector3& /*position*/,
                        const Vector3& /*direction*/) const {
    return true;
  }

  void handleSurfaceReached(State& state, const Vector3& /*position*/,
                            const Vector3& /*direction*/,
                            const Surface& /*surface*/) const {
    if (state.navigationBreak) {
      return;
    }

    ACTS_VERBOSE("DirectNavigator::handleSurfaceReached");

    // Set the current surface
    state.currentSurface = &state.navSurface();
    ACTS_VERBOSE("Current surface set to  "
                 << state.currentSurface->geometryId());
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
