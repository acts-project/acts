// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/AnyIntersection.hpp"
#include "Acts/Propagator/detail/NavigationCandidate.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <limits>
#include <memory>
#include <variant>
#include <vector>

namespace Acts {

class TryAllOverstepNavigator {
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
    BoundaryCheck boundaryCheckSurfaceApproach = BoundaryCheck(true);
  };

  struct IntersectionCandidate {
    detail::AnyIntersection intersection;
    BoundaryCheck boundaryCheck;

    IntersectionCandidate(detail::AnyIntersection _intersection,
                          BoundaryCheck _boundaryCheck)
        : intersection(std::move(_intersection)),
          boundaryCheck(std::move(_boundaryCheck)) {}

    static bool forwardOrder(const IntersectionCandidate& aCandidate,
                             const IntersectionCandidate& bCandidate) {
      return Intersection3D::forwardOrder(
          aCandidate.intersection.intersection(),
          bCandidate.intersection.intersection());
    }
  };

  struct State {
    const Surface* startSurface = nullptr;
    const TrackingVolume* startVolume = nullptr;

    const Surface* targetSurface = nullptr;

    const Surface* currentSurface = nullptr;
    const TrackingVolume* currentVolume = nullptr;

    std::vector<detail::NavigationCandidate> candidates;
    std::vector<IntersectionCandidate> activeCandidates;
    std::size_t activeCandidateIndex = 0;

    std::optional<Vector3> lastPosition;
    std::optional<detail::AnyIntersection> lastIntersection;

    bool targetReached = false;
    bool navigationBreak = false;

    const IntersectionCandidate& activeCandidate() const {
      return activeCandidates.at(activeCandidateIndex);
    }
  };

  TryAllOverstepNavigator(Config cfg,
                          std::unique_ptr<const Logger> _logger =
                              getDefaultLogger("TryAllOverstepNavigator",
                                               Logging::INFO))
      : m_cfg(std::move(cfg)), m_logger{std::move(_logger)} {}

  State makeState(const Surface* startSurface,
                  const Surface* targetSurface) const {
    State result;
    result.startSurface = startSurface;
    result.targetSurface = targetSurface;
    return result;
  }

  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  const TrackingVolume* currentVolume(const State& state) const {
    return state.currentVolume;
  }

  const IVolumeMaterial* currentVolumeMaterial(const State& state) const {
    if (state.currentVolume == nullptr) {
      return nullptr;
    }
    return state.currentVolume->volumeMaterial();
  }

  const Surface* startSurface(const State& state) const {
    return state.startSurface;
  }

  const Surface* targetSurface(const State& state) const {
    return state.targetSurface;
  }

  bool targetReached(const State& state) const { return state.targetReached; }

  bool endOfWorldReached(State& state) const {
    return state.currentVolume == nullptr;
  }

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

  template <typename propagator_state_t, typename stepper_t>
  void initialize(propagator_state_t& state, const stepper_t& stepper) const {
    ACTS_VERBOSE("initialize");

    if (state.navigation.startSurface != nullptr &&
        state.navigation.startSurface->associatedLayer() != nullptr) {
      ACTS_VERBOSE(
          "Fast start initialization through association from Surface.");
      const auto* startLayer = state.navigation.startSurface->associatedLayer();
      state.navigation.startVolume = startLayer->trackingVolume();
    } else if (state.navigation.startVolume != nullptr) {
      ACTS_VERBOSE(
          "Fast start initialization through association from Volume.");
    } else {
      ACTS_VERBOSE("Slow start initialization through search.");
      ACTS_VERBOSE("Starting from position "
                   << toString(stepper.position(state.stepping))
                   << " and direction "
                   << toString(stepper.direction(state.stepping)));
      state.navigation.startVolume =
          m_cfg.trackingGeometry->lowestTrackingVolume(
              state.geoContext, stepper.position(state.stepping));
    }

    // Initialize current volume, layer and surface
    {
      state.navigation.currentVolume = state.navigation.startVolume;
      if (state.navigation.currentVolume != nullptr) {
        ACTS_VERBOSE(volInfo(state) << "Start volume resolved.");
      } else {
        ACTS_ERROR("Start volume not resolved.");
      }

      state.navigation.currentSurface = state.navigation.startSurface;
      if (state.navigation.currentSurface != nullptr) {
        ACTS_VERBOSE(volInfo(state)
                     << "Current surface set to start surface "
                     << state.navigation.currentSurface->geometryId());

        state.navigation.lastIntersection = detail::AnyIntersection(
            state.navigation.currentSurface
                ->intersect(
                    state.geoContext, stepper.position(state.stepping),
                    state.options.direction * stepper.direction(state.stepping),
                    BoundaryCheck(false), state.options.surfaceTolerance)
                .closest());
      } else {
        ACTS_VERBOSE(volInfo(state) << "No start surface set.");
      }
    }

    // Initialize navigation candidates for the start volume
    reinitializeCandidates(state);

    state.navigation.lastPosition.reset();
    state.navigation.lastIntersection.reset();
  }

  template <typename propagator_state_t, typename stepper_t>
  void preStep(propagator_state_t& state, const stepper_t& stepper) const {
    ACTS_VERBOSE(volInfo(state) << "pre step");

    // Navigator preStep always resets the current surface
    state.navigation.currentSurface = nullptr;

    ACTS_VERBOSE(volInfo(state) << "handle active candidates");

    // Check next navigation candidate
    while (state.navigation.activeCandidateIndex !=
           state.navigation.activeCandidates.size()) {
      // Screen output how much is left to try
      ACTS_VERBOSE(volInfo(state)
                   << (state.navigation.activeCandidates.size() -
                       state.navigation.activeCandidateIndex)
                   << " out of " << state.navigation.activeCandidates.size()
                   << " surfaces remain to try.");

      const auto& candidate = state.navigation.activeCandidate();
      const auto& intersection = candidate.intersection;
      const Surface& surface = *intersection.representation();
      BoundaryCheck boundaryCheck = candidate.boundaryCheck;

      // Screen output which surface you are on
      ACTS_VERBOSE(volInfo(state) << "Next surface candidate will be "
                                  << surface.geometryId());

      // Estimate the surface status
      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, surface, intersection.index(),
          state.options.direction, boundaryCheck,
          state.options.surfaceTolerance, logger());

      if (surfaceStatus == IntersectionStatus::onSurface) {
        ACTS_ERROR(volInfo(state) << "We are on surface before trying to reach "
                                     "it. This should not happen. Good luck.");
        ++state.navigation.activeCandidateIndex;
        continue;
      }

      if (surfaceStatus == IntersectionStatus::reachable) {
        ACTS_VERBOSE(volInfo(state)
                     << "Surface reachable, step size updated to "
                     << stepper.outputStepSize(state.stepping));
        break;
      }

      ACTS_WARNING(volInfo(state) << "Surface unreachable, skip.");
      ++state.navigation.activeCandidateIndex;
    }

    if (state.navigation.activeCandidateIndex ==
        state.navigation.activeCandidates.size()) {
      state.navigation.lastPosition = stepper.position(state.stepping);

      stepper.releaseStepSize(state.stepping);

      ACTS_VERBOSE(volInfo(state)
                   << "blindly step forwards. step size updated to "
                   << stepper.outputStepSize(state.stepping));

      return;
    }
  }

  template <typename propagator_state_t, typename stepper_t>
  bool postStep(propagator_state_t& state, const stepper_t& stepper) const {
    ACTS_VERBOSE(volInfo(state) << "post step");

    assert(state.navigation.currentSurface == nullptr &&
           "Current surface must be reset.");

    bool result = true;

    if (state.navigation.activeCandidateIndex ==
        state.navigation.activeCandidates.size()) {
      ACTS_VERBOSE(volInfo(state) << "evaluate blind step");

      state.navigation.activeCandidates.clear();

      assert(state.navigation.lastPosition.has_value() &&
             "last position not set");

      Vector3 start = state.navigation.lastPosition.value();
      Vector3 end = stepper.position(state.stepping);
      Vector3 step = end - start;
      double distance = step.norm();
      Vector3 direction = step.normalized();

      double nearLimit = -distance + state.options.surfaceTolerance;
      double farLimit = state.options.surfaceTolerance;

      // Find intersections with all candidates
      for (const auto& candidate : state.navigation.candidates) {
        auto intersections = candidate.intersect(
            state.geoContext, end, direction, state.options.surfaceTolerance);
        for (const auto& intersection : intersections.split()) {
          // exclude invalid intersections
          if (!intersection) {
            continue;
          }
          // exclude intersection outside of step
          if (intersection.pathLength() < nearLimit ||
              intersection.pathLength() > farLimit) {
            continue;
          }
          // exclude last candidate
          if (state.navigation.lastIntersection.has_value() &&
              state.navigation.lastIntersection->representation() ==
                  intersection.representation() &&
              state.navigation.lastIntersection->index() ==
                  intersection.index()) {
            continue;
          }
          // store candidate
          state.navigation.activeCandidates.emplace_back(
              intersection, candidate.boundaryCheck);
        }
      }

      std::sort(state.navigation.activeCandidates.begin(),
                state.navigation.activeCandidates.end(),
                IntersectionCandidate::forwardOrder);

      state.navigation.activeCandidateIndex = 0;

      ACTS_VERBOSE(volInfo(state)
                   << "found " << state.navigation.activeCandidates.size()
                   << " intersections")

      result = state.navigation.activeCandidateIndex ==
               state.navigation.activeCandidates.size();
    }

    if (state.navigation.activeCandidateIndex !=
        state.navigation.activeCandidates.size()) {
      ACTS_VERBOSE(volInfo(state) << "handle active candidates");

      std::vector<IntersectionCandidate> hitCandidates;

      while (state.navigation.activeCandidateIndex !=
             state.navigation.activeCandidates.size()) {
        const auto& candidate = state.navigation.activeCandidate();
        const auto& intersection = candidate.intersection;
        const Surface& surface = *intersection.representation();

        Intersection3D::Status surfaceStatus = stepper.updateSurfaceStatus(
            state.stepping, surface, intersection.index(),
            state.options.direction, BoundaryCheck(false),
            state.options.surfaceTolerance, logger());

        if (surfaceStatus != IntersectionStatus::onSurface) {
          break;
        }

        hitCandidates.emplace_back(candidate);

        ++state.navigation.activeCandidateIndex;
      }

      if (hitCandidates.empty()) {
        ACTS_VERBOSE(volInfo(state) << "Staying focussed on surface.");
        return result;
      }

      result = true;
      state.navigation.lastIntersection.reset();

      std::vector<IntersectionCandidate> trueHitCandidates;

      for (const auto& candidate : hitCandidates) {
        const auto& intersection = candidate.intersection;
        const Surface& surface = *intersection.representation();

        Intersection3D::Status surfaceStatus = stepper.updateSurfaceStatus(
            state.stepping, surface, intersection.index(),
            state.options.direction, BoundaryCheck(true),
            state.options.surfaceTolerance, logger());

        if (surfaceStatus != IntersectionStatus::onSurface) {
          continue;
        }

        trueHitCandidates.emplace_back(candidate);
      }

      ACTS_VERBOSE(volInfo(state)
                   << "Found " << trueHitCandidates.size()
                   << " intersections on surface with bounds check.");

      if (trueHitCandidates.empty()) {
        ACTS_VERBOSE(volInfo(state)
                     << "Surface successfully hit, but outside bounds.");
        return true;
      }

      if (trueHitCandidates.size() > 1) {
        ACTS_VERBOSE(volInfo(state)
                     << "Only using first intersection within bounds.");
      }

      const auto& intersection = trueHitCandidates.front().intersection;
      const Surface& surface = *intersection.representation();

      state.navigation.lastIntersection = intersection;

      ACTS_VERBOSE(volInfo(state) << "Surface successfully hit, storing it.");
      // Set in navigation state, so actors and aborters can access it
      state.navigation.currentSurface = &surface;
      if (state.navigation.currentSurface) {
        ACTS_VERBOSE(volInfo(state) << "Current surface set to surface "
                                    << surface.geometryId());
      }

      if (intersection.template checkType<SurfaceIntersection>()) {
        ACTS_VERBOSE(volInfo(state) << "This is a surface");
      } else if (intersection.template checkType<LayerIntersection>()) {
        ACTS_VERBOSE(volInfo(state) << "This is a layer");
      } else if (intersection.template checkType<BoundaryIntersection>()) {
        ACTS_VERBOSE(volInfo(state)
                     << "This is a boundary. Reinitialize navigation");

        const auto& boundary =
            *intersection.template object<BoundaryIntersection>();

        state.navigation.currentVolume = boundary.attachedVolume(
            state.geoContext, stepper.position(state.stepping),
            state.options.direction * stepper.direction(state.stepping));

        ACTS_VERBOSE(volInfo(state) << "Switched volume");

        reinitializeCandidates(state);
      } else {
        ACTS_ERROR(volInfo(state) << "Unknown intersection type");
      }
    }

    return result;
  }

 private:
  template <typename propagator_state_t>
  void initializeVolumeCandidates(propagator_state_t& state) const {
    const TrackingVolume* volume = state.navigation.currentVolume;
    ACTS_VERBOSE(volInfo(state) << "Initialize volume");

    if (volume == nullptr) {
      state.navigation.navigationBreak = true;
      ACTS_VERBOSE(volInfo(state) << "No volume set. Good luck.");
      return;
    }

    auto addCandidate = [&](detail::NavigationCandidate::AnyObject object,
                            const Surface* representation,
                            BoundaryCheck boundaryCheck) {
      state.navigation.candidates.emplace_back(object, representation,
                                               boundaryCheck);
    };

    // Find boundary candidates
    {
      ACTS_VERBOSE(volInfo(state) << "Searching for boundaries.");

      const auto& boundaries = volume->boundarySurfaces();

      ACTS_VERBOSE(volInfo(state)
                   << "Found " << boundaries.size() << " boundaries.");

      for (const auto& boundary : boundaries) {
        addCandidate(boundary.get(), &boundary->surfaceRepresentation(),
                     BoundaryCheck(true));
      }
    }

    // Find layer candidates
    {
      ACTS_VERBOSE(volInfo(state) << "Searching for layers.");

      const auto& layers = volume->confinedLayers()->arrayObjects();

      ACTS_VERBOSE(volInfo(state) << "Found " << layers.size() << " layers.");

      for (const auto& layer : layers) {
        if (!layer->resolve(m_cfg.resolveSensitive, m_cfg.resolveMaterial,
                            m_cfg.resolvePassive)) {
          continue;
        }

        addCandidate(layer.get(), &layer->surfaceRepresentation(),
                     BoundaryCheck(true));

        if (layer->approachDescriptor() != nullptr) {
          const auto& approaches =
              layer->approachDescriptor()->containedSurfaces();

          for (const auto& approach : approaches) {
            addCandidate(layer.get(), approach, BoundaryCheck(true));
          }
        }

        // Find surface candidates
        {
          ACTS_VERBOSE(volInfo(state) << "Searching for surfaces.");

          if (layer->surfaceArray() != nullptr) {
            const auto& surfaces = layer->surfaceArray()->surfaces();

            ACTS_VERBOSE(volInfo(state)
                         << "Found " << surfaces.size() << " surfaces.");

            for (const auto& surface : surfaces) {
              addCandidate(surface, surface,
                           m_cfg.boundaryCheckSurfaceApproach);
            }
          }
        }
      }
    }
  }

  template <typename propagator_state_t>
  void reinitializeCandidates(propagator_state_t& state) const {
    state.navigation.candidates.clear();
    state.navigation.activeCandidates.clear();
    state.navigation.activeCandidateIndex = 0;

    initializeVolumeCandidates(state);
  }

  template <typename propagator_state_t>
  std::string volInfo(const propagator_state_t& state) const {
    return (state.navigation.currentVolume != nullptr
                ? state.navigation.currentVolume->volumeName()
                : "No Volume") +
           " | ";
  }

  const Logger& logger() const { return *m_logger; }

  Config m_cfg;

  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts
