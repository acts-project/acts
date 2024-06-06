// This file is part of the Acts project.
//
// Copyright (C) 2023-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/detail/NavigationHelpers.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

namespace Acts {

/// @brief Alternative @c Navigator which tries all possible intersections
///
/// See @c Navigator for more general information about the Navigator concept.
///
/// This Navigator tries all possible intersections with all surfaces in the
/// current volume. It does not use any information about the geometry to
/// optimize the search. It is therefore very slow, but can be used as a
/// reference implementation.
///
class TryAllNavigator {
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

  /// @brief Nested State struct
  ///
  /// It acts as an internal state which is created for every propagation and
  /// meant to keep thread-local navigation information.
  struct State {
    // Starting geometry information of the navigation which should only be set
    // while initialization. NOTE: This information is mostly used by actors to
    // check if we are on the starting surface (e.g. MaterialInteraction).
    const Surface* startSurface = nullptr;
    const TrackingVolume* startVolume = nullptr;

    // Target geometry information of the navigation which should only be set
    // while initialization. NOTE: This information is mostly used by actors to
    // check if we are on the target surface (e.g. MaterialInteraction).
    const Surface* targetSurface = nullptr;

    // Current geometry information of the navigation which is set during
    // initialization and potentially updated after each step.
    const Surface* currentSurface = nullptr;
    const TrackingVolume* currentVolume = nullptr;

    /// The vector of navigation candidates to work through
    std::vector<detail::NavigationObjectCandidate> navigationCandidates;
    /// The vector of intersection candidates to work through
    std::vector<detail::IntersectionCandidate> intersectionCandidates;

    /// Indicator if the target is reached
    bool targetReached = false;
    /// If a break has been detected
    bool navigationBreak = false;
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param _logger a logger instance
  TryAllNavigator(Config cfg,
                  std::unique_ptr<const Logger> _logger =
                      getDefaultLogger("TryAllNavigator", Logging::INFO))
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

  /// @brief Initialize call - start of navigation
  ///
  /// @tparam propagator_state_t The state type of the propagator
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void initialize(propagator_state_t& state, const stepper_t& stepper) const {
    ACTS_VERBOSE("initialize");

    const TrackingVolume* startVolume = nullptr;

    if (state.navigation.startSurface != nullptr &&
        state.navigation.startSurface->associatedLayer() != nullptr) {
      ACTS_VERBOSE(
          "Fast start initialization through association from Surface.");
      const auto* startLayer = state.navigation.startSurface->associatedLayer();
      startVolume = startLayer->trackingVolume();
    } else {
      ACTS_VERBOSE("Slow start initialization through search.");
      ACTS_VERBOSE("Starting from position "
                   << toString(stepper.position(state.stepping))
                   << " and direction "
                   << toString(stepper.direction(state.stepping)));
      startVolume = m_cfg.trackingGeometry->lowestTrackingVolume(
          state.geoContext, stepper.position(state.stepping));
    }

    // Initialize current volume, layer and surface
    {
      state.navigation.currentVolume = startVolume;
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
      } else {
        ACTS_VERBOSE(volInfo(state) << "No start surface set.");
      }
    }

    // Initialize navigation candidates for the start volume
    reinitializeCandidates(state);
  }

  /// @brief Navigator pre step call
  ///
  /// This determines the next surface to be targeted and sets the step length
  /// accordingly.
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void preStep(propagator_state_t& state, const stepper_t& stepper) const {
    ACTS_VERBOSE(volInfo(state) << "pre step");

    // Navigator preStep always resets the current surface
    state.navigation.currentSurface = nullptr;

    ACTS_VERBOSE(volInfo(state) << "intersect candidates");

    Vector3 position = stepper.position(state.stepping);
    Vector3 direction =
        state.options.direction * stepper.direction(state.stepping);

    double nearLimit = state.options.surfaceTolerance;
    double farLimit = std::numeric_limits<double>::max();

    // handle overstepping
    if (!state.navigation.intersectionCandidates.empty()) {
      const detail::IntersectionCandidate& previousIntersection =
          state.navigation.intersectionCandidates.front();

      const Surface& surface = *previousIntersection.intersection.object();
      std::uint8_t index = previousIntersection.intersection.index();
      BoundaryCheck boundaryCheck = previousIntersection.boundaryCheck;

      auto intersection = surface.intersect(
          state.geoContext, position, direction, boundaryCheck,
          state.options.surfaceTolerance)[index];

      if (intersection.pathLength() < 0) {
        ACTS_VERBOSE(volInfo(state) << "handle overstepping");

        nearLimit = std::min(nearLimit, intersection.pathLength() -
                                            state.options.surfaceTolerance);
        farLimit = -state.options.surfaceTolerance;
      }
    }

    std::vector<detail::IntersectionCandidate> intersectionCandidates;

    // Find intersections with all candidates
    for (const auto& candidate : state.navigation.navigationCandidates) {
      auto intersections =
          candidate.intersect(state.geoContext, position, direction,
                              state.options.surfaceTolerance);
      for (const auto& intersection : intersections.first.split()) {
        // exclude invalid intersections
        if (!detail::checkIntersection(intersection, nearLimit, farLimit)) {
          continue;
        }
        // store candidate
        intersectionCandidates.emplace_back(intersection, intersections.second,
                                            candidate.boundaryCheck);
      }
    }

    std::sort(intersectionCandidates.begin(), intersectionCandidates.end(),
              detail::IntersectionCandidate::forwardOrder);

    ACTS_VERBOSE(volInfo(state) << "found " << intersectionCandidates.size()
                                << " intersections");

    for (const auto& candidate : intersectionCandidates) {
      const auto& intersection = candidate.intersection;
      const Surface& surface = *intersection.object();
      BoundaryCheck boundaryCheck = candidate.boundaryCheck;

      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, surface, intersection.index(),
          state.options.direction, boundaryCheck,
          state.options.surfaceTolerance, logger());

      if (surfaceStatus == IntersectionStatus::onSurface) {
        ACTS_ERROR(volInfo(state)
                   << "We are on surface " << surface.geometryId()
                   << " before trying to reach it. This should not happen. "
                      "Good luck.");
        continue;
      }

      ACTS_VERBOSE(volInfo(state) << "aiming at surface "
                                  << surface.geometryId() << ". step size is "
                                  << stepper.outputStepSize(state.stepping));
      break;
    }

    if (intersectionCandidates.empty()) {
      stepper.releaseStepSize(state.stepping, ConstrainedStep::actor);

      ACTS_VERBOSE(volInfo(state) << "no intersections found. advance without "
                                     "constraints. step size is "
                                  << stepper.outputStepSize(state.stepping));
    }

    state.navigation.intersectionCandidates = std::move(intersectionCandidates);
  }

  /// @brief Navigator post step call
  ///
  /// This determines if we hit the next navigation candidate and deals with it
  /// accordingly. It sets the current surface, enters layers and changes
  /// volumes.
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  ///
  /// @return Boolean to indicate if we continue with the actors and
  ///         aborters or if we should target again.
  template <typename propagator_state_t, typename stepper_t>
  void postStep(propagator_state_t& state, const stepper_t& stepper) const {
    ACTS_VERBOSE(volInfo(state) << "post step");

    assert(state.navigation.currentSurface == nullptr &&
           "Current surface must be reset.");

    if (state.navigation.intersectionCandidates.empty()) {
      ACTS_VERBOSE(volInfo(state) << "no intersections.");
      return;
    }

    std::vector<detail::IntersectionCandidate> hitCandidates;

    for (const auto& candidate : state.navigation.intersectionCandidates) {
      const auto& intersection = candidate.intersection;
      const Surface& surface = *intersection.object();

      Intersection3D::Status surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, surface, intersection.index(),
          state.options.direction, BoundaryCheck(false),
          state.options.surfaceTolerance, logger());

      if (surfaceStatus != IntersectionStatus::onSurface) {
        break;
      }

      hitCandidates.emplace_back(candidate);
    }

    if (hitCandidates.empty()) {
      ACTS_VERBOSE(volInfo(state) << "Staying focussed on surface.");
      return;
    }

    state.navigation.intersectionCandidates.clear();

    ACTS_VERBOSE(volInfo(state)
                 << "Found " << hitCandidates.size()
                 << " intersections on surface without bounds check.");

    std::vector<detail::IntersectionCandidate> trueHitCandidates;

    for (const auto& candidate : hitCandidates) {
      const auto& intersection = candidate.intersection;
      const Surface& surface = *intersection.object();

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
      return;
    }

    if (trueHitCandidates.size() > 1) {
      ACTS_VERBOSE(volInfo(state)
                   << "Only using first intersection within bounds.");
    }

    const auto& candidate = trueHitCandidates.front();
    const auto& intersection = candidate.intersection;
    const Surface& surface = *intersection.object();

    ACTS_VERBOSE(volInfo(state) << "Surface successfully hit, storing it.");
    // Set in navigation state, so actors and aborters can access it
    state.navigation.currentSurface = &surface;
    if (state.navigation.currentSurface) {
      ACTS_VERBOSE(volInfo(state) << "Current surface set to surface "
                                  << surface.geometryId());
    }

    if (candidate.template checkType<Surface>()) {
      ACTS_VERBOSE(volInfo(state) << "This is a surface");
    } else if (candidate.template checkType<Layer>()) {
      ACTS_VERBOSE(volInfo(state) << "This is a layer");
    } else if (candidate.template checkType<BoundarySurface>()) {
      ACTS_VERBOSE(volInfo(state)
                   << "This is a boundary. Reinitialize navigation");

      const auto& boundary = *candidate.template object<BoundarySurface>();

      state.navigation.currentVolume = boundary.attachedVolume(
          state.geoContext, stepper.position(state.stepping),
          state.options.direction * stepper.direction(state.stepping));

      ACTS_VERBOSE(volInfo(state) << "Switched volume");

      reinitializeCandidates(state);
    } else {
      ACTS_ERROR(volInfo(state) << "Unknown intersection type");
    }
  }

 private:
  /// Helper method to initialize navigation candidates for the current volume.
  template <typename propagator_state_t>
  void initializeVolumeCandidates(propagator_state_t& state) const {
    const TrackingVolume* volume = state.navigation.currentVolume;
    ACTS_VERBOSE(volInfo(state) << "Initialize volume");

    if (volume == nullptr) {
      state.navigation.navigationBreak = true;
      ACTS_VERBOSE(volInfo(state) << "No volume set. Good luck.");
      return;
    }

    emplaceAllVolumeCandidates(state.navigation.navigationCandidates, *volume,
                               m_cfg.resolveSensitive, m_cfg.resolveMaterial,
                               m_cfg.resolvePassive,
                               m_cfg.boundaryCheckSurfaceApproach, logger());
  }

  /// Helper method to reset and reinitialize the navigation candidates.
  template <typename propagator_state_t>
  void reinitializeCandidates(propagator_state_t& state) const {
    state.navigation.navigationCandidates.clear();
    state.navigation.intersectionCandidates.clear();

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
