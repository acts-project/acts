// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

template <typename Guider>
class GuidedNavigator {
  Guider m_guider;
  std::shared_ptr<const TrackingGeometry> m_tgeo;

 public:
  struct State {
    std::vector<SurfaceIntersection> navSurfaces;
    std::vector<SurfaceIntersection>::const_iterator navSurfaceIter;
    const Surface *currentSurface;
    const Surface *startSurface;

    bool navigationBreak = false;
    bool targetReached = false;

    const TrackingVolume *currentVolume = nullptr;
    const TrackingVolume *worldVolume = nullptr;
    const Surface *targetSurface = nullptr;
  };

  // TODO Must the Navigator really be default constructable?
  // See Propagator.hpp:247
  GuidedNavigator() { throw_assert(false, "Don't use this constructor!"); }

  /// Constructor of the ML Navigator
  GuidedNavigator(Guider &&guider,
                  std::shared_ptr<const Acts::TrackingGeometry> tgeo)
      : m_guider(guider), m_tgeo(tgeo) {}

  /// The whole purpose of this function is to set currentSurface, if possible
  template <typename propagator_state_t, typename stepper_t>
  void status(propagator_state_t &state, const stepper_t &stepper) const {
    const auto &logger = state.options.logger;
    ACTS_VERBOSE(">>>>>>>> STATUS <<<<<<<<<");

    if (state.navigation.navigationBreak)
      return;

    // TODO Is initialization correct for all cases?
    if (state.navigation.navSurfaces.empty()) {
      initialize(state, stepper);
      return;
    }

    // Navigator status always resets the current surface
    state.navigation.currentSurface = nullptr;

    // Establish the surface status
    auto surfaceStatus = stepper.updateSurfaceStatus(
        state.stepping, *state.navigation.navSurfaceIter->representation,
        false);

    if (surfaceStatus == Acts::Intersection3D::Status::onSurface) {
      // Set the current surface
      state.navigation.currentSurface =
          state.navigation.navSurfaceIter->representation;
      ACTS_VERBOSE("Current surface set to  "
                   << state.navigation.currentSurface->geometryId());

      // Release Stepsize
      ACTS_VERBOSE("Release Stepsize");
      stepper.releaseStepSize(state.stepping);

      // Reset state
      state.navigation.navSurfaces.clear();
      state.navigation.navSurfaceIter = state.navigation.navSurfaces.end();

    } else if (surfaceStatus == Acts::Intersection3D::Status::reachable) {
      ACTS_VERBOSE("Next surface reachable at distance  "
                   << stepper.outputStepSize(state.stepping));
    } else {
      ACTS_VERBOSE(
          "Surface unreachable or missed, hopefully the target(...) call "
          "fixes this");
    }
  }

  template <typename propagator_state_t, typename stepper_t>
  void target(propagator_state_t &state, const stepper_t &stepper) const {
    const auto &logger = state.options.logger;
    ACTS_VERBOSE(">>>>>>>> TARGET <<<<<<<<<");
    auto &navstate = state.navigation;

    // TODO Workaround since we don't have a currentVolume yet
    // In the case of the trial&error surface provider this does not happen,
    // since it stops if no surface is reachable
    if (!navstate.worldVolume->inside(stepper.position(state.stepping))) {
      navstate.navigationBreak = true;
      navstate.currentVolume = nullptr;
      ACTS_VERBOSE("not inside worldVolume, stop navigation");
      return;
    }

    // Fetch new targets if there are no candidates in state
    if (navstate.navSurfaceIter == navstate.navSurfaces.end()) {
      // This means also we are currently on a surface
      assert(navstate.currentSurface != nullptr);
      ACTS_VERBOSE(
          "It seems like we are on a surface and must fetch new targets");

      navstate.navSurfaces = fetch_surfaces(state, stepper);
      navstate.navSurfaceIter = navstate.navSurfaces.begin();
    }

    // It seems like we are done
    if (navstate.navigationBreak) {
      ACTS_VERBOSE("Navigation break");
      return;
    }

    // Check if we are in a correct state
    assert(!state.navigation.navSurfaces.empty());
    assert(state.navigation.navSurfaceIter !=
           state.navigation.navSurfaces.end());

    // Navigator target always resets the current surface
    // It is set later by the status call if possible
    navstate.currentSurface = nullptr;

    ACTS_VERBOSE("Ask for SurfaceStatus of currently most probable target");

    // Establish & update the surface status
    auto surfaceStatus = stepper.updateSurfaceStatus(
        state.stepping, *navstate.navSurfaceIter->representation, false);

    // Everything OK
    if (surfaceStatus == Acts::Intersection3D::Status::reachable) {
      ACTS_VERBOSE("Navigation stepSize set to "
                   << stepper.outputStepSize(state.stepping));
    }
    // Try another surface
    else {
      ACTS_VERBOSE(
          "Surface not reachable, search another one which is "
          "reachable");

      // Search in the current surface pool
      for (; navstate.navSurfaceIter != navstate.navSurfaces.end();
           ++navstate.navSurfaceIter) {
        if (stepper.updateSurfaceStatus(
                state.stepping, *navstate.navSurfaceIter->representation,
                false) == Intersection3D::Status::reachable)
          break;
      }

      // If we did not find a surface, stop the navigation
      if (navstate.navSurfaceIter == navstate.navSurfaces.end()) {
        ACTS_VERBOSE("Could not a reachable surface, stop navigation!");
        navstate.currentVolume = nullptr;
        navstate.navigationBreak = true;
      } else {
        ACTS_VERBOSE("Found a reachable surface "
                     << navstate.navSurfaceIter->representation->geometryId());
      }
    }
  }

  template <typename propagator_state_t, typename stepper_t>
  auto initialize(propagator_state_t &state, const stepper_t &) const {
    state.navigation.currentSurface = state.navigation.startSurface;
    state.navigation.worldVolume = m_tgeo->highestTrackingVolume();

    // TODO this is just a workaround, so it is not nullptr...
    state.navigation.currentVolume = m_tgeo->highestTrackingVolume();
  }

  /// @brief Fetches new surfaces from provider and processes the resulting
  /// std::vector of surfaces
  /// @return Sorted std::vector of SurfaceIntersections
  template <typename propagator_state_t, typename stepper_t>
  auto fetch_surfaces(const propagator_state_t &state,
                      const stepper_t &stepper) const {
    const auto &logger = state.options.logger;

    // Do prediction
    const std::vector<const Surface *> predictions = m_guider(state, stepper);

    // Important parameters
    const FreeVector &params = state.stepping.pars;
    const GeometryContext &gctx = state.geoContext;
    //     const double olimit = stepper.overstepLimit(state.stepping);

    // Make intersection objects
    std::vector<SurfaceIntersection> sfis(predictions.size());
    std::transform(
        predictions.begin(), predictions.end(), sfis.begin(), [&](auto s) {
          return s->intersect(
              gctx, params.segment<3>(eFreePos0),
              state.stepping.navDir * params.segment<3>(eFreeDir0), true);
        });

    ACTS_VERBOSE("have " << sfis.size() << " candidate intersections");

    // TODO overstepLimit etc. missing here
    sfis.erase(std::remove_if(sfis.begin(), sfis.end(),
                              [](auto sfi) { return !static_cast<bool>(sfi); }),
               sfis.end());

    ACTS_VERBOSE("after remove " << sfis.size() << " candidates remain");

    if (state.stepping.navDir == forward) {
      std::sort(sfis.begin(), sfis.end());
    } else {
      std::sort(sfis.begin(), sfis.end(), std::greater<>());
    }

    return sfis;
  }
};

}  // namespace Acts
