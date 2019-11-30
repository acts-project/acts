// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {

namespace detail {

using cstep = detail::ConstrainedStep;

/// Update surface status - Single component
///
/// It checks the status to the reference surface & updates
/// the step size accordingly
///
/// @param state [in,out] The stepping state (thread-local cache)
/// @param surface [in] The surface provided
/// @param bcheck [in] The boundary check for this status update
template <typename stepper_t>
Acts::Intersection::Status updateSurfaceStatus_t(
    const stepper_t& stepper, typename stepper_t::State& state,
    const Surface& surface, const BoundaryCheck& bcheck) {
  // Now intersect (should exclude punch-through)
  auto sIntersection =
      surface.intersect(state.geoContext, stepper.position(state),
                        state.navDir * stepper.direction(state), bcheck);

  // The intersection is on surface already
  if (sIntersection.intersection.status == Intersection::Status::onSurface) {
    // Release navigation step size
    state.stepSize.release(cstep::actor);
    return Intersection::Status::onSurface;
  } else if (sIntersection.intersection or sIntersection.alternative) {
    // Path and overstep limit checking
    double pLimit = state.stepSize.value(cstep::aborter);
    double oLimit = stepper.overstepLimit(state);
    auto checkIntersection = [&](const Intersection& intersection) -> bool {
      double cLimit = intersection.pathLength;
      bool accept = (cLimit > oLimit and cLimit * cLimit < pLimit * pLimit);
      if (accept) {
        stepper.setStepSize(state, state.navDir * cLimit);
      }
      return accept;
    };
    // If either of the two intersections are viable return reachable
    if (checkIntersection(sIntersection.intersection) or
        (sIntersection.alternative and
         checkIntersection(sIntersection.alternative))) {
      return Intersection::Status::reachable;
    }
  }
  return Intersection::Status::unreachable;
}

/// Update the Step size - single component
///
/// It takes a (valid) object intersection from the compatibleX(...)
/// calls in the geometry and updates the step size
///
/// @param state [in,out] The stepping state (thread-local cache)
/// @param oIntersection [in] The object that yielded this step size
/// @param release [in] A release flag
template <typename stepper_t, typename object_intersection_t>
void updateStepSize_t(typename stepper_t::State& state,
                      const object_intersection_t& oIntersection,
                      bool release = true) {
  double stepSize = oIntersection.intersection.pathLength;
  state.stepSize.update(stepSize, cstep::actor, release);
}

}  // namespace detail
}  // namespace Acts