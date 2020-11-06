// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {

namespace detail {

/// Update surface status - Single component
///
/// This method intersect the provided surface and update the navigation
/// step estimation accordingly (hence it changes the state). It also
/// returns the status of the intersection to trigger onSurface in case
/// the surface is reached.
///
/// @param state [in,out] The stepping state (thread-local cache)
/// @param surface [in] The surface provided
/// @param bcheck [in] The boundary check for this status update
template <typename stepper_t>
Acts::Intersection3D::Status updateSingleSurfaceStatus(
    const stepper_t& stepper, typename stepper_t::State& state,
    const Surface& surface, const BoundaryCheck& bcheck) {
  auto sIntersection = surface.intersect(
      stepper.geometryContext(state), stepper.position(state),
      stepper.steppingDirection(state) * stepper.direction(state), bcheck);

  // The intersection is on surface already
  if (sIntersection.intersection.status == Intersection3D::Status::onSurface) {
    // Release navigation step size
    stepper.stepControl.release(state, ConstrainedStep::actor);
    return Intersection3D::Status::onSurface;
  } else if (sIntersection.intersection or sIntersection.alternative) {
    // Path and overstep limit checking
    double pLimit = stepper.stepControl.size(state, ConstrainedStep::aborter);
    double oLimit = stepper.overstepLimit(state);
    auto checkIntersection = [&](const Intersection3D& intersection) -> bool {
      double cLimit = intersection.pathLength;
      bool accept = (cLimit > oLimit and cLimit * cLimit < pLimit * pLimit);
      if (accept) {
        stepper.stepControl.set(state,
                                stepper.steppingDirection(state) * cLimit);
      }
      return accept;
    };
    // If either of the two intersections are viable return reachable
    if (checkIntersection(sIntersection.intersection) or
        (sIntersection.alternative and
         checkIntersection(sIntersection.alternative))) {
      return Intersection3D::Status::reachable;
    }
  }
  return Intersection3D::Status::unreachable;
}

}  // namespace detail
}  // namespace Acts