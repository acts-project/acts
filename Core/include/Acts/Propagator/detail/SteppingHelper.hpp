// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"

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
    const Surface& surface, const BoundaryCheck& bcheck, LoggerWrapper logger) {
  ACTS_VERBOSE("Update single surface status for surface: " << &surface);

  auto sIntersection =
      surface.intersect(state.geoContext, stepper.position(state),
                        state.navDir * stepper.direction(state), bcheck);

  // The intersection is on surface already
  if (sIntersection.intersection.status == Intersection3D::Status::onSurface) {
    // Release navigation step size
    state.stepSize.release(ConstrainedStep::actor);
    ACTS_VERBOSE("Intersection: state is ON SURFACE");
    return Intersection3D::Status::onSurface;
  } else if (sIntersection.intersection or sIntersection.alternative) {
    // Path and overstep limit checking
    double pLimit = state.stepSize.value(ConstrainedStep::aborter);
    double oLimit = stepper.overstepLimit(state);
    auto checkIntersection = [&](const Intersection3D& intersection) -> bool {
      double cLimit = intersection.pathLength;
      ACTS_VERBOSE(" -> pLimit, oLimit, cLimit: " << pLimit << ", " << oLimit
                                                  << ", " << cLimit);
      bool accept = (cLimit > oLimit and cLimit * cLimit < pLimit * pLimit);
      if (accept) {
        ACTS_VERBOSE("Intersection is WITHIN limit");
        stepper.setStepSize(state, state.navDir * cLimit);
      }

      else {
        ACTS_VERBOSE("Intersection is OUTSIDE limit because: ");
        if (cLimit <= oLimit) {
          ACTS_VERBOSE("- intersection path length "
                       << cLimit << " <= overstep limit " << oLimit);
        }
        if (cLimit * cLimit > pLimit * pLimit + s_onSurfaceTolerance) {
          ACTS_VERBOSE("- intersection path length "
                       << std::abs(cLimit) << " is over the path limit "
                       << (std::abs(pLimit) + s_onSurfaceTolerance)
                       << " (including tolerance of "
                       << s_curvilinearProjTolerance << ")");
        }
      }

      return accept;
    };
    // If either of the two intersections are viable return reachable
    if (checkIntersection(sIntersection.intersection) or
        (sIntersection.alternative and
         checkIntersection(sIntersection.alternative))) {
      ACTS_VERBOSE("Surface is reachable");
      return Intersection3D::Status::reachable;
    }
  }
  ACTS_VERBOSE("Surface is NOT reachable");
  return Intersection3D::Status::unreachable;
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
void updateSingleStepSize(typename stepper_t::State& state,
                          const object_intersection_t& oIntersection,
                          bool release = true) {
  double stepSize = oIntersection.intersection.pathLength;
  state.stepSize.update(stepSize, ConstrainedStep::actor, release);
}

}  // namespace detail
}  // namespace Acts