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
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/detail/CurvatureIntersection.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <concepts>
#include <limits>

namespace Acts::detail {

/// Calculate geometric curvature from the current state and local field.
/// @param stepper Stepper providing position, direction, q/p and field
/// @param state Current stepping state; its field cache may be updated
/// @return Global curvature in inverse length, or a field lookup error
template <typename Stepper>
Result<Vector3> magneticCurvature(const Stepper& stepper,
                                  typename Stepper::State& state) {
  const double qOverP = stepper.qOverP(state);
  // Zero q/p implies no bending and does not require a field lookup.
  if (qOverP == 0.) {
    return Result<Vector3>::success(Vector3::Zero());
  }
  const auto field = stepper.getField(state, stepper.position(state));
  if (!field.ok()) {
    return Result<Vector3>::failure(field.error());
  }
  return Result<Vector3>::success(qOverP *
                                  stepper.direction(state).cross(*field));
}

/// Prefer a stepper's curvature accessor, with compatibility for older
/// steppers.
/// @param stepper Stepper to query
/// @param state Current stepping state
/// @return Global curvature in inverse length, or a field lookup error
template <typename Stepper>
Result<Vector3> getCurvature(const Stepper& stepper,
                             typename Stepper::State& state) {
  if constexpr (requires {
                  { stepper.curvature(state) } -> std::same_as<Result<Vector3>>;
                }) {
    return stepper.curvature(state);
  } else {
    return magneticCurvature(stepper, state);
  }
}

/// Update surface status - Single component
///
/// This method intersect the provided surface and update the navigation
/// step estimation accordingly (hence it changes the state). It also
/// returns the status of the intersection to trigger onSurface in case
/// the surface is reached.
///
/// @tparam stepper_t The type of stepper used for the propagation
///
/// @param stepper [in] The stepper in use
/// @param state [in,out] The stepping state (thread-local cache)
/// @param surface [in] The surface provided
/// @param index [in] The surface intersection index
/// @param direction [in] The propagation direction
/// @param boundaryTolerance [in] The boundary check for this status update
/// @param surfaceTolerance [in] Surface tolerance used for intersection
/// @param stype [in] The step size type to be set
/// @param logger [in] A @c Logger instance
/// @param curvature Global local curvature for the crossing estimate
template <typename stepper_t>
IntersectionStatus updateSingleSurfaceStatus(
    const stepper_t& stepper, typename stepper_t::State& state,
    const Surface& surface, std::uint8_t index, Direction direction,
    const BoundaryTolerance& boundaryTolerance, double surfaceTolerance,
    ConstrainedStep::Type stype, const Logger& logger,
    const Vector3& curvature = Vector3::Zero()) {
  ACTS_VERBOSE("Update single surface status for surface: "
               << surface.geometryId() << " index " << static_cast<int>(index));

  auto sIntersection =
      detail::intersectWithCurvature(
          surface, state.options.geoContext, stepper.position(state),
          direction * stepper.direction(state), curvature, boundaryTolerance,
          surfaceTolerance)
          .at(index);

  // The intersection is on surface already
  if (sIntersection.status() == IntersectionStatus::onSurface) {
    ACTS_VERBOSE("Intersection: state is ON SURFACE");
    state.stepSize.release(stype);
    stepper.updateStepSize(state, sIntersection.pathLength(), stype);
    return IntersectionStatus::onSurface;
  }

  const double nearLimit = std::numeric_limits<double>::lowest();
  const double farLimit = std::numeric_limits<double>::max();

  if (sIntersection.isValid() &&
      detail::checkPathLength(sIntersection.pathLength(), nearLimit, farLimit,
                              logger)) {
    ACTS_VERBOSE("Surface is reachable");
    stepper.releaseStepSize(state, stype);
    stepper.updateStepSize(state, sIntersection.pathLength(), stype);
    return IntersectionStatus::reachable;
  }

  ACTS_VERBOSE("Surface is NOT reachable");
  return IntersectionStatus::unreachable;
}

/// Update a target with curvature, retaining the old API for custom steppers.
/// Composite steppers query each component; single steppers receive their own
/// current curvature through the extended surface-status method.
/// @param stepper Stepper instance
/// @param state Stepping state
/// @param surface Target surface
/// @param index Ray intersection branch
/// @param direction Propagation direction
/// @param boundaryTolerance Physical surface bounds tolerance
/// @param surfaceTolerance On-surface tolerance
/// @param stype Step constraint type
/// @param logger Logger
/// @return Status or curvature lookup error
template <typename Stepper>
Result<IntersectionStatus> updateSurfaceStatusWithCurvature(
    const Stepper& stepper, typename Stepper::State& state,
    const Surface& surface, std::uint8_t index, Direction direction,
    const BoundaryTolerance& boundaryTolerance, double surfaceTolerance,
    ConstrainedStep::Type stype, const Logger& logger) {
  if constexpr (requires {
                  {
                    stepper.updateSurfaceStatusWithCurvature(
                        state, surface, index, direction, boundaryTolerance,
                        surfaceTolerance, stype, logger)
                  } -> std::same_as<Result<IntersectionStatus>>;
                }) {
    return stepper.updateSurfaceStatusWithCurvature(
        state, surface, index, direction, boundaryTolerance, surfaceTolerance,
        stype, logger);
  } else if constexpr (requires {
                         {
                           stepper.updateSurfaceStatus(
                               state, surface, index, direction,
                               boundaryTolerance, surfaceTolerance, stype,
                               logger, Vector3::Zero())
                         } -> std::same_as<IntersectionStatus>;
                       }) {
    const auto curvature = getCurvature(stepper, state);
    if (!curvature.ok()) {
      return curvature.error();
    }
    return Result<IntersectionStatus>::success(stepper.updateSurfaceStatus(
        state, surface, index, direction, boundaryTolerance, surfaceTolerance,
        stype, logger, *curvature));
  } else {
    return Result<IntersectionStatus>::success(stepper.updateSurfaceStatus(
        state, surface, index, direction, boundaryTolerance, surfaceTolerance,
        stype, logger));
  }
}

}  // namespace Acts::detail
