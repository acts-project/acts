// This file is part of the ACTS project.
//
// Copyright (C) 2026 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Result.hpp"

#include <concepts>

namespace Acts::detail {

/// Calculate geometric curvature from the current state and local field.
/// @param stepper Stepper providing position, direction, charge and field
/// @param state Current stepping state; its field cache may be updated
/// @return Global curvature in inverse length, or a field lookup error
template <typename Stepper>
Result<Vector3> magneticCurvature(const Stepper& stepper,
                                  typename Stepper::State& state) {
  // Neutral parameters may store 1/p in the q/p slot.
  if (stepper.charge(state) == 0.) {
    return Result<Vector3>::success(Vector3::Zero());
  }
  const auto field = stepper.getField(state, stepper.position(state));
  if (!field.ok()) {
    return Result<Vector3>::failure(field.error());
  }
  return Result<Vector3>::success(stepper.qOverP(state) *
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

}  // namespace Acts::detail
