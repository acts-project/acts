// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <numbers>

namespace Acts::detail {

/// Estimate the loop protection limit
template <typename path_aborter_t, typename propagator_state_t,
          typename stepper_t>
void setupLoopProtection(propagator_state_t& state, const stepper_t& stepper,
                         path_aborter_t& pathAborter, bool releaseLimit,
                         const Logger& logger) {
  if (!state.options.loopProtection) {
    return;
  }

  if (releaseLimit) {
    pathAborter.internalLimit =
        state.options.direction * std::numeric_limits<long double>::max();
  }

  // Get the field at the start position
  auto fieldRes =
      stepper.getField(state.stepping, stepper.position(state.stepping));
  if (!fieldRes.ok()) {
    // there's no great way to return the error here, so resort to warning
    // and not applying the loop protection in this case
    ACTS_WARNING("Field lookup was unsuccessful, this is very likely an error");
    return;
  }
  const Vector3 field = *fieldRes;
  const long double B = field.norm();
  if (B == 0) {
    return;
  }

  // Transverse component at start is taken for the loop protection
  const long double p = stepper.absoluteMomentum(state.stepping);
  // Calculate the full helix path
  const long double helixPath =
      state.options.direction * 2 * std::numbers::pi * p / B;
  // And set it as the loop limit if it overwrites the internal limit
  const long double loopLimit = state.options.loopFraction * helixPath;
  const long double previousLimit = pathAborter.internalLimit;
  if (std::abs(loopLimit) < std::abs(previousLimit)) {
    pathAborter.internalLimit = loopLimit;

    ACTS_VERBOSE("Path aborter limit set to "
                 << loopLimit << " (full helix = " << helixPath
                 << ", previous limit = " << previousLimit << ")");
  } else {
    ACTS_VERBOSE("Path aborter limit not updated to "
                 << loopLimit << " (full helix = " << helixPath
                 << ", previous limit = " << previousLimit << ")");
  }
}

}  // namespace Acts::detail
