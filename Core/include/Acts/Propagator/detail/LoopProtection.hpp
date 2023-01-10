// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {
namespace detail {

/// Estimate the loop protection limit
template <typename path_aborter_t, typename propagator_state_t,
          typename stepper_t>
void setupLoopProtection(propagator_state_t& state, const stepper_t& stepper,
                         path_aborter_t& pathAborter, const Logger& logger) {
  if (!state.options.loopProtection) {
    return;
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
  Vector3 field = *fieldRes;
  const double B = field.norm();
  if (B == 0) {
    return;
  }

  // Transverse component at start is taken for the loop protection
  const double p = stepper.momentum(state.stepping);
  // Calculate the full helix path
  const double helixPath = state.stepping.navDir * 2 * M_PI * p / B;
  // And set it as the loop limit if it overwrites the internal limit
  double loopLimit = state.options.loopFraction * helixPath;
  double pathLimit = pathAborter.internalLimit;
  if (std::abs(loopLimit) < std::abs(pathLimit)) {
    pathAborter.internalLimit = loopLimit;

    ACTS_VERBOSE("Path aborter limit set to "
                 << loopLimit << " (full helix =  " << helixPath << ")");
  }
}

}  // namespace detail
}  // namespace Acts
