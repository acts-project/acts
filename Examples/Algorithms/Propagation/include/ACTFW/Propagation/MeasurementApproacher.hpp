// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// This is a simple3D measurement generator
struct MeasurementApproacher {
  /// Simple result struct to be returned
  /// It collects the generated measurements
  struct this_result {
    std::vector<Vector3D> approaches = {};
    std::vector<Vector3D> mts = {};
    std::vector<Vector3D>::iterator mt = mts.end();
  };

  using result_type = this_result;

  std::vector<Vector3D> measurements = {};

  /// Actor that remembers the volumes passed
  ///
  /// @tparam propagator_state_t is the type of the Propagator state
  /// it is not used in this stepper
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param state is the mutable propagator state object
  /// @param result is the mutable result state object
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    // Initialize the radii if necessary
    if (result.mts.empty() and result.mt == result.mts.end()) {
      result.mts = measurements;
      result.mt = result.mts.begin();
    } else if (result.mt == result.mts.end()) {
      return;
    }

    auto measurement = (*result.mt);

    // Get the position & direction
    auto position = stepper.position(state.stepping);
    auto direction = stepper.direction(state.stepping);

    // Let's create a hyperplane at the measurement point
    using Plane = Eigen::Hyperplane<double, 3>;
    using Line = Eigen::ParametrizedLine<double, 3>;

    // Make the plane and intersect
    auto plane = Plane(direction, measurement);
    auto line = Line::Through(position, position + direction);
    double path = line.intersection(plane);
    if (std::abs(path) < s_onSurfaceTolerance) {
      result.approaches.push_back(position);
      ++result.mt;
      return;
    }

    stepper.setStepSize(state.stepping, path);
    return;
  }

  /// Pure observer interface
  /// - this does not apply to the output collector
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& /*state*/,
                  const stepper_t& /*unused*/) const {}
};

}  // namespace Acts