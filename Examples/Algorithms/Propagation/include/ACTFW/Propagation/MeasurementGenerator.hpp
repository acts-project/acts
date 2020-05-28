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
struct MeasurementGenerator {
  /// Simple result struct to be returned
  /// It collects the generated measurements
  struct this_result {
    std::vector<Vector3D> spacepoints = {};
    std::vector<double> rts = {};
    std::vector<double>::iterator rt = rts.end();
  };

  using result_type = this_result;

  std::vector<double> radii = {};

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
    if (result.rts.empty() and result.rt == result.rts.end()) {
      result.rts = radii;
      result.rt = result.rts.begin();
    } else if (result.rt == result.rts.end()) {
      return;
    }

    // Get the position
    auto position = stepper.position(state.stepping);
    double r = VectorHelpers::perp(position);
    // Check if you are on a radius
    double R = (*result.rt);
    auto dr = R - r;
    if (std::abs(dr) < s_onSurfaceTolerance) {
      result.spacepoints.push_back(position);
      ++result.rt;
      return;
    }

    // Calculate and Adapt the step size
    auto direction = stepper.direction(state.stepping);

    // Get the transformation matrtix
    Vector3D caxis(0., 0., 1);

    // Check documentation for explanation
    Vector3D pcXcd = position.cross(caxis);
    Vector3D ldXcd = direction.cross(caxis);
    double a = ldXcd.dot(ldXcd);
    double b = 2. * (ldXcd.dot(pcXcd));
    double c = pcXcd.dot(pcXcd) - (R * R);
    // And solve the qaudratic equation
    auto rqe = detail::RealQuadraticEquation(a, b, c);

    double path = 0.;
    if (result.rt == result.rts.begin()) {
      // always take the positive one for the first step
      path = rqe.first < 0 ? rqe.second : rqe.first;
    } else {
      // take the closest in case of overstepping
      path = rqe.first * rqe.first < rqe.second * rqe.second ? rqe.first
                                                             : rqe.second;
    }

    Vector3D solution = position + path * direction;
    double dS = std::copysign((solution - position).norm(), path);
    stepper.setStepSize(state.stepping, dS);
    return;
  }

  /// Pure observer interface
  /// - this does not apply to the output collector
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& /*state*/,
                  const stepper_t& /*unused*/) const {}
};

}  // namespace Acts