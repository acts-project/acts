// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

namespace Acts {

/// This actor allows to include the scattering contributions to the covariance
/// matrix in dense environment.
struct DenseMaterialInteractor {
  /// Storage variables to track the contribution
  struct Result {
    /// The path length in the environment
    ActsScalar s = 0.;
    /// The q/p value at the entrance
    ActsScalar qop0 = 0.;
    /// The time at the entrance
    ActsScalar t0 = 0.;
    /// The current volume
    TrackingVolume const* currentVolume = nullptr;
    /// The accumulated path length
    ActsScalar sAccumulated = 0.;
    /// The path length in the environment in terms of X0
    ActsScalar sInX0 = 0.;
    /// The direction at the entrance
    Vector3 dir0 = Vector3::Zero();
    /// The position of the previous step
    Vector3 previousPos = Vector3::Zero();
  };
  using result_type = Result;

  /// Main function for the recording and evaluation of the scattering
  /// contribution
  ///
  /// @tparam propagator_state_t Type of the propagation state
  /// @tparam stepper_t Type of the stepper
  /// @param [in, out] state State of the propagation
  /// @param [in] stepper The stepper class
  /// @param [in, out] result The state for tracking the material contribution
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    // The volume changed
    if (result.currentVolume != state.navigation.currentVolume) {
      // If it was a dense environment, the contributions will be evaluated and
      // applied
      if (result.s != 0.) {
        stepper.covarianceTransport(state.stepping);
        BoundSymMatrix msCovariance = evaluateMultipleScattering(
            state.options.absPdgCode, state.options.mass,
            stepper.charge(state.stepping), stepper.momentum(state.stepping),
            stepper.time(state.stepping), result);
        state.stepping.cov += msCovariance;
      }

      // Set the current volume
      result.currentVolume = state.navigation.currentVolume;

      // Reset the tracker variables
      reset(stepper.charge(state.stepping), stepper.momentum(state.stepping),
            stepper.time(state.stepping), stepper.direction(state.stepping),
            result);
    } else {
      // Same volume with material: accumulate path
      if (result.currentVolume != nullptr &&
          result.currentVolume->volumeMaterial() != nullptr) {
        ActsScalar deltaS =
            state.stepping.pathAccumulated - result.sAccumulated;
        result.s += deltaS;
        result.sInX0 += deltaS / result.currentVolume->volumeMaterial()
                                     ->material(result.previousPos)
                                     .X0();
      }

      // Test for recorded path length
      if (result.s != 0.) {
        // Transport and add the contributions if the target is reached
        if (state.navigation.targetReached) {
          stepper.covarianceTransport(state.stepping);
        }
        // Add the contribution whenever a transport occured
        // TODO: Replace this check, depends on #763
        if (state.stepping.jacTransport == FreeMatrix::Identity()) {
          BoundSymMatrix msCovariance = evaluateMultipleScattering(
              state.options.absPdgCode, state.options.mass,
              stepper.charge(state.stepping), stepper.momentum(state.stepping),
              stepper.time(state.stepping), result);
          state.stepping.cov += msCovariance;
          reset(stepper.charge(state.stepping),
                stepper.momentum(state.stepping), stepper.time(state.stepping),
                stepper.direction(state.stepping), result);
        }
      }
    }
    result.sAccumulated = state.stepping.pathAccumulated;
    result.previousPos = stepper.position(state.stepping);
  }

 private:
  /// Resets the tracking parameters
  ///
  /// @param [in] q Charge of the particle
  /// @param [in] momentum Momentum of the particle
  /// @param [in] time Timestamp of the particle
  /// @param [in] direction Direction of the particle
  /// @param [in, out] result Storage of scattering contribution variables
  void reset(ActsScalar q, ActsScalar momentum, ActsScalar time,
             Vector3 direction, result_type& result) const;

  /// Evaluates the scattering contribution
  ///
  /// @param [in] pdg PDG ID of the particle
  /// @param [in] mass Mass of the particle
  /// @param [in] q Charge of the particle
  /// @param [in] momentum Momentum of the particle
  /// @param [in] time Timestamp of the particle
  /// @param [in, out] result Storage of scattering contribution variables
  ///
  /// @return Matrix containing the additivie contributions to the covariance
  /// matrix
  BoundSymMatrix evaluateMultipleScattering(int pdg, float mass, float q,
                                            ActsScalar momentum,
                                            ActsScalar time,
                                            result_type& result) const;
};
}  // namespace Acts
