// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

/// @brief Actor as configurator of the Stepper for working with material
/// interaction. It sets up a tracking of the material and the particles mass
/// and adds further configuration properties.
struct StepActor
{
  // Configurations for Stepper
  /// Toggle between mean and mode evaluation of energy loss
  bool m_meanEnergyLoss = true;
  /// Tolerance for the error of the integration
  double m_tolerance = 5e-5;
  /// Boolean flag for inclusion of d(dEds)d(q/p) into energy loss
  bool m_includeGgradient = true;
  /// Cut-off value for the momentum in SI units
  double m_momentumCutOff = 0.;
  /// Cut-off value for the step size
  double m_stepSizeCutOff = 0.;

  /// @brief Main call operator for setting up stepper properties
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param [in, out] state State of the propagator
  template <typename propagator_state_t>
  void
  operator()(propagator_state_t& state) const
  {
    // If we are on target, everything should have been done
    if (state.navigation.targetReached) {
      return;
    }

    // Initialize all parameters
    if (state.stepping.pathAccumulated == 0.) {
      // Let the stepper track the volume and particles mass
      state.stepping.extension.volume = &state.navigation.currentVolume;
      state.stepping.extension.mass   = state.options.mass;
      state.stepping.extension.pdg    = state.options.absPdgCode;

      // Initialize user defined parameters
      state.stepping.extension.meanEnergyLoss   = m_meanEnergyLoss;
      state.stepping.extension.includeGgradient = m_includeGgradient;
      state.stepping.extension.momentumCutOff   = m_momentumCutOff;
      state.stepping.tolerance                  = m_tolerance;
      state.stepping.stepSizeCutOff             = m_stepSizeCutOff;
    }
  }
};
}
