// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once


namespace Acts {
namespace detail {

struct DenseEnvironmentExtension
{
  /// @brief This struct serves as data container to keep track of all
  /// parameters that are related to an energy loss of a particle in matter.
  struct EnergyLossData
  {
    /// Particles momentum at k1
    double initialMomentum = 0.;
    /// Particles mass in SI units
    double massSI = 0.;
    /// Material that will be passed
    std::shared_ptr<const Material> material;
    /// Derivatives dLambda''dlambda at each sub-step point
    std::array<double, 4> dLdl;
    /// q/p at each sub-step
    std::array<double, 4> qop;
    /// Derivatives dPds at each sub-step
    std::array<double, 4> dPds;
    /// Propagation of derivatives of dLambda''dlambda at each sub-step
    std::array<double, 4> jdL;
    /// Derivative d(dEds)d(q/p) evaluated at the initial point
    double dgdqopValue = 0.;
    /// Derivative dEds at the initial point
    double g = 0.;
  };
  
	    /// Mass
    double mass = 0.;

    /// PDG code
    int pdg = 0;

    /// Volume with material that is passed
    TrackingVolume const* const* volume = nullptr;

    /// Boolean flag for energy loss while stepping
    bool energyLossFlag = true;

    /// Toggle between mean and mode evaluation of energy loss
    bool meanEnergyLoss = true;

    /// Tolerance for the error of the integration
    double tolerance = 5e-5;

    /// Boolean flag for inclusion of d(dEds)d(q/p) into energy loss
    bool includeGgradient = true;

    /// Cut-off value for the momentum in SI units
    double momentumCutOff = 0.;

    /// Cut-off value for the step size
    double stepSizeCutOff = 0.;

	/// Data container for the energy loss
    EnergyLossData elData;
};

}  // namespace detail
}  // namespace Acts
