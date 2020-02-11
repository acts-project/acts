// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/Interactions.hpp"
#include "ActsFatras/Utilities/LandauDistribution.hpp"

namespace ActsFatras {

/// @brief The struct for the EnergyLoss physics list
///
/// This generates the energy loss according to the Bethe-Bloch
/// description, applying a landau generated enery loss
///
/// It follows the interface of EnergyLoss samplers in Fatras
/// that could return radiated photons for further processing,
/// however, for the Bethe-Bloch application the return vector
/// is always 0.
struct BetheBloch {
  /// The flag to include BetheBloch process or not
  bool betheBloch = true;

  /// Scaling for most probable value
  double scaleFactorMPV = 1.;

  /// Scaling for Sigma
  double scaleFactorSigma = 1.;

  /// @brief Call operator for the Bethe Bloch energy loss
  ///
  /// @tparam generator_t is a random number generator type
  /// @tparam detector_t is the detector information type
  /// @tparam particle_t is the particle information type
  ///
  /// @param[in] generator is the random number generator
  /// @param[in] detector the detector information
  /// @param[in] particle the particle which is being scattered
  ///
  /// @return empty vector for BetheBloch - no secondaries created
  template <typename generator_t, typename detector_t, typename particle_t>
  std::vector<particle_t> operator()(generator_t &generator,
                                     const detector_t &detector,
                                     particle_t &particle) const {
    // Do nothing if the flag is set to false
    if (not betheBloch) {
      return {};
    }

    // Create a random landau distribution between in the intervall [0,1]
    LandauDistribution landauDist(0., 1.);
    double landau = landauDist(generator);
    double qop = particle.q() / particle.p();

    // @TODO Double investigate if we could do one call
    double energyLoss = Acts::computeEnergyLossLandau(
        detector, particle.pdg(), particle.m(), qop, particle.q());
    double energyLossSigma = Acts::computeEnergyLossLandauSigma(
        detector, particle.pdg(), particle.m(), qop, particle.q());

    // Simulate the energy loss
    double sampledEnergyLoss = scaleFactorMPV * std::fabs(energyLoss) +
                               scaleFactorSigma * energyLossSigma * landau;

    // Apply the energy loss
    particle.energyLoss(sampledEnergyLoss);

    // return empty children
    return {};
  }
};

}  // namespace ActsFatras
