// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>

#include "Acts/Material/Interactions.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Utilities/LandauDistribution.hpp"

namespace ActsFatras {

/// Simulate energy loss using the Bethe-Bloch/Landau description.
///
/// Energy loss is computed using the most probable value and appropriate
/// fluctuations from a Landau distribution. No secondaries are generated
/// for the removed energy.
struct BetheBloch {
  /// Scaling for most probable value
  double scaleFactorMPV = 1.;
  /// Scaling for Sigma
  double scaleFactorSigma = 1.;

  /// Simulate energy loss and update the particle parameters.
  ///
  /// @param[in]     generator is the random number generator
  /// @param[in]     slab      defines the passed material
  /// @param[in,out] particle  is the particle being updated
  /// @return Empty secondaries containers.
  ///
  /// @tparam generator_t is a RandomNumberEngine
  template <typename generator_t>
  std::array<Particle, 0> operator()(generator_t &generator,
                                     const Acts::MaterialProperties &slab,
                                     Particle &particle) const {
    // compute energy loss distribution parameters
    const auto pdg = particle.pdg();
    const auto m = particle.mass();
    const auto qOverP = particle.charge() / particle.absMomentum();
    const auto q = particle.charge();
    // most probable value
    const auto energyLoss =
        Acts::computeEnergyLossLandau(slab, pdg, m, qOverP, q);
    // Gaussian-equivalent sigma
    const auto energyLossSigma =
        Acts::computeEnergyLossLandauSigma(slab, pdg, m, qOverP, q);

    // Simulate the energy loss
    // TODO landau location and scale parameters are not identical to the most
    //      probable value and the Gaussian-equivalent sigma
    LandauDistribution lossDistribution(scaleFactorMPV * energyLoss,
                                        scaleFactorSigma * energyLossSigma);
    const auto loss = lossDistribution(generator);

    // Apply the energy loss
    particle.correctEnergy(-loss);

    // Generates no new particles
    return {};
  }
};

}  // namespace ActsFatras
