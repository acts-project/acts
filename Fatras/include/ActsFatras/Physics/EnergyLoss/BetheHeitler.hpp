// This file is part of the Acts project.
//
// Copyright (C) 2018-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialSlab.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <array>
#include <random>
#include <iostream>
namespace ActsFatras {

/// Simulate electron energy loss using the Bethe-Heitler description.
///
/// Bethe-Heitler for electron bremsstrahlung description as described here:
/// "A Gaussian-mixture approximation of the Bethe–Heitler model of electron
/// energy loss by bremsstrahlung" R. Frühwirth
struct BetheHeitler {
  /// A scaling factor to
  double scaleFactor = 1.;

  /// Simulate energy loss and update the particle parameters.
  ///
  /// @param[in]     generator is the random number generator
  /// @param[in]     slab      defines the passed material
  /// @param[in,out] particle  is the particle being updated
  /// @return Produced photon.
  ///
  /// @tparam generator_t is a RandomNumberEngine
  template <typename generator_t>
  std::array<Particle, 1> operator()(generator_t &generator,
                                     const Acts::MaterialSlab &slab,
                                     Particle &particle) const {    
    // Take a random gamma-distributed value - depending on t/X0
    std::gamma_distribution<double> gDist(slab.thicknessInX0() / std::log(2.0),
                                          1.0);

    const auto u = gDist(generator);
    const auto z = std::exp(-u);
    const auto sampledEnergyLoss =
        std::abs(scaleFactor * particle.energy() * (z - 1.));
    const Particle::Scalar momentum = particle.absoluteMomentum();
	
    // apply the energy loss
    particle.correctEnergy(-sampledEnergyLoss);

	// Build the produced photon
    Particle photon(particle.particleId().makeDescendant(0), Acts::PdgParticle::eGamma);
    photon.setProcess(ProcessType::eBremsstrahlung).setPosition4(particle.fourPosition())
		.setDirection(particle.unitDirection()).setAbsoluteMomentum(std::min(sampledEnergyLoss, momentum));
	
    return {photon};
  }
};

}  // namespace ActsFatras
