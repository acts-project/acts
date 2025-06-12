// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <array>
#include <cfloat>
#include <cmath>
#include <numbers>
#include <random>

namespace ActsFatras {

/// Simulate electron energy loss using the Bethe-Heitler description.
///
/// Bethe-Heitler for electron bremsstrahlung description as described here:
/// "A Gaussian-mixture approximation of the Bethe–Heitler model of electron
/// energy loss by bremsstrahlung" R. Frühwirth
struct BetheHeitler {
  /// A scaling factor to
  double scaleFactor = 1.;

  // Simplified angle evaluation
  bool uniformHertzDipoleAngle = false;

  /// Simulate the photon emission
  ///
  /// @param [in] particle The unmodified electron
  /// @param [in] gammaE Energy of the photon
  /// @param [in] rndPsi Random number for the azimuthal angle
  /// @param [in] rndTheta1 Random number for the polar angle
  /// @param [in] rndTheta2 Random number for the polar angle
  /// @param [in] rndTheta3 Random number for the polar angle
  Particle bremPhoton(const Particle &particle, double gammaE, double rndPsi,
                      double rndTheta1, double rndTheta2,
                      double rndTheta3) const;

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
    std::uniform_real_distribution<double> uDist(0., 1.);

    // Take a random gamma-distributed value - depending on t/X0
    const double alpha = slab.thicknessInX0() / std::numbers::ln2;
    double u = 0.;
    if (alpha > 1.) {
      std::gamma_distribution<double> gDist(alpha, 1.);
      u = gDist(generator);
    } else {
      std::gamma_distribution<double> gDistPlus1(alpha + 1., 1.);
      // Get a random number using alpha+1
      const auto uPlus1 = gDistPlus1(generator);
      double u2;
      do {
        u2 = uDist(generator);
      } while (u2 == 0.);
      // Check if this would cause underflow
      if (std::log(u2) / alpha + std::log(uPlus1) < std::log(DBL_MIN)) {
        // This would be underflow - get 0.0 right away
        u = 0.;
      } else {
        u = uPlus1 * std::pow(u2, 1 / alpha);
      }
    }

    const auto z = std::exp(-u);
    const auto sampledEnergyLoss =
        std::abs(scaleFactor * particle.energy() * (z - 1.));

    // Build the produced photon
    Particle photon =
        bremPhoton(particle, sampledEnergyLoss, uDist(generator),
                   uDist(generator), uDist(generator), uDist(generator));
    // Recoil input momentum
    particle.setDirection(particle.direction() * particle.absoluteMomentum() -
                          photon.energy() * photon.direction());

    // apply the energy loss
    particle.correctEnergy(-sampledEnergyLoss);

    return {photon};
  }
};

}  // namespace ActsFatras
