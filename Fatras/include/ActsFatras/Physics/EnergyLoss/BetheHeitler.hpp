// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/Kernel/detail/RandomNumberDistributions.hpp"

namespace ActsFatras {

const double log_2 = std::log(2.);

/// The struct for the EnergyLoss physics list
///
/// Bethe-Heitler for electron brem description as described here:
/// "A Gaussian-mixture approximation of the Bethe–Heitler model of electron
/// energy loss by bremsstrahlung" R. Frühwirth
///
struct BetheHeitler {
  /// The flag to include BetheHeitler process or not
  bool betheHeitler = true;

  /// A scaling factor to
  double scaleFactor = 1.;

  /// @brief Call operator for the Bethe-Heitler energy loss
  ///
  /// @tparam generator_t is a random number generator type
  /// @tparam detector_t is the detector information type
  /// @tparam particle_t is the particle information type
  ///
  /// @param[in] generator is the random number generator
  /// @param[in] detector the detector information
  /// @param[in] particle the particle which is being scattered
  ///
  /// @return eventually produced photons
  template <typename generator_t, typename detector_t, typename particle_t>
  std::vector<particle_t> operator()(generator_t &generator,
                                     const detector_t &detector,
                                     particle_t &particle) const {
    // Do nothing if the flag is set to false
    if (not betheHeitler) {
      return {};
    }

    double tInX0 = detector.thickness() / detector.material().X0();

    // Take a random gamma-distributed value - depending on t/X0
    GammaDist gDist = GammaDist(tInX0 / log_2, 1.);

    double u = gDist(generator);
    double z = std::exp(-1. * u);
    double sampledEnergyLoss = std::abs(scaleFactor * particle.E() * (z - 1.));

    // apply the energy loss
    particle.energyLoss(sampledEnergyLoss);

    // @TODO return photons, needs particle_creator_t
    return {};
  }
};

}  // namespace ActsFatras
