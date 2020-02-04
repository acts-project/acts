// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <random>

#include "Acts/Material/Interactions.hpp"

namespace ActsFatras {

/// @The struct for the Highland-based scattering
///
/// This will scatter particles with a single gaussian distribution
/// according to the highland formula.
struct Highland {
  /// @brief Call operator to perform this scattering
  ///
  /// @tparam generator_t is a random number generator type
  /// @tparam detector_t is the detector information type
  /// @tparam particle_t is the particle information type
  ///
  /// @param[in] generator is the random number generator
  /// @param[in] detector the detector information
  /// @param[in] particle the particle which is being scattered
  ///
  /// @return a scattering angle in 3D
  template <typename generator_t, typename detector_t, typename particle_t>
  double operator()(generator_t &generator, const detector_t &detector,
                    particle_t &particle) const {
    // Gauss distribution, will be sampled sampled with generator
    std::normal_distribution<double> gaussDist(0., 1.);

    double qop = particle.q() / particle.p();
    double theta0 = Acts::computeMultipleScatteringTheta0(
        detector, particle.pdg(), particle.m(), qop, particle.q());
    // Return projection factor times sigma times grauss random
    return M_SQRT2 * theta0 * gaussDist(generator);
  }
};

}  // namespace ActsFatras
