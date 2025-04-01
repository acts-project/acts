// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/Interactions.hpp"

#include <numbers>
#include <random>

namespace ActsFatras::detail {

/// Generate scattering angles using a Gaussian mixture model.
struct GaussianMixture {
  /// Steering parameter
  bool optGaussianMixtureG4 = false;
  long double gausMixSigma1_a0 = 8.471e-1;
  long double gausMixSigma1_a1 = 3.347e-2;
  long double gausMixSigma1_a2 = -1.843e-3;
  long double gausMixEpsilon_a0 = 4.841e-2;
  long double gausMixEpsilon_a1 = 6.348e-3;
  long double gausMixEpsilon_a2 = 6.096e-4;
  long double gausMixEpsilon_b0 = -1.908e-2;
  long double gausMixEpsilon_b1 = 1.106e-1;
  long double gausMixEpsilon_b2 = -5.729e-3;

  /// Generate a single 3D scattering angle.
  ///
  /// @param[in]     generator is the random number generator
  /// @param[in]     slab      defines the passed material
  /// @param[in,out] particle  is the particle being scattered
  /// @return a 3d scattering angle
  ///
  /// @tparam generator_t is a RandomNumberEngine
  template <typename generator_t>
  long double operator()(generator_t &generator, const Acts::MaterialSlab &slab,
                         Particle &particle) const {
    /// Calculate the highland formula first
    long double sigma = Acts::computeMultipleScatteringTheta0(
        slab, particle.absolutePdg(), particle.mass(), particle.qOverP(),
        particle.absoluteCharge());
    long double sigma2 = sigma * sigma;

    // Gauss distribution, will be sampled with generator
    std::normal_distribution<long double> gaussDist(0., 1.);
    // Uniform distribution, will be sampled with generator
    std::uniform_real_distribution<long double> uniformDist(0., 1.);

    // Now correct for the tail fraction
    // d_0'
    // beta² = (p/E)² = p²/(p² + m²) = 1/(1 + (m/p)²)
    // 1/beta² = 1 + (m/p)²
    long double mOverP = particle.mass() / particle.absoluteMomentum();
    long double beta2inv = 1 + mOverP * mOverP;
    long double dprime = slab.thicknessInX0() * beta2inv;
    long double log_dprime = std::log(dprime);
    // d_0''
    long double log_dprimeprime =
        std::log(std::pow(slab.material().Z(), 2.0 / 3.0) * dprime);

    // get epsilon
    long double epsilon =
        log_dprimeprime < 0.5
            ? gausMixEpsilon_a0 + gausMixEpsilon_a1 * log_dprimeprime +
                  gausMixEpsilon_a2 * log_dprimeprime * log_dprimeprime
            : gausMixEpsilon_b0 + gausMixEpsilon_b1 * log_dprimeprime +
                  gausMixEpsilon_b2 * log_dprimeprime * log_dprimeprime;

    // the standard sigma
    long double sigma1square = gausMixSigma1_a0 +
                               gausMixSigma1_a1 * log_dprime +
                               gausMixSigma1_a2 * log_dprime * log_dprime;

    // G4 optimised / native long double Gaussian model
    if (optGaussianMixtureG4) {
      sigma2 = 225. * dprime /
               (particle.absoluteMomentum() * particle.absoluteMomentum());
    }
    // throw the random number core/tail
    if (uniformDist(generator) < epsilon) {
      sigma2 *= (1. - (1. - epsilon) * sigma1square) / epsilon;
    }
    // return back to the
    return std::numbers::sqrt2 * std::sqrt(sigma2) * gaussDist(generator);
  }
};

}  // namespace ActsFatras::detail
