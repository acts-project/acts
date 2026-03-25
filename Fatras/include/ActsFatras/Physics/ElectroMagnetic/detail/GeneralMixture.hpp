// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Material/Interactions.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <numbers>
#include <random>

namespace ActsFatras::detail {

/// Generate scattering angles using a general mixture model.
///
/// Emulates core and tail scattering as described in
///
///     General mixture model Fruehwirth, M. Liendl.
///     Comp. Phys. Comm. 141 (2001) 230-246
///
struct GeneralMixture {
  /// Steering parameter
  bool logInclude = true;
  /// Scale the mixture level
  double genMixtureScalor = 1.;

  /// Generate a single 3D scattering angle.
  ///
  /// @param[in]     generator is the random number generator
  /// @param[in]     slab      defines the passed material
  /// @param[in,out] particle  is the particle being scattered
  /// @return a 3d scattering angle
  ///
  /// @tparam generator_t is a RandomNumberEngine
  template <typename generator_t>
  double operator()(generator_t &generator, const Acts::MaterialSlab &slab,
                    Particle &particle) const {
    double theta = 0.0;

    if (particle.absolutePdg() != Acts::PdgParticle::eElectron) {
      //----------------------------------------------------------------------------
      // see Mixture models of multiple scattering: computation and simulation.
      // -
      // R.Fruehwirth, M. Liendl. -
      // Computer Physics Communications 141 (2001) 230â246
      //----------------------------------------------------------------------------

      // Decide which mixture is best
      //   beta² = (p/E)² = p²/(p² + m²) = 1/(1 + (m/p)²)
      // 1/beta² = 1 + (m/p)²
      //    beta = 1/sqrt(1 + (m/p)²)
      const double mOverP = particle.mass() / particle.absoluteMomentum();
      const double beta2Inv = 1 + mOverP * mOverP;
      const double beta = 1 / std::sqrt(beta2Inv);
      const double tInX0 = slab.thicknessInX0();
      const double tob2 = tInX0 * beta2Inv;
      if (tob2 > 0.6 / std::pow(slab.material().Z(), 0.6)) {
        std::array<double, 4> scatteringParams{};
        // Gaussian mixture or pure Gaussian
        if (tob2 > 10) {
          scatteringParams = createGaussian(beta, particle.absoluteMomentum(),
                                            tInX0, genMixtureScalor);
        } else {
          scatteringParams =
              createGaussmix(beta, particle.absoluteMomentum(), tInX0,
                             slab.material().Z(), genMixtureScalor);
        }
        // Simulate
        theta = gaussmix(generator, scatteringParams);
      } else {
        // Semigaussian mixture - get parameters
        const std::array<double, 6> scatteringParamsSg =
            createSemigauss(beta, particle.absoluteMomentum(), tInX0,
                            slab.material().Z(), genMixtureScalor);
        // Simulate
        theta = semigauss(generator, scatteringParamsSg);
      }
    } else {
      // for electrons we fall back to the Highland (extension)
      // return projection factor times sigma times gauss random
      const auto theta0 = Acts::computeMultipleScatteringTheta0(
          slab, particle.absolutePdg(), particle.mass(), particle.qOverP(),
          particle.absoluteCharge());
      theta = std::normal_distribution<double>(0.0, theta0)(generator);
    }
    // scale from planar to 3d angle
    return std::numbers::sqrt2 * theta;
  }

  // helper methods for getting parameters and simulating

  std::array<double, 4> createGaussian(double beta, double p, double tInX0,
                                       double scale) const {
    std::array<double, 4> scatteringParams{};
    // Total standard deviation of mixture
    scatteringParams[0] = 15. / beta / p * std::sqrt(tInX0) * scale;
    // Variance of core
    scatteringParams[1] = 1.0;
    // Variance of tails
    scatteringParams[2] = 1.0;
    // Mixture weight of tail component
    scatteringParams[3] = 0.5;
    return scatteringParams;
  }

  std::array<double, 4> createGaussmix(double beta, double p, double tInX0,
                                       double Z, double scale) const {
    std::array<double, 4> scatteringParams{};
    // Total standard deviation of mixture
    scatteringParams[0] = 15. / beta / p * std::sqrt(tInX0) * scale;
    const double d1 = std::log(tInX0 / (beta * beta));
    const double d2 = std::log(std::pow(Z, 2.0 / 3.0) * tInX0 / (beta * beta));
    const double var1 = (-1.843e-3 * d1 + 3.347e-2) * d1 + 8.471e-1;
    double epsi = 0;
    if (d2 < 0.5) {
      epsi = (6.096e-4 * d2 + 6.348e-3) * d2 + 4.841e-2;
    } else {
      epsi = (-5.729e-3 * d2 + 1.106e-1) * d2 - 1.908e-2;
    }
    // Variance of core
    scatteringParams[1] = var1;
    // Variance of tails
    scatteringParams[2] = (1 - (1 - epsi) * var1) / epsi;
    // Mixture weight of tail component
    scatteringParams[3] = epsi;
    return scatteringParams;
  }

  std::array<double, 6> createSemigauss(double beta, double p, double tInX0,
                                        double Z, double scale) const {
    std::array<double, 6> scatteringParams{};
    const double N = tInX0 * 1.587E7 * std::pow(Z, 1.0 / 3.0) / (beta * beta) /
                     (Z + 1) / std::log(287 / std::sqrt(Z));
    // Total standard deviation of mixture
    scatteringParams[4] = 15. / beta / p * std::sqrt(tInX0) * scale;
    const double rho = 41000 / std::pow(Z, 2.0 / 3.0);
    const double b = rho / std::sqrt(N * (std::log(rho) - 0.5));
    const double n = std::pow(Z, 0.1) * std::log(N);
    const double var1 = (5.783E-4 * n + 3.803E-2) * n + 1.827E-1;
    const double a =
        (((-4.590E-5 * n + 1.330E-3) * n - 1.355E-2) * n + 9.828E-2) * n +
        2.822E-1;
    const double epsi = (1 - var1) / (a * a * (std::log(b / a) - 0.5) - var1);
    // Mixture weight of tail component
    scatteringParams[3] = (epsi > 0) ? epsi : 0.0;
    // Parameter 1 of tails
    scatteringParams[0] = a;
    // Parameter 2 of tails
    scatteringParams[1] = b;
    // Variance of core
    scatteringParams[2] = var1;
    // Average number of scattering processes
    scatteringParams[5] = N;
    return scatteringParams;
  }

  /// @brief Retrieve the gaussian mixture
  ///
  /// @tparam generator_t Type of the generator
  ///
  /// @param udist The uniform distribution handed over by the call operator
  /// @param scattering_params the tuned parameters for the generation
  ///
  /// @return a double value that represents the gaussian mixture
  template <typename generator_t>
  double gaussmix(generator_t &generator,
                  const std::array<double, 4> &scatteringParams) const {
    std::uniform_real_distribution<double> udist(0.0, 1.0);
    const double sigmaTot = scatteringParams[0];
    const double var1 = scatteringParams[1];
    const double var2 = scatteringParams[2];
    const double epsi = scatteringParams[3];
    const bool ind = udist(generator) > epsi;
    const double u = udist(generator);
    if (ind) {
      return std::sqrt(var1) * std::sqrt(-2 * std::log(u)) * sigmaTot;
    } else {
      return std::sqrt(var2) * std::sqrt(-2 * std::log(u)) * sigmaTot;
    }
  }

  /// @brief Retrieve the semi-gaussian mixture
  ///
  /// @tparam generator_t Type of the generator
  ///
  /// @param udist The uniform distribution handed over by the call operator
  /// @param scattering_params the tuned parameters for the generation
  ///
  /// @return a double value that represents the gaussian mixture
  template <typename generator_t>
  double semigauss(generator_t &generator,
                   const std::array<double, 6> &scattering_params) const {
    std::uniform_real_distribution<double> udist(0.0, 1.0);
    const double a = scattering_params[0];
    const double b = scattering_params[1];
    const double var1 = scattering_params[2];
    const double epsi = scattering_params[3];
    const double sigmaTot = scattering_params[4];
    const bool ind = udist(generator) > epsi;
    const double u = udist(generator);
    if (ind) {
      return std::sqrt(var1) * std::sqrt(-2 * std::log(u)) * sigmaTot;
    } else {
      return a * b * std::sqrt((1 - u) / (u * b * b + a * a)) * sigmaTot;
    }
  }
};

}  // namespace ActsFatras::detail
