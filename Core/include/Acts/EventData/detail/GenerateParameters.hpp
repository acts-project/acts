// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/AngleHelpers.hpp"

#include <cmath>
#include <numbers>
#include <random>
#include <utility>

namespace Acts::detail::Test {

template <typename scalar_t, std::size_t kSize, typename generator_t>
inline auto generateParameters(generator_t& rng)
    -> Eigen::Matrix<scalar_t, kSize, 1> {
  using Scalar = scalar_t;
  using ParametersVector = Eigen::Matrix<scalar_t, kSize, 1>;

  std::normal_distribution<Scalar> standardNormal(0, 1);

  ParametersVector params;
  for (auto i = 0u; i < kSize; ++i) {
    params[i] = standardNormal(rng);
  }

  return params;
}

template <typename scalar_t, std::size_t kSize, typename generator_t>
inline auto generateCovariance(generator_t& rng)
    -> Eigen::Matrix<scalar_t, kSize, kSize> {
  using Scalar = scalar_t;
  using ParametersVector = Eigen::Matrix<scalar_t, kSize, 1>;
  using CovarianceMatrix = Eigen::Matrix<scalar_t, kSize, kSize>;

  std::normal_distribution<Scalar> standardNormal(0, 1);
  std::uniform_real_distribution<Scalar> distCorr(-1, 1);

  // generate standard deviations
  ParametersVector stddev;
  for (auto i = 0u; i < kSize; ++i) {
    stddev[i] = std::abs(standardNormal(rng));
  }
  // generate correlation matrix
  CovarianceMatrix corr;
  for (auto i = 0u; i < kSize; ++i) {
    corr(i, i) = 1;
    // only need generate the sub-diagonal elements
    for (auto j = 0u; j < i; ++j) {
      corr(i, j) = corr(j, i) = distCorr(rng);
    }
  }
  // construct the covariance matrix
  CovarianceMatrix cov = stddev.asDiagonal() * corr * stddev.asDiagonal();

  return cov;
}

/// Generate a random parameters vector and covariance matrix.
///
/// @return std:::pair<ParametersVector, CovarianceMatrix>
template <typename scalar_t, std::size_t kSize, typename generator_t>
inline auto generateParametersCovariance(generator_t& rng)
    -> std::pair<Eigen::Matrix<scalar_t, kSize, 1>,
                 Eigen::Matrix<scalar_t, kSize, kSize>> {
  auto params = generateParameters<scalar_t, kSize, generator_t>(rng);
  auto cov = generateCovariance<scalar_t, kSize, generator_t>(rng);
  return {params, cov};
}

struct GenerateBoundDirectionOptions {
  /// Low, high (exclusive) for the transverse direction angle.
  long double phiMin = -std::numbers::pi;
  long double phiMax = std::numbers::pi;

  /// Low, high (inclusive) for  the longitudinal direction angle.
  ///
  /// This intentionally uses theta instead of eta so it can represent the
  /// full direction space with finite values.
  ///
  /// @note This is the standard generation, for detector performance
  /// classification, where a flat distribution in eta can be useful,
  /// this can be set by the etaUniform flag;
  ///
  long double thetaMin = AngleHelpers::thetaFromEta(6.0);
  long double thetaMax = AngleHelpers::thetaFromEta(-6.0);

  bool etaUniform = true;
};

template <typename generator_t>
inline std::pair<long double, long double> generateBoundDirection(
    generator_t& rng, const GenerateBoundDirectionOptions& options) {
  using UniformReal = std::uniform_real_distribution<long double>;
  assert(options.thetaMin >= 0.f);
  assert(options.thetaMax <= std::numbers::pi);
  assert(options.thetaMin <= options.thetaMax);

  // since we want to draw the direction uniform on the unit sphere, we must
  // draw from cos(theta) instead of theta. see e.g.
  // https://mathworld.wolfram.com/SpherePointPicking.html
  // Get cosThetaMin from thetaMax and vice versa, because cos is
  // monothonical decreasing between [0, pi]
  long double cosThetaMin = std::cos(options.thetaMax);
  // ensure upper bound is included. see e.g.
  // https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
  long double cosThetaMax = std::nextafter(
      std::cos(options.thetaMin), std::numeric_limits<long double>::max());

  // in case we force uniform eta generation
  long double etaMin = Acts::AngleHelpers::etaFromTheta(options.thetaMax);
  long double etaMax = Acts::AngleHelpers::etaFromTheta(options.thetaMin);

  UniformReal phiDist(options.phiMin, options.phiMax);
  UniformReal cosThetaDist(cosThetaMin, cosThetaMax);
  UniformReal etaDist(etaMin, etaMax);

  // draw parameters
  long double phi = phiDist(rng);

  long double theta = 0;
  if (!options.etaUniform) {
    const long double cosTheta = cosThetaDist(rng);
    theta = std::acos(cosTheta);
  } else {
    const long double eta = etaDist(rng);
    theta = AngleHelpers::thetaFromEta(eta);
  }

  return {phi, theta};
}

struct GenerateQoverPOptions {
  /// Low, high (exclusive) for absolute/transverse momentum.
  long double pMin = 1 * UnitConstants::GeV;
  long double pMax = 100 * UnitConstants::GeV;

  /// Indicate if the momentum referse to transverse momentum
  bool pTransverse = true;

  /// Indicate if the momentum should be uniformly distributed in log space.
  bool pLogUniform = false;

  /// Charge of the parameters.
  long double charge = 1;

  /// Randomize the charge and flip the PDG particle number sign accordingly.
  bool randomizeCharge = true;
};

template <typename generator_t>
inline long double generateQoverP(generator_t& rng,
                                  const GenerateQoverPOptions& options,
                                  long double theta) {
  using UniformIndex = std::uniform_int_distribution<std::uint8_t>;
  using UniformReal = std::uniform_real_distribution<long double>;

  auto drawP = [&options](generator_t& rng_,
                          long double theta_) -> long double {
    const long double pTransverseScaling =
        options.pTransverse ? 1. / std::sin(theta_) : 1.;

    if (options.pLogUniform) {
      UniformReal pLogDist(std::log(options.pMin), std::log(options.pMax));
      return std::exp(pLogDist(rng_)) * pTransverseScaling;
    }

    UniformReal pDist(options.pMin, options.pMax);
    return pDist(rng_) * pTransverseScaling;
  };

  // choose between particle/anti-particle if requested
  // the upper limit of the distribution is inclusive
  UniformIndex particleTypeChoice(0u, options.randomizeCharge ? 1u : 0u);
  // (anti-)particle choice is one random draw but defines two properties
  const long double qChoices[] = {
      options.charge,
      -options.charge,
  };

  // draw parameters
  const std::uint8_t type = particleTypeChoice(rng);
  const long double q = qChoices[type];

  const long double p = drawP(rng, theta);
  const long double qOverP = (q != 0) ? q / p : 1 / p;

  return qOverP;
}

struct GenerateBoundParametersOptions {
  struct {
    long double loc0Mean = 0 * UnitConstants::mm;
    long double loc0Std = 1 * UnitConstants::mm;

    long double loc1Mean = 0 * UnitConstants::mm;
    long double loc1Std = 1 * UnitConstants::mm;

    long double timeMean = 0 * UnitConstants::ns;
    long double timeStd = 1 * UnitConstants::ns;
  } position;

  GenerateBoundDirectionOptions direction;

  GenerateQoverPOptions qOverP;
};

inline BoundVector someBoundParametersA() {
  return {0.0, 0.0, 0.0, std::numbers::pi / 2, 1.0, 0.0};
}

template <typename generator_t>
inline BoundVector generateBoundParameters(
    generator_t& rng, const GenerateBoundParametersOptions& options) {
  std::normal_distribution<long double> standardNormal(0, 1);

  const long double loc0 = options.position.loc0Mean +
                           options.position.loc0Std * standardNormal(rng);
  const long double loc1 = options.position.loc1Mean +
                           options.position.loc1Std * standardNormal(rng);

  auto [phi, theta] = generateBoundDirection(rng, options.direction);
  auto qOverP = generateQoverP(rng, options.qOverP, theta);

  const long double time = options.position.timeMean +
                           options.position.timeStd * standardNormal(rng);

  return {loc0, loc1, phi, theta, qOverP, time};
}

template <typename generator_t>
inline std::pair<BoundVector, BoundMatrix> generateBoundParametersCovariance(
    generator_t& rng, const GenerateBoundParametersOptions& options) {
  auto params = generateBoundParameters(rng, options);
  auto cov = generateCovariance<long double, eBoundSize>(rng);
  return {params, cov};
}

struct GenerateFreeParametersOptions {
  struct {
    long double xMean = 0 * UnitConstants::mm;
    long double xStd = 1 * UnitConstants::mm;

    long double yMean = 0 * UnitConstants::mm;
    long double yStd = 1 * UnitConstants::mm;

    long double zMean = 0 * UnitConstants::mm;
    long double zStd = 1 * UnitConstants::mm;

    long double timeMean = 0 * UnitConstants::ns;
    long double timeStd = 1 * UnitConstants::ns;
  } position;

  GenerateBoundDirectionOptions direction;

  GenerateQoverPOptions qOverP;
};

inline FreeVector someFreeParametersA() {
  return {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
}

template <typename generator_t>
inline FreeVector generateFreeParameters(
    generator_t& rng, const GenerateFreeParametersOptions& options) {
  std::normal_distribution<long double> standardNormal(0, 1);

  const long double x =
      options.position.xMean + options.position.xStd * standardNormal(rng);
  const long double y =
      options.position.yMean + options.position.yStd * standardNormal(rng);
  const long double z =
      options.position.zMean + options.position.zStd * standardNormal(rng);
  const long double time = options.position.timeMean +
                           options.position.timeStd * standardNormal(rng);

  auto [phi, theta] = generateBoundDirection(rng, options.direction);

  Vector3 direction = makeDirectionFromPhiTheta(phi, theta);

  auto qOverP = generateQoverP(rng, options.qOverP, theta);

  FreeVector freeParams;
  freeParams << x, y, z, time, direction, qOverP;
  return freeParams;
}

template <typename generator_t>
inline std::pair<FreeVector, FreeMatrix> generateFreeParametersCovariance(
    generator_t& rng, const GenerateFreeParametersOptions& options) {
  auto params = generateFreeParameters(rng, options);
  auto cov = generateCovariance<long double, eFreeSize>(rng);
  return {params, cov};
}

}  // namespace Acts::detail::Test
