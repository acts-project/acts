// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "ActsTests/CommonHelpers/BenchmarkTools.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

using namespace Acts;
using namespace ActsTests;

// Previous scaled-first implementation, retained for same-run comparisons.
namespace {
[[gnu::noinline]] Result<void> scaledFirstRemoval(
    BoundMatrix& covariance, const Matrix<2, 2>& directionNoise,
    double varianceQoverP) {
  if (directionNoise.isZero(0.) && varianceQoverP == 0.) {
    return Result<void>::success();
  }

  BoundMatrix updated = covariance;
  updated.block<2, 2>(eBoundDir0, eBoundDir0) -= directionNoise;
  updated(eBoundQOverP, eBoundQOverP) -= varianceQoverP;

  // Use diagonal variance scales so the tolerance does not compare position,
  // angle, time and inverse-momentum units directly. Include the proposed
  // diagonal to detect even tiny negative variances when starting from zero.
  BoundVector scale;
  for (int i = 0; i < eBoundSize; ++i) {
    const double variance = std::max(covariance(i, i), std::abs(updated(i, i)));
    scale[i] = variance > 0. ? 1. / std::sqrt(variance) : 1.;
  }
  BoundMatrix normalized = updated;
  for (int i = 0; i < eBoundSize; ++i) {
    normalized.row(i) *= scale[i];
    normalized.col(i) *= scale[i];
  }
  if (!normalized.allFinite()) {
    return MaterialInteractionError::InvalidCovarianceAfterNoiseRemoval;
  }
  // Positive-definite results need only a Cholesky factorization. Fall back to
  // eigenvalues for semidefinite boundaries and the existing roundoff policy.
  const Eigen::LLT<BoundMatrix> cholesky(normalized);
  if (cholesky.info() == Eigen::Success) {
    covariance = updated;
    return Result<void>::success();
  }
  const Eigen::SelfAdjointEigenSolver<BoundMatrix> solver(
      normalized, Eigen::EigenvaluesOnly);
  if (solver.info() != Eigen::Success) {
    return MaterialInteractionError::InvalidCovarianceAfterNoiseRemoval;
  }
  const double tolerance =
      64. * std::numeric_limits<double>::epsilon() *
      std::max(1., solver.eigenvalues().cwiseAbs().maxCoeff());
  if (solver.eigenvalues().minCoeff() < -tolerance) {
    return MaterialInteractionError::InvalidCovarianceAfterNoiseRemoval;
  }
  covariance = updated;
  return Result<void>::success();
}

// Isolate the unscaled-first change from in-place factorization/copy reduction.
[[gnu::noinline]] Result<void> unscaledCopyingRemoval(BoundMatrix& covariance,
                                                      const Matrix<2, 2>& noise,
                                                      double varianceQoverP) {
  BoundMatrix updated = covariance;
  updated.block<2, 2>(eBoundDir0, eBoundDir0) -= noise;
  updated(eBoundQOverP, eBoundQOverP) -= varianceQoverP;
  if (updated.allFinite()) {
    const Eigen::LLT<BoundMatrix> factor(updated);
    if (factor.info() == Eigen::Success && factor.matrixLLT().allFinite()) {
      covariance = updated;
      return Result<void>::success();
    }
  }
  return scaledFirstRemoval(covariance, noise, varianceQoverP);
}
}  // namespace

int main() {
  // Fixed inputs, generated outside the timed region. Return the entire result
  // so the benchmark barriers keep all covariance writes observable.
  std::mt19937 generator(42);
  std::normal_distribution<double> normal;
  std::vector<BoundMatrix> inputs;
  for (int n = 0; n < 256; ++n) {
    BoundMatrix factor;
    for (int i = 0; i < factor.size(); ++i) {
      factor.data()[i] = normal(generator);
    }
    inputs.push_back(factor * factor.transpose() + BoundMatrix::Identity());
  }
  Matrix<2, 2> noise = Matrix<2, 2>::Zero();
  noise.diagonal() << 0.02, 0.01;
  const auto run = [&](const char* label, auto&& operation) {
    std::cout << label << ": "
              << microBenchmark(operation, inputs, 20000,
                                std::chrono::milliseconds(200))
              << '\n';
  };
  for (const auto mode :
       {NoiseUpdateMode::addNoise, NoiseUpdateMode::removeNoise}) {
    const double sign = toUnderlying(mode);
    run(mode == NoiseUpdateMode::addNoise ? "legacy add" : "legacy remove",
        [&](const BoundMatrix& input) {
          BoundMatrix covariance = input;
          covariance(eBoundDir0, eBoundDir0) = std::max(
              0., covariance(eBoundDir0, eBoundDir0) + sign * noise(0, 0));
          covariance(eBoundDir1, eBoundDir1) = std::max(
              0., covariance(eBoundDir1, eBoundDir1) + sign * noise(1, 1));
          covariance(eBoundQOverP, eBoundQOverP) = std::max(
              0., covariance(eBoundQOverP, eBoundQOverP) + sign * 0.001);
          return covariance;
        });
    run(mode == NoiseUpdateMode::addNoise ? "checked add" : "checked remove",
        [&](const BoundMatrix& input) {
          BoundMatrix covariance = input;
          const auto result =
              detail::updateMaterialCovariance(covariance, noise, 0.001, mode);
          if (!result.ok()) {
            std::abort();
          }
          return covariance;
        });
  }
  run("scaled-first checked remove", [&](const BoundMatrix& input) {
    BoundMatrix covariance = input;
    const auto result = scaledFirstRemoval(covariance, noise, 0.001);
    if (!result.ok()) {
      std::abort();
    }
    return covariance;
  });
  run("unscaled-first copying remove", [&](const BoundMatrix& input) {
    BoundMatrix covariance = input;
    const auto result = unscaledCopyingRemoval(covariance, noise, 0.001);
    if (!result.ok()) {
      std::abort();
    }
    return covariance;
  });
  const MaterialSlab slab(makeSilicon(), 1.);
  const auto particle = ParticleHypothesis::pion();
  const Vector3 direction = Vector3(3., 4., 5.).normalized();
  run("legacy physics+noise add", [&](const BoundMatrix& input) {
    const auto effects =
        detail::computeMaterialEffects(slab, particle, 0.1f, true, true, true);
    const double sigmaPhi = std::sqrt(effects.scatteringVariance) *
                            direction.norm() / direction.head<2>().norm();
    BoundMatrix covariance = input;
    covariance(eBoundDir0, eBoundDir0) += sigmaPhi * sigmaPhi;
    covariance(eBoundDir1, eBoundDir1) += effects.scatteringVariance;
    covariance(eBoundQOverP, eBoundQOverP) += effects.varianceQoverP;
    return std::pair{covariance, effects.eLoss};
  });
  run("checked physics+noise add", [&](const BoundMatrix& input) {
    const auto effects =
        detail::computeMaterialEffects(slab, particle, 0.1f, true, true, true);
    const auto mapped =
        detail::scatteringNoise(BoundParametrization::globalSpherical(),
                                direction, effects.scatteringVariance);
    if (!mapped.ok()) {
      std::abort();
    }
    BoundMatrix covariance = input;
    const auto result = detail::updateMaterialCovariance(
        covariance, *mapped, effects.varianceQoverP, NoiseUpdateMode::addNoise);
    if (!result.ok()) {
      std::abort();
    }
    return std::pair{covariance, effects.eLoss};
  });

  for (auto& input : inputs) {
    input.setZero();
    input.block<2, 2>(eBoundDir0, eBoundDir0) = noise;
    input(eBoundQOverP, eBoundQOverP) = 0.001;
  }
  run("checked semidefinite removal", [&](const BoundMatrix& input) {
    BoundMatrix covariance = input;
    const auto result = detail::updateMaterialCovariance(
        covariance, noise, 0.001, NoiseUpdateMode::removeNoise);
    if (!result.ok()) {
      std::abort();
    }
    return covariance;
  });
  run("scaled-first semidefinite removal", [&](const BoundMatrix& input) {
    BoundMatrix covariance = input;
    const auto result = scaledFirstRemoval(covariance, noise, 0.001);
    if (!result.ok()) {
      std::abort();
    }
    return covariance;
  });
}
