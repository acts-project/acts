// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"

#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/Interactions.hpp"

#include <cassert>
#include <limits>

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

namespace Acts {

MaterialUpdateMode detail::determineMaterialUpdateMode(
    const Surface& surface, const Surface* startSurface,
    const Surface* targetSurface, MaterialUpdateMode requestedMode) {
  // Avoid double counting of the material effects at start and target surfaces
  // by restricting the update mode accordingly
  MaterialUpdateMode updateMode = requestedMode;
  if (&surface == startSurface) {
    updateMode &= MaterialUpdateMode::PostUpdate;
  } else if (&surface == targetSurface) {
    updateMode &= MaterialUpdateMode::PreUpdate;
  }

  return updateMode;
}

Result<MaterialSlab> detail::evaluateMaterialSlab(
    const GeometryContext& geoContext, const Surface& surface,
    Direction propagationDirection, const Vector3& position,
    const Vector3& direction, MaterialUpdateMode updateMode) {
  const ISurfaceMaterial* material = surface.surfaceMaterial();
  if (material == nullptr) {
    return MaterialSlab::Nothing();
  }

  const double pathCorrection =
      surface.pathCorrection(geoContext, position, direction);
  auto lposition = surface.globalToLocal(geoContext, position, direction);
  if (!lposition.ok()) {
    return lposition.error();
  }

  MaterialSlab slab =
      surface.materialSlab(lposition.value(), propagationDirection, updateMode);
  slab.scaleThickness(pathCorrection);
  return slab;
}

detail::PointwiseMaterialEffects detail::computeMaterialEffects(
    const MaterialSlab& slab, const ParticleHypothesis& particleHypothesis,
    float qOverP, bool multipleScattering, bool energyLoss, bool covTransport) {
  PointwiseMaterialEffects result;

  const double mass = particleHypothesis.mass();
  const PdgParticle absPdg = particleHypothesis.absolutePdg();
  const double absQ = particleHypothesis.absoluteCharge();

  if (energyLoss) {
    result.eLoss = computeEnergyLossBethe(slab, mass, qOverP, absQ);
  }

  if (covTransport) {
    if (multipleScattering) {
      //! [scattering variance]
      const double theta0 =
          computeMultipleScatteringTheta0(slab, absPdg, mass, qOverP, absQ);
      result.scatteringVariance = theta0 * theta0;
      //! [scattering variance]
    }
    if (energyLoss) {
      //! [energy loss variance]
      const double sigmaQoverP =
          computeEnergyLossLandauSigmaQOverP(slab, mass, qOverP, absQ);
      result.varianceQoverP = sigmaQoverP * sigmaQoverP;
      //! [energy loss variance]
    }
  }

  return result;
}

Result<void> detail::removeMaterialCovarianceNoise(
    BoundMatrix& covariance, const Matrix<2, 2>& directionNoise,
    double varianceQoverP) {
  if (directionNoise.isZero(0.) && varianceQoverP == 0.) {
    return Result<void>::success();
  }

  // Factor a single scratch matrix in place. The input stays untouched until
  // validation succeeds; committing needs only the five changed entries.
  BoundMatrix scratch;
  const auto prepare = [&] {
    scratch = covariance;
    scratch.block<2, 2>(eBoundDir0, eBoundDir0) -= directionNoise;
    scratch(eBoundQOverP, eBoundQOverP) -= varianceQoverP;
  };
  const auto commit = [&] {
    covariance.block<2, 2>(eBoundDir0, eBoundDir0) -= directionNoise;
    covariance(eBoundQOverP, eBoundQOverP) -= varianceQoverP;
    return Result<void>::success();
  };
  prepare();
  if (scratch.allFinite()) {
    const Eigen::LLT<Eigen::Ref<BoundMatrix>> cholesky(scratch);
    if (cholesky.info() == Eigen::Success && scratch.allFinite()) {
      return commit();
    }
  }

  // A failed in-place factorization overwrites scratch. Recreate the candidate
  // and use the established unit-scaled check for difficult numerical cases.
  prepare();
  BoundVector scale;
  for (int i = 0; i < eBoundSize; ++i) {
    const double variance = std::max(covariance(i, i), std::abs(scratch(i, i)));
    scale[i] = variance > 0. ? 1. / std::sqrt(variance) : 1.;
  }
  for (int i = 0; i < eBoundSize; ++i) {
    scratch.row(i) *= scale[i];
    scratch.col(i) *= scale[i];
  }
  if (!scratch.allFinite()) {
    return MaterialInteractionError::InvalidCovarianceAfterNoiseRemoval;
  }
  const Eigen::LLT<BoundMatrix> cholesky(scratch);
  if (cholesky.info() == Eigen::Success) {
    return commit();
  }
  const Eigen::SelfAdjointEigenSolver<BoundMatrix> solver(
      scratch, Eigen::EigenvaluesOnly);
  if (solver.info() != Eigen::Success) {
    return MaterialInteractionError::InvalidCovarianceAfterNoiseRemoval;
  }
  const double tolerance =
      64. * std::numeric_limits<double>::epsilon() *
      std::max(1., solver.eigenvalues().cwiseAbs().maxCoeff());
  if (solver.eigenvalues().minCoeff() < -tolerance) {
    return MaterialInteractionError::InvalidCovarianceAfterNoiseRemoval;
  }
  return commit();
}

Vector2 detail::globalAngularScatteringVariances(const Vector3& direction,
                                                 double scatteringVariance) {
  assert(direction.allFinite() &&
         std::abs(direction.squaredNorm() - 1.) < 1e-6);
  assert(std::isfinite(scatteringVariance) && scatteringVariance >= 0.);
  if (scatteringVariance == 0.) {
    return Vector2::Zero();
  }
  const double transverseSquared = direction.head<2>().squaredNorm();
  return {transverseSquared == 0. ? std::numeric_limits<double>::infinity()
                                  : scatteringVariance / transverseSquared,
          scatteringVariance};
}

}  // namespace Acts
