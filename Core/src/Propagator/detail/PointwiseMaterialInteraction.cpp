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

Result<Matrix<2, 2>> detail::scatteringNoise(
    const BoundParametrization& parametrization, const Vector3& direction,
    double scatteringVariance) {
  assert(direction.allFinite() &&
         std::abs(direction.squaredNorm() - 1.) < 1e-6);
  assert(std::isfinite(scatteringVariance) && scatteringVariance >= 0.);
  if (scatteringVariance == 0.) {
    return Result<Matrix<2, 2>>::success(Matrix<2, 2>::Zero());
  }
  if (parametrization.isGlobalSpherical()) {
    const Vector2 variances =
        globalAngularScatteringVariances(direction, scatteringVariance);
    if (!variances.allFinite()) {
      return Result<Matrix<2, 2>>::failure(
          BoundParametrizationError::SingularDirection);
    }
    return Result<Matrix<2, 2>>::success(Matrix<2, 2>(variances.asDiagonal()));
  }
  const auto jacobian =
      parametrization.directionCoordinatesJacobianFromUnitDirection(direction,
                                                                    0.);
  if (!jacobian.ok()) {
    return Result<Matrix<2, 2>>::failure(jacobian.error());
  }
  // The encoding derivative already includes the tangent projection I - d*d^T.
  return Result<Matrix<2, 2>>::success(scatteringVariance * (*jacobian) *
                                       jacobian->transpose());
}

}  // namespace Acts
