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

MaterialSlab detail::evaluateMaterialSlab(const GeometryContext& geoContext,
                                          const Surface& surface,
                                          Direction propagationDirection,
                                          const Vector3& position,
                                          const Vector3& direction,
                                          MaterialUpdateMode updateMode) {
  const ISurfaceMaterial* material = surface.surfaceMaterial();
  if (material == nullptr) {
    return MaterialSlab::Nothing();
  }

  const double pathCorrection =
      surface.pathCorrection(geoContext, position, direction);
  MaterialSlab slab =
      material->materialSlab(position, propagationDirection, updateMode);
  slab.scaleThickness(pathCorrection);

  return slab;
}

detail::PointwiseMaterialEffects detail::computeMaterialEffects(
    const MaterialSlab& slab, const ParticleHypothesis& particleHypothesis,
    const Vector3& direction, float qOverP, bool multipleScattering,
    bool energyLoss, bool covTransport) {
  PointwiseMaterialEffects result;

  const double mass = particleHypothesis.mass();
  const PdgParticle absPdg = particleHypothesis.absolutePdg();
  const double absQ = particleHypothesis.absoluteCharge();

  if (energyLoss) {
    result.eLoss = computeEnergyLossBethe(slab, mass, qOverP, absQ);
  }

  if (covTransport) {
    if (multipleScattering) {
      const double theta0 =
          computeMultipleScatteringTheta0(slab, absPdg, mass, qOverP, absQ);
      // sigmaPhi = theta0 / sin(theta)
      const double sigmaPhi =
          theta0 * (direction.norm() / VectorHelpers::perp(direction));
      result.variancePhi = sigmaPhi * sigmaPhi;
      // sigmaTheta = theta0
      result.varianceTheta = theta0 * theta0;
    }
    if (energyLoss) {
      const double sigmaQoverP =
          computeEnergyLossLandauSigmaQOverP(slab, mass, qOverP, absQ);
      result.varianceQoverP = sigmaQoverP * sigmaQoverP;
    }
  }

  return result;
}

}  // namespace Acts
