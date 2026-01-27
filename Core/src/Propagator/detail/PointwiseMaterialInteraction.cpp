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
#include "Acts/Utilities/Helpers.hpp"

namespace Acts::detail {

bool PointwiseMaterialInteraction::evaluateMaterialSlab(
    MaterialUpdateMode requestedMode) {
  updateMode = requestedMode;
  if (surface == startSurface) {
    updateMode &= MaterialUpdateMode::PostUpdate;
  } else if (surface == targetSurface) {
    updateMode &= MaterialUpdateMode::PreUpdate;
  }

  // Retrieve the material properties
  const ISurfaceMaterial* material = surface->surfaceMaterial();
  slab = material->materialSlab(pos, propDir, updateMode);

  // Correct the material properties for non-zero incidence
  slab.scaleThickness(pathCorrection);

  return !slab.isVacuum();
}

void PointwiseMaterialInteraction::evaluatePointwiseMaterialInteraction(
    bool multipleScattering, bool energyLoss) {
  if (energyLoss) {
    Eloss = computeEnergyLossBethe(slab, mass, qOverP, absQ);
  }
  if (performCovarianceTransport) {
    evaluateCovarianceContributions(multipleScattering, energyLoss);
  }
}

void PointwiseMaterialInteraction::evaluateCovarianceContributions(
    bool multipleScattering, bool energyLoss) {
  if (multipleScattering) {
    const double theta0 =
        computeMultipleScatteringTheta0(slab, absPdg, mass, qOverP, absQ);
    // sigmaPhi = theta0 / sin(theta)
    const double sigmaPhi = theta0 * (dir.norm() / VectorHelpers::perp(dir));
    variancePhi = sigmaPhi * sigmaPhi;
    // sigmaTheta = theta0
    varianceTheta = theta0 * theta0;
  }
  if (energyLoss) {
    const double sigmaQoverP =
        computeEnergyLossLandauSigmaQOverP(slab, mass, qOverP, absQ);
    varianceQoverP = sigmaQoverP * sigmaQoverP;
  }
}

double PointwiseMaterialInteraction::updateVariance(
    double variance, double change, NoiseUpdateMode updateMode) {
  // Add/Subtract the change
  // Protect the variance against becoming negative
  return std::max(0.,
                  variance + std::copysign(change, toUnderlying(updateMode)));
}

}  // namespace Acts::detail
