// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

namespace Acts::detail {

struct MaterialAccumulator {
  double maxXOverX0Step = 0;

  ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();
  double initialMomentum = 0;

  MaterialSlab accumulatedMaterial;

  double varAngle = 0;
  double varPosition = 0;
  double covAnglePosition = 0;
  double molarElectronDensity = 0;

  bool isVacuum() const { return accumulatedMaterial.isVacuum(); }

  void reset() { *this = MaterialAccumulator(); }

  void initialize(double maxXOverX0Step_,
                  const ParticleHypothesis& particleHypothesis_,
                  double initialMomentum_) {
    reset();
    maxXOverX0Step = maxXOverX0Step_;
    particleHypothesis = particleHypothesis_;
    initialMomentum = initialMomentum_;
  }

  void accumulate(const MaterialSlab& slab, double qOverPin, double qOverPout) {
    double mass = particleHypothesis.mass();
    double absQ = particleHypothesis.absoluteCharge();
    PdgParticle absPdg = particleHypothesis.absolutePdg();

    double momentumIn = particleHypothesis.extractMomentum(qOverPin);
    double momentumOut = particleHypothesis.extractMomentum(qOverPout);

    std::size_t substepCount =
        slab.isVacuum() ? 1
                        : static_cast<std::size_t>(
                              std::ceil(slab.thicknessInX0() / maxXOverX0Step));
    double substep = slab.thickness() / substepCount;
    MaterialSlab subslab(slab.material(), substep);

    for (std::size_t i = 0; i < substepCount; ++i) {
      double momentumMean =
          momentumIn + (momentumOut - momentumIn) * (i + 0.5) / substepCount;
      double qOverPmean = particleHypothesis.qOverP(momentumMean, absQ);

      double theta0in = computeMultipleScatteringTheta0(
          accumulatedMaterial, absPdg, mass, qOverPmean, absQ);

      molarElectronDensity =
          (molarElectronDensity * accumulatedMaterial.thickness() +
           subslab.material().molarElectronDensity() * subslab.thickness()) /
          (accumulatedMaterial.thickness() + subslab.thickness());
      accumulatedMaterial =
          MaterialSlab::combineLayers(accumulatedMaterial, subslab);

      double theta0out = computeMultipleScatteringTheta0(
          accumulatedMaterial, absPdg, mass, qOverPmean, absQ);

      double deltaVarTheta = square(theta0out) - square(theta0in);
      double deltaVarPos = varAngle * square(substep) +
                           2 * covAnglePosition * substep +
                           deltaVarTheta * (square(substep) / 3);
      double deltaCovAnglePosition =
          varAngle * substep + deltaVarTheta * substep / 2;
      varAngle += deltaVarTheta;
      varPosition += deltaVarPos;
      covAnglePosition += deltaCovAnglePosition;
    }
  }

  std::optional<FreeMatrix> computeAdditionalFreeCovariance(
      const Vector3& direction) {
    if (isVacuum()) {
      return std::nullopt;
    }

    FreeMatrix additionalFreeCovariance = FreeMatrix::Zero();

    // handle multiple scattering
    {
      // for derivation see
      // https://github.com/andiwand/cern-scripts/blob/5f0ebf1bef35db65322f28c2e840c1db1aaaf9a7/notebooks/2023-12-07_qp-dense-nav.ipynb
      //
      SquareMatrix3 directionProjection =
          (ActsSquareMatrix<3>::Identity() - direction * direction.transpose());

      additionalFreeCovariance.block<3, 3>(eFreeDir0, eFreeDir0) =
          varAngle * directionProjection;
      additionalFreeCovariance.block<3, 3>(eFreePos0, eFreePos0) =
          varPosition * directionProjection;
      additionalFreeCovariance.block<3, 3>(eFreePos0, eFreeDir0) =
          covAnglePosition * directionProjection;
      additionalFreeCovariance.block<3, 3>(eFreeDir0, eFreePos0) =
          additionalFreeCovariance.block<3, 3>(eFreePos0, eFreeDir0);
    }

    // handle energy loss covariance
    {
      double mass = particleHypothesis.mass();
      double absQ = particleHypothesis.absoluteCharge();
      double qOverP = particleHypothesis.qOverP(
          initialMomentum, particleHypothesis.absoluteCharge());

      Material other = accumulatedMaterial.material();
      Material tmp = Material::fromMolarDensity(
          other.X0(), other.L0(), other.Ar(),
          molarElectronDensity / other.molarDensity(), other.molarDensity());
      MaterialSlab tmpslab(tmp, accumulatedMaterial.thickness());

      float qOverPSigma =
          computeEnergyLossLandauSigmaQOverP(tmpslab, mass, qOverP, absQ);

      additionalFreeCovariance(eFreeQOverP, eFreeQOverP) =
          qOverPSigma * qOverPSigma;

      // in principle the energy loss uncertainty also affects the time
      // uncertainty continuously. these terms are not included here.
    }

    return additionalFreeCovariance;
  }
};

}  // namespace Acts::detail
