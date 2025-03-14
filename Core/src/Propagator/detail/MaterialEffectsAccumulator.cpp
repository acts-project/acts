// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/MaterialEffectsAccumulator.hpp"

#include "Acts/Material/Interactions.hpp"

namespace Acts::detail {

void MaterialEffectsAccumulator::initialize(
    double maxXOverX0Step, const ParticleHypothesis& particleHypothesis,
    double initialMomentum) {
  reset();
  m_maxXOverX0Step = maxXOverX0Step;
  m_particleHypothesis = particleHypothesis;
  m_initialMomentum = initialMomentum;
}

void MaterialEffectsAccumulator::accumulate(const MaterialSlab& slab,
                                            double qOverPin, double qOverPout) {
  double mass = m_particleHypothesis.mass();
  double absQ = m_particleHypothesis.absoluteCharge();
  PdgParticle absPdg = m_particleHypothesis.absolutePdg();

  double momentumIn = m_particleHypothesis.extractMomentum(qOverPin);
  double momentumOut = m_particleHypothesis.extractMomentum(qOverPout);

  std::size_t substepCount =
      slab.isVacuum() ? 1
                      : static_cast<std::size_t>(
                            std::ceil(slab.thicknessInX0() / m_maxXOverX0Step));
  double substep = slab.thickness() / substepCount;
  MaterialSlab subslab(slab.material(), substep);

  for (std::size_t i = 0; i < substepCount; ++i) {
    double momentumMean =
        momentumIn + (momentumOut - momentumIn) * (i + 0.5) / substepCount;
    double qOverPmean = m_particleHypothesis.qOverP(momentumMean, absQ);

    double theta0in = Acts::computeMultipleScatteringTheta0(
        m_accumulatedMaterial, absPdg, mass, qOverPmean, absQ);

    m_molarElectronDensity =
        (m_molarElectronDensity * m_accumulatedMaterial.thickness() +
         subslab.material().molarElectronDensity() * subslab.thickness()) /
        (m_accumulatedMaterial.thickness() + subslab.thickness());
    m_accumulatedMaterial =
        MaterialSlab::combineLayers(m_accumulatedMaterial, subslab);

    double theta0out = Acts::computeMultipleScatteringTheta0(
        m_accumulatedMaterial, absPdg, mass, qOverPmean, absQ);

    double deltaVarTheta = square(theta0out) - square(theta0in);
    double deltaVarPos = m_varAngle * square(substep) +
                         2 * m_covAnglePosition * substep +
                         deltaVarTheta * (square(substep) / 3);
    double deltaCovAnglePosition =
        m_varAngle * substep + deltaVarTheta * substep / 2;
    m_varAngle += deltaVarTheta;
    m_varPosition += deltaVarPos;
    m_covAnglePosition += deltaCovAnglePosition;
  }
}

std::optional<FreeMatrix>
MaterialEffectsAccumulator::computeAdditionalFreeCovariance(
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
        m_varAngle * directionProjection;
    additionalFreeCovariance.block<3, 3>(eFreePos0, eFreePos0) =
        m_varPosition * directionProjection;
    additionalFreeCovariance.block<3, 3>(eFreePos0, eFreeDir0) =
        m_covAnglePosition * directionProjection;
    additionalFreeCovariance.block<3, 3>(eFreeDir0, eFreePos0) =
        additionalFreeCovariance.block<3, 3>(eFreePos0, eFreeDir0);
  }

  // handle energy loss covariance
  {
    double mass = m_particleHypothesis.mass();
    double absQ = m_particleHypothesis.absoluteCharge();
    double qOverP = m_particleHypothesis.qOverP(
        m_initialMomentum, m_particleHypothesis.absoluteCharge());

    Material other = m_accumulatedMaterial.material();
    Material tmp = Material::fromMolarDensity(
        other.X0(), other.L0(), other.Ar(),
        m_molarElectronDensity / other.molarDensity(), other.molarDensity());
    MaterialSlab tmpslab(tmp, m_accumulatedMaterial.thickness());

    float qOverPSigma =
        computeEnergyLossLandauSigmaQOverP(tmpslab, mass, qOverP, absQ);

    additionalFreeCovariance(eFreeQOverP, eFreeQOverP) =
        qOverPSigma * qOverPSigma;

    // in principle the energy loss uncertainty also affects the time
    // uncertainty continuously. these terms are not included here.
  }

  return additionalFreeCovariance;
}

}  // namespace Acts::detail
