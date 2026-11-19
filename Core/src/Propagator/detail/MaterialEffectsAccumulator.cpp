// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/MaterialEffectsAccumulator.hpp"

#include "Acts/Material/Interactions.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <cmath>

namespace Acts::detail {

namespace {

/// Remove `material` with the thickness `-substep` from the end of
/// `accumulated`. `combine` is additive in thickness, x/X0, l/L0 and molar
/// amount, so a negative thickness removes exactly what a positive one adds.
/// If `material` does not match what was accumulated, e.g. for binned volume
/// material, the subtraction can leave no or negative matter. Then remove the
/// average material instead.
MaterialSlab removeMaterial(const MaterialSlab& accumulated,
                            const Material& material, double substep) {
  const double thickness =
      static_cast<double>(accumulated.thickness()) + substep;
  if (thickness <= 0) {
    return MaterialSlab::Nothing();
  }

  // The check is on X0 and L0, because `MaterialSlab` sets x/X0 to zero for a
  // non-positive X0.
  const MaterialSlab combined =
      MaterialSlab::combine(accumulated, material, substep);
  if (!combined.isVacuum() && combined.material().X0() > 0 &&
      combined.material().L0() > 0) {
    return combined;
  }

  return MaterialSlab(accumulated.material(), static_cast<float>(thickness));
}

}  // namespace

void MaterialEffectsAccumulator::initialize(
    double maxXOverX0Step, const ParticleHypothesis& particleHypothesis,
    double initialMomentum) {
  reset();
  m_maxXOverX0Step = maxXOverX0Step;
  m_particleHypothesis = particleHypothesis;
  m_initialMomentum = initialMomentum;
}

void MaterialEffectsAccumulator::accumulate(const Material& material,
                                            double pathLength, double qOverPin,
                                            double qOverPout) {
  // A negative path length steps back along the path, e.g. to correct an
  // overstep. It removes material from the end of the path, but never more
  // than the path holds.
  // Removing more is not an error: the float thickness rounds below the double
  // path length, and the vacuum path before the first material step is not
  // accumulated. Nothing is left to remove, so the clamp is the correct result.
  if (pathLength < 0 &&
      -pathLength >= static_cast<double>(m_accumulatedMaterial.thickness())) {
    m_accumulatedMaterial = MaterialSlab::Nothing();
    m_varAngle = 0;
    m_varPosition = 0;
    m_covAnglePosition = 0;
    return;
  }

  const MaterialSlab slab(material, std::abs(pathLength));

  const float mass = m_particleHypothesis.mass();
  const float absQ = m_particleHypothesis.absoluteCharge();
  const PdgParticle absPdg = m_particleHypothesis.absolutePdg();

  const double momentumIn = m_particleHypothesis.extractMomentum(qOverPin);
  const double momentumOut = m_particleHypothesis.extractMomentum(qOverPout);

  const std::size_t substepCount =
      material.isVacuum() ? 1
                          : static_cast<std::size_t>(std::ceil(
                                slab.thicknessInX0() / m_maxXOverX0Step));
  const double substep = pathLength / substepCount;

  for (std::size_t i = 0; i < substepCount; ++i) {
    const double momentumMean =
        momentumIn + (momentumOut - momentumIn) * (i + 0.5) / substepCount;
    const double qOverPmean = m_particleHypothesis.qOverP(momentumMean, absQ);

    const double theta0in = computeMultipleScatteringTheta0(
        m_accumulatedMaterial, absPdg, mass, qOverPmean, absQ);

    m_accumulatedMaterial =
        substep >= 0
            ? MaterialSlab::combine(m_accumulatedMaterial, material, substep)
            : removeMaterial(m_accumulatedMaterial, material, substep);

    const double theta0out = computeMultipleScatteringTheta0(
        m_accumulatedMaterial, absPdg, mass, qOverPmean, absQ);

    // Integrate the scattering moments with a signed substep. A negative
    // substep integrates backwards, so no extra sign is needed.
    const double deltaVarTheta = square(theta0out) - square(theta0in);
    const double deltaVarPos = m_varAngle * square(substep) +
                               2 * m_covAnglePosition * substep +
                               deltaVarTheta * (square(substep) / 3);
    const double deltaCovAnglePosition =
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
        (SquareMatrix<3>::Identity() - direction * direction.transpose());

    additionalFreeCovariance.template block<3, 3>(eFreeDir0, eFreeDir0) =
        m_varAngle * directionProjection;
    additionalFreeCovariance.template block<3, 3>(eFreePos0, eFreePos0) =
        m_varPosition * directionProjection;
    additionalFreeCovariance.template block<3, 3>(eFreePos0, eFreeDir0) =
        m_covAnglePosition * directionProjection;
    additionalFreeCovariance.template block<3, 3>(eFreeDir0, eFreePos0) =
        additionalFreeCovariance.template block<3, 3>(eFreePos0, eFreeDir0);
  }

  // handle energy loss covariance
  {
    const double mass = m_particleHypothesis.mass();
    const double absQ = m_particleHypothesis.absoluteCharge();
    const double qOverP = m_particleHypothesis.qOverP(
        m_initialMomentum, m_particleHypothesis.absoluteCharge());

    const double qOverPSigma = computeEnergyLossLandauSigmaQOverP(
        m_accumulatedMaterial, mass, qOverP, absQ);

    additionalFreeCovariance(eFreeQOverP, eFreeQOverP) =
        qOverPSigma * qOverPSigma;

    // in principle the energy loss uncertainty also affects the time
    // uncertainty continuously. these terms are not included here.
  }

  return additionalFreeCovariance;
}

}  // namespace Acts::detail
