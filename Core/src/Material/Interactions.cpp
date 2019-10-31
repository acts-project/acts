// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/Interactions.hpp"

#include <cassert>
#include <cmath>

using namespace Acts::UnitLiterals;

namespace {

// PDG particle type identifiers
enum ParticleType : int {
  Electron = 11,
  Positron = -11,
  Muon = 13,
  AntiMuon = -13,
};

// values from RPP2018 table 33.1
// electron mass
static constexpr float Me = 0.5109989461_MeV;
// Bethe formular prefactor
static constexpr float K = 0.307075_MeV * 1_cm * 1_cm;
// Energy scale for plasma energy.
static constexpr float PlasmaEnergyScale = 28.816_eV;

// defined in ATL-SOFT-PUB-2008-003
static constexpr float MeanExitationPotentialScale = 16_eV;
// 4 * scaling factor for Landau FWHM to equivalent Gaussian sigma
// values is 2 / sqrt(2 * log(2))
static constexpr float Landau2Gauss = 2.577567883f;

/// Additional derived relativistic quantities.
struct RelativisticQuantities {
  float q2OverBeta2 = 0.0f;
  float beta2 = 0.0f;
  float betaGamma = 0.0f;
  float gamma = 0.0f;

  RelativisticQuantities(float m, float qOverP, float q) {
    // beta²/q² = (p/E)²/q² = p²/(q²m² + q²p²) = 1/(q² + (m²(q/p)²)
    // q²/beta² = q² + m²(q/p)²
    q2OverBeta2 = q * q + (m * qOverP) * (m * qOverP);
    // 1/p = q/(qp) = (q/p)/q
    const auto mOverP = m * std::abs(qOverP / q);
    const auto pOverM = 1.0f / mOverP;
    // beta² = p²/E² = p²/(m² + p²) = 1/(1 + (m/p)²)
    beta2 = 1.0f / (1.0f + mOverP * mOverP);
    // beta*gamma = (p/sqrt(m² + p²))*(sqrt(m² + p²)/m) = p/m
    betaGamma = pOverM;
    // gamma = sqrt(m² + p²)/m = sqrt(1 + (p/m)²)
    gamma = std::sqrt(1.0f + pOverM * pOverM);
  }
};

/// Compute q/p derivative of beta².
inline float deriveBeta2(float qOverP, const RelativisticQuantities& rq) {
  return -2 / (qOverP * rq.gamma * rq.gamma);
}

/// Compute the 2 * m * (beta * gamma)² mass term.
inline float computeMassTerm(float m, const RelativisticQuantities& rq) {
  return 2 * m * rq.betaGamma * rq.betaGamma;
}
/// Compute mass term logarithmic derivative w/ respect to q/p.
inline float logDeriveMassTerm(float qOverP) {
  // only need to compute d((beta*gamma)²)/(beta*gamma)²; rest cancels.
  return -2 / qOverP;
}

/// Compute the maximum energy transfer in a single collision.
///
/// Uses RPP2018 eq. 33.4.
inline float computeWMax(float m, const RelativisticQuantities& rq) {
  const auto mfrac = Me / m;
  const auto nominator = 2 * Me * rq.betaGamma * rq.betaGamma;
  const auto denonimator = 1.0f + 2 * rq.gamma * mfrac + mfrac * mfrac;
  return nominator / denonimator;
}
/// Compute WMax logarithmic derivative w/ respect to q/p.
inline float logDeriveWMax(float m, float qOverP,
                           const RelativisticQuantities& rq) {
  // this is (q/p) * (beta/q).
  // both quantities have the same sign and the product must always be
  // positive. we can thus reuse the known (unsigned) quantity (q/beta)².
  const auto a = std::abs(qOverP / std::sqrt(rq.q2OverBeta2));
  // (m² + me²) / me = me (1 + (m/me)²)
  const auto b = Me * (1.0f + (m / Me) * (m / Me));
  return -2 * (a * b - 2 + rq.beta2) / (qOverP * (a * b + 2));
}

/// Compute the mean excitation potential.
///
/// See ATL-SOFT-PUB-2008-003 equation (4)
inline float computeMeanExcitationPotential(float nuclearCharge) {
  return MeanExitationPotentialScale * std::pow(nuclearCharge, 0.9f);
}

/// Compute epsilon energy pre-factor for RPP2018 eq. 33.11.
///
/// Defined as
///
///     (K/2) * (Z/A)*rho * x * (q²/beta²)
///
/// where (Z/A)*rho it the electron density in the material and x is the
/// traversed length (thickness) of the material.
inline float computeEpsilon(float electronDensity, float thickness,
                            const RelativisticQuantities& rq) {
  return 0.5f * K * electronDensity * thickness * rq.q2OverBeta2;
}
/// Compute epsilon logarithmic derivative w/ respect to q/p.
inline float logDeriveEpsilon(float qOverP, const RelativisticQuantities& rq) {
  // only need to compute d(q²/beta²)/(q²/beta²); everything else cancels.
  return 2 / (qOverP * rq.gamma * rq.gamma);
}

/// Compute the density correction factor delta/2.
///
/// Uses RPP2018 eq. 33.6 which is only valid for high energies.
///
/// @todo Should we use RPP2018 eq. 33.7 instead w/ tabulated constants?
inline float computeDeltaHalf(float meanExitationPotential,
                              float electronDensity,
                              const RelativisticQuantities& rq) {
  // only relevant for very high ernergies; use arbitrary cutoff
  if (rq.betaGamma < 10.0f) {
    return 0.0f;
  }
  // pre-factor according to RPP2019 table 33.1
  const auto plasmaEnergy = PlasmaEnergyScale * std::sqrt(electronDensity);
  return std::log(rq.betaGamma) +
         std::log(plasmaEnergy / meanExitationPotential) - 0.5f;
}
/// Compute derivative w/ respect to q/p for the density correction.
inline float deriveDeltaHalf(float qOverP, const RelativisticQuantities& rq) {
  // original equation is of the form
  //     log(beta*gamma) + log(eplasma/I) - 1/2
  // which the resulting derivative as
  //     d(beta*gamma)/(beta*gamma)
  return (rq.betaGamma < 10.0f) ? 0.0f : (-1.0f / qOverP);
}

}  // namespace

#define ASSERT_INPUTS(thickness, m, qOverP, q)            \
  assert(0 < thickness and "Thickness must be positive"); \
  assert(0 < m and "Mass must be positive");              \
  assert(0 < (qOverP * q) and "Inconsistent q/p and q signs");

std::pair<float, float> Acts::computeIonisationLossMean(
    const Material& material, float thickness, float m, float qOverP, float q) {
  ASSERT_INPUTS(thickness, m, qOverP, q)

  // return early in case of vacuum
  if (not material) {
    return {0.0f, 0.0f};
  }

  const auto I = computeMeanExcitationPotential(material.Z());
  const auto rq = RelativisticQuantities(m, qOverP, q);
  const auto eps = computeEpsilon(material.zOverAtimesRho(), thickness, rq);
  const auto dhalf = computeDeltaHalf(I, material.zOverAtimesRho(), rq);
  const auto u = computeMassTerm(Me, rq);
  const auto wmax = computeWMax(m, rq);
  // uses RPP2018 eq. 33.5 scaled from mass stopping power to linear stopping
  // power and multiplied with the material thickness to get a total energy loss
  // instead of an energy loss per length.
  // the required modification only change the prefactor which becomes
  // identical to the prefactor epsilon for the most probable value.
  const auto running =
      0.5f * std::log(u / I) + 0.5f * std::log(wmax / I) - rq.beta2 - dhalf;
  return {eps * running, eps * Landau2Gauss};
}

float Acts::deriveIonisationLossMeanQOverP(const Material& material,
                                           float thickness, float m,
                                           float qOverP, float q) {
  ASSERT_INPUTS(thickness, m, qOverP, q)

  // return early in case of vacuum
  if (not material) {
    return 0.0f;
  }

  const auto I = computeMeanExcitationPotential(material.Z());
  const auto rq = RelativisticQuantities(m, qOverP, q);
  const auto eps = computeEpsilon(material.zOverAtimesRho(), thickness, rq);
  const auto dhalf = computeDeltaHalf(I, material.zOverAtimesRho(), rq);
  const auto u = computeMassTerm(Me, rq);
  const auto wmax = computeWMax(m, rq);
  // original equation is of the form
  //
  //     eps * (log(u/I)/2 + log(wmax/I)/2 - beta² - delta/2)
  //     = eps * (log(u)/2 + log(wmax/I)/2 - log(I) - beta² - delta/2)
  //
  // with the resulting derivative as
  //
  //     d(eps) * (log(u/I)/2 + log(wmax/I)/2 - beta² - delta/2)
  //     + eps * (d(u)/(2*u) + d(wmax)/(2*wmax) - d(beta²) - d(delta/2))
  //
  // where we can use d(eps) = eps * (d(eps)/eps) for further simplification.
  const auto logDerEps = logDeriveEpsilon(qOverP, rq);
  const auto derDHalf = deriveDeltaHalf(qOverP, rq);
  const auto logDerU = logDeriveMassTerm(qOverP);
  const auto logDerWmax = logDeriveWMax(m, qOverP, rq);
  const auto derBeta2 = deriveBeta2(qOverP, rq);
  const auto rel = logDerEps * (0.5f * std::log(u / I) +
                                0.5f * std::log(wmax / I) - rq.beta2 - dhalf) +
                   0.5f * logDerU + 0.5f * logDerWmax - derBeta2 - derDHalf;
  return eps * rel;
}

std::pair<float, float> Acts::computeIonisationLossMode(
    const Material& material, float thickness, float m, float qOverP, float q) {
  ASSERT_INPUTS(thickness, m, qOverP, q)

  // return early in case of vacuum
  if (not material) {
    return {0.0f, 0.0f};
  }

  const auto I = computeMeanExcitationPotential(material.Z());
  const auto rq = RelativisticQuantities(m, qOverP, q);
  const auto eps = computeEpsilon(material.zOverAtimesRho(), thickness, rq);
  const auto dhalf = computeDeltaHalf(I, material.zOverAtimesRho(), rq);
  const auto t = computeMassTerm(m, rq);
  // uses RPP2018 eq. 33.11
  const auto running =
      std::log(t / I) + std::log(eps / I) + 0.2f - rq.beta2 - 2 * dhalf;
  return {eps * running, eps * Landau2Gauss};
}

float Acts::deriveIonisationLossModeQOverP(const Material& material,
                                           float thickness, float m,
                                           float qOverP, float q) {
  ASSERT_INPUTS(thickness, m, qOverP, q)

  // return early in case of vacuum
  if (not material) {
    return 0.0f;
  }

  const auto I = computeMeanExcitationPotential(material.Z());
  const auto rq = RelativisticQuantities(m, qOverP, q);
  const auto eps = computeEpsilon(material.zOverAtimesRho(), thickness, rq);
  const auto dhalf = computeDeltaHalf(I, material.zOverAtimesRho(), rq);
  const auto t = computeMassTerm(m, rq);
  // original equation is of the form
  //
  //     eps * (log(t/I) - log(eps/I) - 0.2 - beta² - delta)
  //     = eps * (log(t) - log(eps) - 2*log(I) - 0.2 - beta² - delta)
  //
  // with the resulting derivative as
  //
  //     d(eps) * (log(t/I) - log(eps/I) - 0.2 - beta² - delta)
  //     + eps * (d(t)/t + d(eps)/eps - d(beta²) - 2*d(delta/2))
  //
  // where we can use d(eps) = eps * (d(eps)/eps) for further simplification.
  const auto logDerEps = logDeriveEpsilon(qOverP, rq);
  const auto derDHalf = deriveDeltaHalf(qOverP, rq);
  const auto logDerT = logDeriveMassTerm(qOverP);
  const auto derBeta2 = deriveBeta2(qOverP, rq);
  const auto rel = logDerEps * (std::log(t / I) + std::log(eps / I) - 0.2f -
                                rq.beta2 - 2 * dhalf) +
                   logDerT + logDerEps - derBeta2 - 2 * derDHalf;
  return eps * rel;
}

namespace {
/// Compute mean energy loss from bremsstrahlung per radiation length.
inline float computeBremsstrahlungLossMean(float m, float e) {
  return e * (Me / m) * (Me / m);
}
/// Derivative of the bremsstrahlung loss per rad length with respect to energy.
inline float deriveBremsstrahlungLossMeanE(float m) {
  return (Me / m) * (Me / m);
}

/// expansion coefficients for the muon radiation loss as a function of energy
///
/// taken from ATL-SOFT-PUB-2008-003 eq. 7,8 where the expansion is expressed
/// with terms e^n/X0 with fixed units [e] = MeV and [X0] = mm. the evaluated
/// expansion has units MeV/mm. In this implementation, the X0 dependence is
/// factored out and the coefficients must be scaled to the native units such
/// that the evaluated expansion with terms e^n has dimension energy in
/// native units.
constexpr float MuonHighLowThreshold = 1_TeV;
// [low0 / X0] = MeV / mm -> [low0] = MeV
constexpr double MuonLow0 = 0.5345_MeV;
// [low1 * E / X0] = MeV / mm -> [low1] = 1
constexpr double MuonLow1 = -6.803e-5;
// [low2 * E^2 / X0] = MeV / mm -> [low2] = 1/MeV
constexpr double MuonLow2 = -2.278e-11 / 1_MeV;
// [low3 * E^3 / X0] = MeV / mm -> [low3] = 1/MeV^2
constexpr double MuonLow3 = 9.899e-18 / (1_MeV * 1_MeV);
// same as low0
constexpr double MuonHigh0 = 2.986_MeV;
// same as low1
constexpr double MuonHigh1 = 9.253e-5;

/// Compute additional radiation energy loss for muons per radiation length.
inline float computeMuonDirectPairPhotoNuclearLossMean(double e) {
  if (e < MuonHighLowThreshold) {
    return MuonLow0 + MuonLow1 * e + MuonLow2 * e * e + MuonLow3 * e * e * e;
  } else {
    return MuonHigh0 + MuonHigh1 * e;
  }
}
/// Derivative of the additional rad loss per rad length with respect to energy.
inline float deriveMuonDirectPairPhotoNuclearLossMeanE(double e) {
  if (e < MuonHighLowThreshold) {
    return MuonLow1 + 2 * MuonLow2 * e + 3 * MuonLow3 * e * e;
  } else {
    return MuonHigh1;
  }
}
}  // namespace

float Acts::computeRadiationLoss(const Material& material, float thickness,
                                 int pdg, float m, float qOverP, float q) {
  ASSERT_INPUTS(thickness, m, qOverP, q)

  // return early in case of vacuum
  if (not material) {
    return 0.0f;
  }

  // relative radiation length
  const auto x = thickness / material.X0();
  // particle momentum and energy
  const auto p = q / qOverP;
  const auto e = std::sqrt(m * m + p * p);

  auto dEdx = computeBremsstrahlungLossMean(m, e);
  if (((pdg == ParticleType::Muon) or (pdg == ParticleType::AntiMuon)) and
      (8_GeV < e)) {
    dEdx += computeMuonDirectPairPhotoNuclearLossMean(e);
  }
  // scale from energy loss per unit radiation length to total energy
  return dEdx * x;
}

float Acts::deriveRadiationLossQOverP(const Material& material, float thickness,
                                      int pdg, float m, float qOverP, float q) {
  ASSERT_INPUTS(thickness, m, qOverP, q)

  // return early in case of vacuum
  if (not material) {
    return 0.0f;
  }

  // relative radiation length
  const auto x = thickness / material.X0();
  // particle momentum and energy
  const auto p = q / qOverP;
  const auto e = std::sqrt(m * m + p * p);

  // compute derivative w/ respect to energy.
  auto derE = deriveBremsstrahlungLossMeanE(m);
  if (((pdg == ParticleType::Muon) or (pdg == ParticleType::AntiMuon)) and
      (8_GeV < e)) {
    derE += deriveMuonDirectPairPhotoNuclearLossMeanE(e);
  }
  // compute derivative w/ respect to q/p by using the chain rule
  //     df(e)/d(q/p) = df(e)/de de/d(q/p)
  // with
  //     e = sqrt(m² + p²) = sqrt(m² + q²/(q/p)²)
  // and the resulting derivative
  //     de/d(q/p) = -q² / ((q/p)³ * e)
  const auto derQOverP = -(q * q) / (qOverP * qOverP * qOverP * e);
  return derE * derQOverP * x;
}

float Acts::computeEnergyLossMean(const Material& material, float thickness,
                                  int pdg, float m, float qOverP, float q) {
  return computeIonisationLossMean(material, thickness, m, qOverP, q).first +
         computeRadiationLoss(material, thickness, pdg, m, qOverP, q);
}

float Acts::deriveEnergyLossMeanQOverP(const Material& material,
                                       float thickness, int pdg, float m,
                                       float qOverP, float q) {
  return deriveIonisationLossMeanQOverP(material, thickness, m, qOverP, q) +
         deriveRadiationLossQOverP(material, thickness, pdg, m, qOverP, q);
}

float Acts::computeEnergyLossMode(const Material& material, float thickness,
                                  int pdg, float m, float qOverP, float q) {
  // see ATL-SOFT-PUB-2008-003 section 3 for the relative fractions
  // TODO this is inconsistent with the text of the note
  return 0.9f * computeIonisationLossMode(material, thickness, m, qOverP, q)
                    .first +
         0.15f * computeRadiationLoss(material, thickness, pdg, m, qOverP, q);
}

float Acts::deriveEnergyLossModeQOverP(const Material& material,
                                       float thickness, int pdg, float m,
                                       float qOverP, float q) {
  // see ATL-SOFT-PUB-2008-003 section 3 for the relative fractions
  // TODO this is inconsistent with the text of the note
  return 0.9f *
             deriveIonisationLossModeQOverP(material, thickness, m, qOverP, q) +
         0.15f *
             deriveRadiationLossQOverP(material, thickness, pdg, m, qOverP, q);
}
