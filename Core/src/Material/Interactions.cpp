// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/Interactions.hpp"

#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/PdgParticle.hpp"

#include <cassert>
#include <cmath>

using namespace Acts::UnitLiterals;

namespace {

// values from RPP2018 table 33.1
// electron mass
constexpr float Me = 0.5109989461_MeV;
// Bethe formular prefactor. 1/mol unit is just a factor 1 here.
constexpr float K = 0.307075_MeV * 1_cm * 1_cm;
// Energy scale for plasma energy.
constexpr float PlasmaEnergyScale = 28.816_eV;

/// Compute q/p derivative of beta².
inline float deriveBeta2(float qOverP, const Acts::RelativisticQuantities& rq) {
  return -2 / (qOverP * rq.gamma * rq.gamma);
}

/// Compute the 2 * mass * (beta * gamma)² mass term.
inline float computeMassTerm(float mass,
                             const Acts::RelativisticQuantities& rq) {
  return 2 * mass * rq.betaGamma * rq.betaGamma;
}

/// Compute mass term logarithmic derivative w/ respect to q/p.
inline float logDeriveMassTerm(float qOverP) {
  // only need to compute d((beta*gamma)²)/(beta*gamma)²; rest cancels.
  return -2 / qOverP;
}

/// Compute the maximum energy transfer in a single collision.
///
/// Uses RPP2018 eq. 33.4.
inline float computeWMax(float mass, const Acts::RelativisticQuantities& rq) {
  const auto mfrac = Me / mass;
  const auto nominator = 2 * Me * rq.betaGamma * rq.betaGamma;
  const auto denonimator = 1.0f + 2 * rq.gamma * mfrac + mfrac * mfrac;
  return nominator / denonimator;
}

/// Compute WMax logarithmic derivative w/ respect to q/p.
inline float logDeriveWMax(float mass, float qOverP,
                           const Acts::RelativisticQuantities& rq) {
  // this is (q/p) * (beta/q).
  // both quantities have the same sign and the product must always be
  // positive. we can thus reuse the known (unsigned) quantity (q/beta)².
  const auto a = std::abs(qOverP / std::sqrt(rq.q2OverBeta2));
  // (m² + me²) / me = me (1 + (m/me)²)
  const auto b = Me * (1.0f + (mass / Me) * (mass / Me));
  return -2 * (a * b - 2 + rq.beta2) / (qOverP * (a * b + 2));
}

/// Compute epsilon energy pre-factor for RPP2018 eq. 33.11.
///
/// Defined as
///
///     (K/2) * (Z/A)*rho * x * (q²/beta²)
///
/// where (Z/A)*rho is the electron density in the material and x is the
/// traversed length (thickness) of the material.
inline float computeEpsilon(float molarElectronDensity, float thickness,
                            const Acts::RelativisticQuantities& rq) {
  return 0.5f * K * molarElectronDensity * thickness * rq.q2OverBeta2;
}

/// Compute epsilon logarithmic derivative w/ respect to q/p.
inline float logDeriveEpsilon(float qOverP,
                              const Acts::RelativisticQuantities& rq) {
  // only need to compute d(q²/beta²)/(q²/beta²); everything else cancels.
  return 2 / (qOverP * rq.gamma * rq.gamma);
}

/// Compute the density correction factor delta/2.
inline float computeDeltaHalf(float meanExitationPotential,
                              float molarElectronDensity,
                              const Acts::RelativisticQuantities& rq) {
  /// Uses RPP2018 eq. 33.6 which is only valid for high energies.
  // only relevant for very high ernergies; use arbitrary cutoff
  if (rq.betaGamma < 10.0f) {
    return 0.0f;
  }
  // pre-factor according to RPP2019 table 33.1
  const auto plasmaEnergy =
      PlasmaEnergyScale * std::sqrt(1000.f * molarElectronDensity);
  return std::log(rq.betaGamma) +
         std::log(plasmaEnergy / meanExitationPotential) - 0.5f;
}

/// Compute derivative w/ respect to q/p for the density correction.
inline float deriveDeltaHalf(float qOverP,
                             const Acts::RelativisticQuantities& rq) {
  // original equation is of the form
  //     log(beta*gamma) + log(eplasma/I) - 1/2
  // which the resulting derivative as
  //     d(beta*gamma)/(beta*gamma)
  return (rq.betaGamma < 10.0f) ? 0.0f : (-1.0f / qOverP);
}

}  // namespace

#define ASSERT_INPUTS(mass, qOverP, q)              \
  assert((0 < mass) and "Mass must be positive");   \
  assert((qOverP != 0) and "q/p must be non-zero"); \
  assert((q != 0) and "Charge must be non-zero");

float Acts::computeEnergyLossBethe(const MaterialSlab& slab, int /* unused */,
                                   float m, float qOverP, float q) {
  ASSERT_INPUTS(m, qOverP, q)

  // return early in case of vacuum or zero thickness
  if (not slab) {
    return 0.0f;
  }

  const auto I = slab.material().meanExcitationEnergy();
  const auto Ne = slab.material().molarElectronDensity();
  const auto thickness = slab.thickness();
  const auto rq = Acts::RelativisticQuantities(m, qOverP, q);
  const auto eps = computeEpsilon(Ne, thickness, rq);
  const auto dhalf = computeDeltaHalf(I, Ne, rq);
  const auto u = computeMassTerm(Me, rq);
  const auto wmax = computeWMax(m, rq);
  // uses RPP2018 eq. 33.5 scaled from mass stopping power to linear stopping
  // power and multiplied with the material thickness to get a total energy loss
  // instead of an energy loss per length.
  // the required modification only change the prefactor which becomes
  // identical to the prefactor epsilon for the most probable value.
  const auto running =
      std::log(u / I) + std::log(wmax / I) - 2.0f * rq.beta2 - 2.0f * dhalf;
  return eps * running;
}

float Acts::deriveEnergyLossBetheQOverP(const MaterialSlab& slab,
                                        int /* unused */, float m, float qOverP,
                                        float q) {
  ASSERT_INPUTS(m, qOverP, q)

  // return early in case of vacuum or zero thickness
  if (not slab) {
    return 0.0f;
  }

  const auto I = slab.material().meanExcitationEnergy();
  const auto Ne = slab.material().molarElectronDensity();
  const auto thickness = slab.thickness();
  const auto rq = Acts::RelativisticQuantities(m, qOverP, q);
  const auto eps = computeEpsilon(Ne, thickness, rq);
  const auto dhalf = computeDeltaHalf(I, Ne, rq);
  const auto u = computeMassTerm(Me, rq);
  const auto wmax = computeWMax(m, rq);
  // original equation is of the form
  //
  //     eps * (log(u/I) + log(wmax/I) - 2 * beta² - delta)
  //     = eps * (log(u) + log(wmax) - 2 * log(I) - 2 * beta² - delta)
  //
  // with the resulting derivative as
  //
  //     d(eps) * (log(u/I) + log(wmax/I) - 2 * beta² - delta)
  //     + eps * (d(u)/(u) + d(wmax)/(wmax) - 2 * d(beta²) - d(delta))
  //
  // where we can use d(eps) = eps * (d(eps)/eps) for further simplification.
  const auto logDerEps = logDeriveEpsilon(qOverP, rq);
  const auto derDHalf = deriveDeltaHalf(qOverP, rq);
  const auto logDerU = logDeriveMassTerm(qOverP);
  const auto logDerWmax = logDeriveWMax(m, qOverP, rq);
  const auto derBeta2 = deriveBeta2(qOverP, rq);
  const auto rel = logDerEps * (std::log(u / I) + std::log(wmax / I) -
                                2.0f * rq.beta2 - 2.0f * dhalf) +
                   logDerU + logDerWmax - 2.0f * derBeta2 - 2.0f * derDHalf;
  return eps * rel;
}

float Acts::computeEnergyLossLandau(const MaterialSlab& slab, int /* unused */,
                                    float m, float qOverP, float q) {
  ASSERT_INPUTS(m, qOverP, q)

  // return early in case of vacuum or zero thickness
  if (not slab) {
    return 0.0f;
  }

  const auto I = slab.material().meanExcitationEnergy();
  const auto Ne = slab.material().molarElectronDensity();
  const auto thickness = slab.thickness();
  const auto rq = Acts::RelativisticQuantities(m, qOverP, q);
  const auto eps = computeEpsilon(Ne, thickness, rq);
  const auto dhalf = computeDeltaHalf(I, Ne, rq);
  const auto t = computeMassTerm(Me, rq);
  // uses RPP2018 eq. 33.11
  const auto running =
      std::log(t / I) + std::log(eps / I) + 0.2f - rq.beta2 - 2 * dhalf;
  return eps * running;
}

float Acts::deriveEnergyLossLandauQOverP(const MaterialSlab& slab,
                                         int /* unused */, float m,
                                         float qOverP, float q) {
  ASSERT_INPUTS(m, qOverP, q)

  // return early in case of vacuum or zero thickness
  if (not slab) {
    return 0.0f;
  }

  const auto I = slab.material().meanExcitationEnergy();
  const auto Ne = slab.material().molarElectronDensity();
  const auto thickness = slab.thickness();
  const auto rq = Acts::RelativisticQuantities(m, qOverP, q);
  const auto eps = computeEpsilon(Ne, thickness, rq);
  const auto dhalf = computeDeltaHalf(I, Ne, rq);
  const auto t = computeMassTerm(Me, rq);
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

/// Convert Landau full-width-half-maximum to an equivalent Gaussian sigma,
///
/// Full-width-half-maximum for a Gaussian is given as
///
///     fwhm = 2 * sqrt(2 * log(2)) * sigma
/// -> sigma = fwhm / (2 * sqrt(2 * log(2)))
///
inline float convertLandauFwhmToGaussianSigma(float fwhm) {
  return fwhm / (2 * std::sqrt(2 * std::log(2.0f)));
}

}  // namespace

float Acts::computeEnergyLossLandauSigma(const MaterialSlab& slab,
                                         int /* unused */, float m,
                                         float qOverP, float q) {
  ASSERT_INPUTS(m, qOverP, q)

  // return early in case of vacuum or zero thickness
  if (not slab) {
    return 0.0f;
  }

  const auto Ne = slab.material().molarElectronDensity();
  const auto thickness = slab.thickness();
  const auto rq = Acts::RelativisticQuantities(m, qOverP, q);
  // the Landau-Vavilov fwhm is 4*eps (see RPP2018 fig. 33.7)
  const auto fwhm = 4 * computeEpsilon(Ne, thickness, rq);
  return convertLandauFwhmToGaussianSigma(fwhm);
}

float Acts::computeEnergyLossLandauFwhm(
    const MaterialSlab& slab, const Acts::RelativisticQuantities& rq) {
  // return early in case of vacuum or zero thickness
  if (not slab) {
    return 0.0f;
  }

  const auto Ne = slab.material().molarElectronDensity();
  const auto thickness = slab.thickness();
  // the Landau-Vavilov fwhm is 4*eps (see RPP2018 fig. 33.7)
  return 4 * computeEpsilon(Ne, thickness, rq);
}

float Acts::computeEnergyLossLandauSigmaQOverP(const MaterialSlab& slab,
                                               int /* unused */, float m,
                                               float qOverP, float q) {
  ASSERT_INPUTS(m, qOverP, q)

  const auto rq = Acts::RelativisticQuantities(m, qOverP, q);
  const auto fwhm = computeEnergyLossLandauFwhm(slab, rq);
  const auto sigmaE = convertLandauFwhmToGaussianSigma(fwhm);
  //  var(q/p) = (d(q/p)/dE)² * var(E)
  // d(q/p)/dE = d/dE (q/sqrt(E²-m²))
  //           = q * -(1/2) * 1/p³ * 2E
  //           = -q/p² E/p = -(q/p)² * 1/(q*beta) = -(q/p)² * (q/beta) / q²
  //  var(q/p) = (q/p)^4 * (q/beta)² * (1/q)^4 * var(E)
  //           = (1/p)^4 * (q/beta)² * var(E)
  // do not need to care about the sign since it is only used squared
  const auto pInv = qOverP / q;
  assert(rq.q2OverBeta2 >= 0 && "Negative q2OverBeta2");
  return clampValue<float>(std::sqrt(rq.q2OverBeta2) * pInv * pInv * sigmaE);
}

namespace {

/// Compute mean energy loss from bremsstrahlung per radiation length.
inline float computeBremsstrahlungLossMean(float mass, float energy) {
  return energy * (Me / mass) * (Me / mass);
}

/// Derivative of the bremsstrahlung loss per rad length with respect to energy.
inline float deriveBremsstrahlungLossMeanE(float mass) {
  return (Me / mass) * (Me / mass);
}

/// Expansion coefficients for the muon radiative loss as a function of energy
///
/// Taken from ATL-SOFT-PUB-2008-003 eq. 7,8 where the expansion is expressed
/// with terms E^n/X0 with fixed units [E] = MeV and [X0] = mm. The evaluated
/// expansion has units MeV/mm. In this implementation, the X0 dependence is
/// factored out and the coefficients must be scaled to the native units such
/// that the evaluated expansion with terms E^n has dimension energy in
/// native units.
constexpr float MuonHighLowThreshold = 1_TeV;
// [low0 / X0] = MeV / mm -> [low0] = MeV
constexpr double MuonLow0 = -0.5345_MeV;
// [low1 * E / X0] = MeV / mm -> [low1] = 1
constexpr double MuonLow1 = 6.803e-5;
// [low2 * E^2 / X0] = MeV / mm -> [low2] = 1/MeV
constexpr double MuonLow2 = 2.278e-11 / 1_MeV;
// [low3 * E^3 / X0] = MeV / mm -> [low3] = 1/MeV^2
constexpr double MuonLow3 = -9.899e-18 / (1_MeV * 1_MeV);
// units are the same as low0
constexpr double MuonHigh0 = -2.986_MeV;
// units are the same as low1
constexpr double MuonHigh1 = 9.253e-5;

/// Compute additional radiation energy loss for muons per radiation length.
inline float computeMuonDirectPairPhotoNuclearLossMean(double energy) {
  if (energy < MuonHighLowThreshold) {
    return MuonLow0 +
           (MuonLow1 + (MuonLow2 + MuonLow3 * energy) * energy) * energy;
  } else {
    return MuonHigh0 + MuonHigh1 * energy;
  }
}

/// Derivative of the additional rad loss per rad length with respect to energy.
inline float deriveMuonDirectPairPhotoNuclearLossMeanE(double energy) {
  if (energy < MuonHighLowThreshold) {
    return MuonLow1 + (2 * MuonLow2 + 3 * MuonLow3 * energy) * energy;
  } else {
    return MuonHigh1;
  }
}

}  // namespace

float Acts::computeEnergyLossRadiative(const MaterialSlab& slab, int pdg,
                                       float m, float qOverP, float q) {
  ASSERT_INPUTS(m, qOverP, q)

  // return early in case of vacuum or zero thickness
  if (not slab) {
    return 0.0f;
  }

  // relative radiation length
  const auto x = slab.thicknessInX0();
  // particle momentum and energy
  // do not need to care about the sign since it is only used squared
  const auto momentum = q / qOverP;
  const auto energy = std::hypot(m, momentum);

  auto dEdx = computeBremsstrahlungLossMean(m, energy);
  if (((pdg == PdgParticle::eMuon) or (pdg == PdgParticle::eAntiMuon)) and
      (8_GeV < energy)) {
    dEdx += computeMuonDirectPairPhotoNuclearLossMean(energy);
  }
  // scale from energy loss per unit radiation length to total energy
  return dEdx * x;
}

float Acts::deriveEnergyLossRadiativeQOverP(const MaterialSlab& slab, int pdg,
                                            float m, float qOverP, float q) {
  ASSERT_INPUTS(m, qOverP, q)

  // return early in case of vacuum or zero thickness
  if (not slab) {
    return 0.0f;
  }

  // relative radiation length
  const auto x = slab.thicknessInX0();
  // particle momentum and energy
  // do not need to care about the sign since it is only used squared
  const auto momentum = q / qOverP;
  const auto energy = std::hypot(m, momentum);

  // compute derivative w/ respect to energy.
  auto derE = deriveBremsstrahlungLossMeanE(m);
  if (((pdg == PdgParticle::eMuon) or (pdg == PdgParticle::eAntiMuon)) and
      (8_GeV < energy)) {
    derE += deriveMuonDirectPairPhotoNuclearLossMeanE(energy);
  }
  // compute derivative w/ respect to q/p by using the chain rule
  //     df(E)/d(q/p) = df(E)/dE dE/d(q/p)
  // with
  //     E = sqrt(m² + p²) = sqrt(m² + q²/(q/p)²)
  // and the resulting derivative
  //     dE/d(q/p) = -q² / ((q/p)³ * E)
  const auto derQOverP = -(q * q) / (qOverP * qOverP * qOverP * energy);
  return derE * derQOverP * x;
}

float Acts::computeEnergyLossMean(const MaterialSlab& slab, int pdg, float m,
                                  float qOverP, float q) {
  return computeEnergyLossBethe(slab, pdg, m, qOverP, q) +
         computeEnergyLossRadiative(slab, pdg, m, qOverP, q);
}

float Acts::deriveEnergyLossMeanQOverP(const MaterialSlab& slab, int pdg,
                                       float m, float qOverP, float q) {
  return deriveEnergyLossBetheQOverP(slab, pdg, m, qOverP, q) +
         deriveEnergyLossRadiativeQOverP(slab, pdg, m, qOverP, q);
}

float Acts::computeEnergyLossMode(const MaterialSlab& slab, int pdg, float m,
                                  float qOverP, float q) {
  // see ATL-SOFT-PUB-2008-003 section 3 for the relative fractions
  // TODO this is inconsistent with the text of the note
  return 0.9f * computeEnergyLossLandau(slab, pdg, m, qOverP, q) +
         0.15f * computeEnergyLossRadiative(slab, pdg, m, qOverP, q);
}

float Acts::deriveEnergyLossModeQOverP(const MaterialSlab& slab, int pdg,
                                       float m, float qOverP, float q) {
  // see ATL-SOFT-PUB-2008-003 section 3 for the relative fractions
  // TODO this is inconsistent with the text of the note
  return 0.9f * deriveEnergyLossLandauQOverP(slab, pdg, m, qOverP, q) +
         0.15f * deriveEnergyLossRadiativeQOverP(slab, pdg, m, qOverP, q);
}

namespace {

/// Multiple scattering theta0 for minimum ionizing particles.
inline float theta0Highland(float xOverX0, float momentumInv,
                            float q2OverBeta2) {
  // RPP2018 eq. 33.15 (treats beta and q² consistenly)
  const auto t = std::sqrt(xOverX0 * q2OverBeta2);
  // log((x/X0) * (q²/beta²)) = log((sqrt(x/X0) * (q/beta))²)
  //                          = 2 * log(sqrt(x/X0) * (q/beta))
  return 13.6_MeV * momentumInv * t * (1.0f + 0.038f * 2 * std::log(t));
}

/// Multiple scattering theta0 for electrons.
inline float theta0RossiGreisen(float xOverX0, float momentumInv,
                                float q2OverBeta2) {
  // TODO add source paper/ resource
  const auto t = std::sqrt(xOverX0 * q2OverBeta2);
  return 17.5_MeV * momentumInv * t *
         (1.0f + 0.125f * std::log10(10.0f * xOverX0));
}

}  // namespace

float Acts::computeMultipleScatteringTheta0(const MaterialSlab& slab, int pdg,
                                            float m, float qOverP, float q) {
  ASSERT_INPUTS(m, qOverP, q)

  // return early in case of vacuum or zero thickness
  if (not slab) {
    return 0.0f;
  }

  // relative radiation length
  const auto xOverX0 = slab.thicknessInX0();
  // 1/p = q/(pq) = (q/p)/q
  const auto momentumInv = std::abs(qOverP / q);
  // q²/beta²; a smart compiler should be able to remove the unused computations
  const auto q2OverBeta2 =
      Acts::RelativisticQuantities(m, qOverP, q).q2OverBeta2;

  if ((pdg == PdgParticle::eElectron) or (pdg == PdgParticle::ePositron)) {
    return theta0RossiGreisen(xOverX0, momentumInv, q2OverBeta2);
  } else {
    return theta0Highland(xOverX0, momentumInv, q2OverBeta2);
  }
}
