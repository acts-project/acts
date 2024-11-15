// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/Interactions.hpp"

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <cassert>
#include <cmath>

using namespace Acts::UnitLiterals;

namespace {

// values from RPP2018 table 33.1
// electron mass
constexpr double Me = 0.5109989461_MeV;
// Bethe formular prefactor. 1/mol unit is just a factor 1 here.
constexpr double K = 0.307075_MeV * 1_cm * 1_cm;
// Energy scale for plasma energy.
constexpr double PlasmaEnergyScale = 28.816_eV;

/// Additional derived relativistic quantities.
struct RelativisticQuantities {
  double q2OverBeta2 = 0;
  double beta2 = 0;
  double betaGamma = 0;
  double gamma = 0;

  RelativisticQuantities(double mass, double qOverP, double absQ) {
    assert((0 < mass) && "Mass must be positive");
    assert((qOverP != 0) && "q/p must be non-zero");
    assert((absQ > 0) && "Absolute charge must be non-zero and positive");
    // beta²/q² = (p/E)²/q² = p²/(q²m² + q²p²) = 1/(q² + (m²(q/p)²)
    // q²/beta² = q² + m²(q/p)²
    q2OverBeta2 = absQ * absQ + (mass * qOverP) * (mass * qOverP);
    assert((q2OverBeta2 >= 0) && "Negative q2OverBeta2");
    // 1/p = q/(qp) = (q/p)/q
    const double mOverP = mass * std::abs(qOverP / absQ);
    const double pOverM = 1 / mOverP;
    // beta² = p²/E² = p²/(m² + p²) = 1/(1 + (m/p)²)
    beta2 = 1 / (1 + mOverP * mOverP);
    assert((beta2 >= 0) && "Negative beta2");
    // beta*gamma = (p/sqrt(m² + p²))*(sqrt(m² + p²)/m) = p/m
    betaGamma = pOverM;
    assert((betaGamma >= 0) && "Negative betaGamma");
    // gamma = sqrt(m² + p²)/m = sqrt(1 + (p/m)²)
    gamma = Acts::fastHypot(1.0f, pOverM);
  }
};

/// Compute q/p derivative of beta².
inline double deriveBeta2(double qOverP, const RelativisticQuantities& rq) {
  return -2 / (qOverP * rq.gamma * rq.gamma);
}

/// Compute the 2 * mass * (beta * gamma)² mass term.
inline double computeMassTerm(double mass, const RelativisticQuantities& rq) {
  return 2 * mass * rq.betaGamma * rq.betaGamma;
}

/// Compute mass term logarithmic derivative w/ respect to q/p.
inline double logDeriveMassTerm(double qOverP) {
  // only need to compute d((beta*gamma)²)/(beta*gamma)²; rest cancels.
  return -2 / qOverP;
}

/// Compute the maximum energy transfer in a single collision.
///
/// Uses RPP2018 eq. 33.4.
inline double computeWMax(double mass, const RelativisticQuantities& rq) {
  const double mfrac = Me / mass;
  const double nominator = 2 * Me * rq.betaGamma * rq.betaGamma;
  const double denonimator = 1 + 2 * rq.gamma * mfrac + mfrac * mfrac;
  return nominator / denonimator;
}

/// Compute WMax logarithmic derivative w/ respect to q/p.
inline double logDeriveWMax(double mass, double qOverP,
                            const RelativisticQuantities& rq) {
  // this is (q/p) * (beta/q).
  // both quantities have the same sign and the product must always be
  // positive. we can thus reuse the known (unsigned) quantity (q/beta)².
  const double a = std::abs(qOverP / std::sqrt(rq.q2OverBeta2));
  // (m² + me²) / me = me (1 + (m/me)²)
  const double b = Me * (1 + (mass / Me) * (mass / Me));
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
inline double computeEpsilon(double molarElectronDensity, double thickness,
                             const RelativisticQuantities& rq) {
  return 0.5 * K * molarElectronDensity * thickness * rq.q2OverBeta2;
}

/// Compute epsilon logarithmic derivative w/ respect to q/p.
inline double logDeriveEpsilon(double qOverP,
                               const RelativisticQuantities& rq) {
  // only need to compute d(q²/beta²)/(q²/beta²); everything else cancels.
  return 2 / (qOverP * rq.gamma * rq.gamma);
}

/// Compute the density correction factor delta/2.
inline double computeDeltaHalf(double meanExitationPotential,
                               double molarElectronDensity,
                               const RelativisticQuantities& rq) {
  /// Uses RPP2018 eq. 33.6 which is only valid for high energies.
  // only relevant for very high ernergies; use arbitrary cutoff
  if (rq.betaGamma < 10.0) {
    return 0;
  }
  // pre-factor according to RPP2019 table 33.1
  const double plasmaEnergy =
      PlasmaEnergyScale * std::sqrt(1000 * molarElectronDensity);
  return std::log(rq.betaGamma) +
         std::log(plasmaEnergy / meanExitationPotential) - 0.5;
}

/// Compute derivative w/ respect to q/p for the density correction.
inline double deriveDeltaHalf(double qOverP, const RelativisticQuantities& rq) {
  // original equation is of the form
  //     log(beta*gamma) + log(eplasma/I) - 1/2
  // which the resulting derivative as
  //     d(beta*gamma)/(beta*gamma)
  return (rq.betaGamma < 10.0) ? 0 : (-1 / qOverP);
}

/// Convert Landau full-width-half-maximum to an equivalent Gaussian sigma,
///
/// Full-width-half-maximum for a Gaussian is given as
///
///     fwhm = 2 * sqrt(2 * log(2)) * sigma
/// -> sigma = fwhm / (2 * sqrt(2 * log(2)))
///
inline double convertLandauFwhmToGaussianSigma(double fwhm) {
  // return fwhm / (2 * std::sqrt(2 * std::log(2.0)));
  return fwhm * 0.42466090014400953;
}

namespace detail {

inline double computeEnergyLossLandauFwhm(const Acts::MaterialSlab& slab,
                                          const RelativisticQuantities& rq) {
  // return early in case of vacuum or zero thickness
  if (!slab.isValid()) {
    return 0;
  }

  const double Ne = slab.material().molarElectronDensity();
  const double thickness = slab.thickness();
  // the Landau-Vavilov fwhm is 4*eps (see RPP2018 fig. 33.7)
  return 4 * computeEpsilon(Ne, thickness, rq);
}

}  // namespace detail

}  // namespace

double Acts::computeEnergyLossBethe(const MaterialSlab& slab, double m,
                                    double qOverP, double absQ) {
  // return early in case of vacuum or zero thickness
  if (!slab.isValid()) {
    return 0;
  }

  const RelativisticQuantities rq{m, qOverP, absQ};
  const double I = slab.material().meanExcitationEnergy();
  const double Ne = slab.material().molarElectronDensity();
  const double thickness = slab.thickness();
  const double eps = computeEpsilon(Ne, thickness, rq);
  const double dhalf = computeDeltaHalf(I, Ne, rq);
  const double u = computeMassTerm(Me, rq);
  const double wmax = computeWMax(m, rq);
  // uses RPP2018 eq. 33.5 scaled from mass stopping power to linear stopping
  // power and multiplied with the material thickness to get a total energy loss
  // instead of an energy loss per length.
  // the required modification only change the prefactor which becomes
  // identical to the prefactor epsilon for the most probable value.
  const double running =
      std::log(u / I) + std::log(wmax / I) - 2 * rq.beta2 - 2 * dhalf;
  return eps * running;
}

double Acts::deriveEnergyLossBetheQOverP(const MaterialSlab& slab, double m,
                                         double qOverP, double absQ) {
  // return early in case of vacuum or zero thickness
  if (!slab.isValid()) {
    return 0;
  }

  const RelativisticQuantities rq{m, qOverP, absQ};
  const double I = slab.material().meanExcitationEnergy();
  const double Ne = slab.material().molarElectronDensity();
  const double thickness = slab.thickness();
  const double eps = computeEpsilon(Ne, thickness, rq);
  const double dhalf = computeDeltaHalf(I, Ne, rq);
  const double u = computeMassTerm(Me, rq);
  const double wmax = computeWMax(m, rq);
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
  const double logDerEps = logDeriveEpsilon(qOverP, rq);
  const double derDHalf = deriveDeltaHalf(qOverP, rq);
  const double logDerU = logDeriveMassTerm(qOverP);
  const double logDerWmax = logDeriveWMax(m, qOverP, rq);
  const double derBeta2 = deriveBeta2(qOverP, rq);
  const double rel = logDerEps * (std::log(u / I) + std::log(wmax / I) -
                                  2 * rq.beta2 - 2 * dhalf) +
                     logDerU + logDerWmax - 2 * derBeta2 - 2 * derDHalf;
  return eps * rel;
}

double Acts::computeEnergyLossLandau(const MaterialSlab& slab, double m,
                                     double qOverP, double absQ) {
  // return early in case of vacuum or zero thickness
  if (!slab.isValid()) {
    return 0;
  }

  const RelativisticQuantities rq{m, qOverP, absQ};
  const double I = slab.material().meanExcitationEnergy();
  const double Ne = slab.material().molarElectronDensity();
  const double thickness = slab.thickness();
  const double eps = computeEpsilon(Ne, thickness, rq);
  const double dhalf = computeDeltaHalf(I, Ne, rq);
  const double u = computeMassTerm(Me, rq);
  // uses RPP2018 eq. 33.12
  const double running =
      std::log(u / I) + std::log(eps / I) + 0.2 - rq.beta2 - 2 * dhalf;
  return eps * running;
}

double Acts::deriveEnergyLossLandauQOverP(const MaterialSlab& slab, double m,
                                          double qOverP, double absQ) {
  // return early in case of vacuum or zero thickness
  if (!slab.isValid()) {
    return 0;
  }

  const RelativisticQuantities rq{m, qOverP, absQ};
  const double I = slab.material().meanExcitationEnergy();
  const double Ne = slab.material().molarElectronDensity();
  const double thickness = slab.thickness();
  const double eps = computeEpsilon(Ne, thickness, rq);
  const double dhalf = computeDeltaHalf(I, Ne, rq);
  const double t = computeMassTerm(Me, rq);
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
  const double logDerEps = logDeriveEpsilon(qOverP, rq);
  const double derDHalf = deriveDeltaHalf(qOverP, rq);
  const double logDerT = logDeriveMassTerm(qOverP);
  const double derBeta2 = deriveBeta2(qOverP, rq);
  const double rel = logDerEps * (std::log(t / I) + std::log(eps / I) - 0.2 -
                                  rq.beta2 - 2 * dhalf) +
                     logDerT + logDerEps - derBeta2 - 2 * derDHalf;
  return eps * rel;
}

double Acts::computeEnergyLossLandauSigma(const MaterialSlab& slab, double m,
                                          double qOverP, double absQ) {
  // return early in case of vacuum or zero thickness
  if (!slab.isValid()) {
    return 0;
  }

  const RelativisticQuantities rq{m, qOverP, absQ};
  const double Ne = slab.material().molarElectronDensity();
  const double thickness = slab.thickness();
  // the Landau-Vavilov fwhm is 4*eps (see RPP2018 fig. 33.7)
  const double fwhm = 4 * computeEpsilon(Ne, thickness, rq);
  return convertLandauFwhmToGaussianSigma(fwhm);
}

double Acts::computeEnergyLossLandauFwhm(const MaterialSlab& slab, double m,
                                         double qOverP, double absQ) {
  const RelativisticQuantities rq{m, qOverP, absQ};
  return detail::computeEnergyLossLandauFwhm(slab, rq);
}

double Acts::computeEnergyLossLandauSigmaQOverP(const MaterialSlab& slab,
                                                double m, double qOverP,
                                                double absQ) {
  const RelativisticQuantities rq{m, qOverP, absQ};
  const double fwhm = detail::computeEnergyLossLandauFwhm(slab, rq);
  const double sigmaE = convertLandauFwhmToGaussianSigma(fwhm);
  //  var(q/p) = (d(q/p)/dE)² * var(E)
  // d(q/p)/dE = d/dE (q/sqrt(E²-m²))
  //           = q * -(1/2) * 1/p³ * 2E
  //           = -q/p² E/p = -(q/p)² * 1/(q*beta) = -(q/p)² * (q/beta) / q²
  //  var(q/p) = (q/p)^4 * (q/beta)² * (1/q)^4 * var(E)
  //           = (1/p)^4 * (q/beta)² * var(E)
  // do not need to care about the sign since it is only used squared
  const double pInv = qOverP / absQ;
  const double qOverBeta = std::sqrt(rq.q2OverBeta2);
  return qOverBeta * pInv * pInv * sigmaE;
}

namespace {

/// Compute mean energy loss from bremsstrahlung per radiation length.
inline double computeBremsstrahlungLossMean(double mass, double energy) {
  return energy * (Me / mass) * (Me / mass);
}

/// Derivative of the bremsstrahlung loss per rad length with respect to energy.
inline double deriveBremsstrahlungLossMeanE(double mass) {
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
constexpr double MuonHighLowThreshold = 1_TeV;
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
inline double computeMuonDirectPairPhotoNuclearLossMean(double energy) {
  if (energy < MuonHighLowThreshold) {
    return MuonLow0 +
           (MuonLow1 + (MuonLow2 + MuonLow3 * energy) * energy) * energy;
  } else {
    return MuonHigh0 + MuonHigh1 * energy;
  }
}

/// Derivative of the additional rad loss per rad length with respect to energy.
inline double deriveMuonDirectPairPhotoNuclearLossMeanE(double energy) {
  if (energy < MuonHighLowThreshold) {
    return MuonLow1 + (2 * MuonLow2 + 3 * MuonLow3 * energy) * energy;
  } else {
    return MuonHigh1;
  }
}

}  // namespace

double Acts::computeEnergyLossRadiative(const MaterialSlab& slab,
                                        PdgParticle absPdg, double m,
                                        double qOverP, double absQ) {
  assert((absPdg == Acts::makeAbsolutePdgParticle(absPdg)) &&
         "pdg is not absolute");

  // return early in case of vacuum or zero thickness
  if (!slab.isValid()) {
    return 0;
  }

  // relative radiation length
  const double x = slab.thicknessInX0();
  // particle momentum and energy
  // do not need to care about the sign since it is only used squared
  const double momentum = absQ / qOverP;
  const double energy = fastHypot(m, momentum);

  double dEdx = computeBremsstrahlungLossMean(m, energy);

  // muon- or muon+
  // TODO magic number 8_GeV
  if ((absPdg == PdgParticle::eMuon) && (8_GeV < energy)) {
    dEdx += computeMuonDirectPairPhotoNuclearLossMean(energy);
  }
  // scale from energy loss per unit radiation length to total energy
  return dEdx * x;
}

double Acts::deriveEnergyLossRadiativeQOverP(const MaterialSlab& slab,
                                             PdgParticle absPdg, double m,
                                             double qOverP, double absQ) {
  assert((absPdg == Acts::makeAbsolutePdgParticle(absPdg)) &&
         "pdg is not absolute");

  // return early in case of vacuum or zero thickness
  if (!slab.isValid()) {
    return 0;
  }

  // relative radiation length
  const double x = slab.thicknessInX0();
  // particle momentum and energy
  // do not need to care about the sign since it is only used squared
  const double momentum = absQ / qOverP;
  const double energy = fastHypot(m, momentum);

  // compute derivative w/ respect to energy.
  double derE = deriveBremsstrahlungLossMeanE(m);

  // muon- or muon+
  // TODO magic number 8_GeV
  if ((absPdg == PdgParticle::eMuon) && (8_GeV < energy)) {
    derE += deriveMuonDirectPairPhotoNuclearLossMeanE(energy);
  }
  // compute derivative w/ respect to q/p by using the chain rule
  //     df(E)/d(q/p) = df(E)/dE dE/d(q/p)
  // with
  //     E = sqrt(m² + p²) = sqrt(m² + q²/(q/p)²)
  // and the resulting derivative
  //     dE/d(q/p) = -q² / ((q/p)³ * E)
  const double derQOverP = -(absQ * absQ) / (qOverP * qOverP * qOverP * energy);
  return derE * derQOverP * x;
}

double Acts::computeEnergyLossMean(const MaterialSlab& slab, PdgParticle absPdg,
                                   double m, double qOverP, double absQ) {
  return computeEnergyLossBethe(slab, m, qOverP, absQ) +
         computeEnergyLossRadiative(slab, absPdg, m, qOverP, absQ);
}

double Acts::deriveEnergyLossMeanQOverP(const MaterialSlab& slab,
                                        PdgParticle absPdg, double m,
                                        double qOverP, double absQ) {
  return deriveEnergyLossBetheQOverP(slab, m, qOverP, absQ) +
         deriveEnergyLossRadiativeQOverP(slab, absPdg, m, qOverP, absQ);
}

double Acts::computeEnergyLossMode(const MaterialSlab& slab, PdgParticle absPdg,
                                   double m, double qOverP, double absQ) {
  // see ATL-SOFT-PUB-2008-003 section 3 for the relative fractions
  // TODO this is inconsistent with the text of the note
  return 0.9 * computeEnergyLossLandau(slab, m, qOverP, absQ) +
         0.15 * computeEnergyLossRadiative(slab, absPdg, m, qOverP, absQ);
}

double Acts::deriveEnergyLossModeQOverP(const MaterialSlab& slab,
                                        PdgParticle absPdg, double m,
                                        double qOverP, double absQ) {
  // see ATL-SOFT-PUB-2008-003 section 3 for the relative fractions
  // TODO this is inconsistent with the text of the note
  return 0.9 * deriveEnergyLossLandauQOverP(slab, m, qOverP, absQ) +
         0.15 * deriveEnergyLossRadiativeQOverP(slab, absPdg, m, qOverP, absQ);
}

namespace {

/// Multiple scattering theta0 for minimum ionizing particles.
inline double theta0Highland(double xOverX0, double momentumInv,
                             double q2OverBeta2) {
  // RPP2018 eq. 33.15 (treats beta and q² consistently)
  const double t = std::sqrt(xOverX0 * q2OverBeta2);
  // log((x/X0) * (q²/beta²)) = log((sqrt(x/X0) * (q/beta))²)
  //                          = 2 * log(sqrt(x/X0) * (q/beta))
  return 13.6_MeV * momentumInv * t * (1 + 0.038 * 2 * std::log(t));
}

/// Multiple scattering theta0 for electrons.
inline double theta0RossiGreisen(double xOverX0, double momentumInv,
                                 double q2OverBeta2) {
  // TODO add source paper/ resource
  const double t = std::sqrt(xOverX0 * q2OverBeta2);
  return 17.5_MeV * momentumInv * t * (1 + 0.125 * std::log10(10.0 * xOverX0));
}

}  // namespace

double Acts::computeMultipleScatteringTheta0(const MaterialSlab& slab,
                                             PdgParticle absPdg, double m,
                                             double qOverP, double absQ) {
  assert((absPdg == Acts::makeAbsolutePdgParticle(absPdg)) &&
         "pdg is not absolute");

  // return early in case of vacuum or zero thickness
  if (!slab.isValid()) {
    return 0;
  }

  // relative radiation length
  const double xOverX0 = slab.thicknessInX0();
  // 1/p = q/(pq) = (q/p)/q
  const double momentumInv = std::abs(qOverP / absQ);
  // q²/beta²; a smart compiler should be able to remove the unused computations
  const double q2OverBeta2 =
      RelativisticQuantities(m, qOverP, absQ).q2OverBeta2;

  // electron or positron
  if (absPdg == PdgParticle::eElectron) {
    return theta0RossiGreisen(xOverX0, momentumInv, q2OverBeta2);
  } else {
    return theta0Highland(xOverX0, momentumInv, q2OverBeta2);
  }
}
