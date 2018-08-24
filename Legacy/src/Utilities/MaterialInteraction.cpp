// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/MaterialInteraction.hpp"

namespace Acts {

static std::pair<double, double>
ionizationEnergyLoss(bool            mean,
                     double          p,
                     double          m,
                     const Material& mat,
                     double          path)
{
  // the return value
  double dE = 0.;
  // kinetic variables
  // and the electron mass in MeV

  double me    = particleMasses.mass[electron];
  double mfrac = me / m;
  double E     = std::sqrt(p * p + m * m);
  double beta  = p / E;
  double gamma = E / m;

  // Ionization - Bethe-Bloch
  // See ATL-SOFT-PUB-2008-003 equation (4)
  // 16 eV * Z**0.9 - bring to MeV
  double I = 16 * units::_eV * std::pow(mat.Z(), 0.9);

  // See (1) table 32.1
  // K/A*Z = 0.5 * 30.7075MeV/(g/mm2) * Z/A * rho[g/mm3]  / scale to mm by this
  double kaz = 0.5 * 30.7075 * (units::_MeV * units::_mm * units::_mm)
      * mat.zOverAtimesRho();
  double eta2 = beta * gamma;
  eta2 *= eta2;
  // density effect, only valid for high energies (gamma > 10 -> p > 1GeV for
  // muons)
  double delta = 0.;
  if (gamma > 10.) {
    // See (1) table 32.1
    double eplasma
        = 28.816 * units::_eV * std::sqrt(1000. * mat.zOverAtimesRho());
    // See (1) formula 32.6
    delta = 2. * std::log(eplasma / I) + std::log(eta2) - 1.;
  }

  //@todo implement energy loss for electrons

  // divide by beta^2 for non-electrons
  kaz /= beta * beta;
  double kazL = kaz * path;

  // The landau width (FWHM) is 4.*kazL
  // The factor is the conversion factor from FWHM to sigma for
  // gaussian curve: 1. / (2. * std::sqrt(2. * std::log(2.))).
  double sigma = 2. * kazL * 1. / (sqrt(2. * std::log(2.)));

  if (mean) {
    // Calculate the mean value for reconstruction
    // See ATL-SOFT-PUB-2008-003 equation (2)
    double tMax = 2. * eta2 * me / (1. + 2. * gamma * mfrac + mfrac * mfrac);
    // See ATL-SOFT-PUB-2008-003 equation (1)
    // or
    // http://pdg.lbl.gov/2014/reviews/rpp2014-rev-passage-particles-matter.pdf
    // PDG formula 32.5
    dE = -kaz * 0.5
        * (log(2. * me * eta2 * tMax / (I * I)) - (beta * beta) - delta * 0.5);
    dE *= path;
  } else {
    // Calculate the most probably value for simulation
    //
    // the landau sigmaL is path length dependent
    //    PDG formula 32.11 for MOP value from
    //    http://pdg.lbl.gov/2014/reviews/rpp2014-rev-passage-particles-matter.pdf
    //
    dE = kazL
        * (log(2. * m * eta2 / I) + std::log(kazL / I) + 0.2 - (beta * beta)
           - delta);
  }

  return std::make_pair(dE, sigma);
}

std::pair<double, double>
ionizationEnergyLossMean(double          p,
                         const Material& mat,
                         ParticleType    particle,
                         double          path)
{
  return ionizationEnergyLoss(
      true, p, particleMasses.mass[particle], mat, path);
}

std::pair<double, double>
ionizationEnergyLossMpv(double                p,
                        const Material&       mat,
                        ParticleType          particle,
                        const ParticleMasses& pMasses,
                        double                path)
{
  return ionizationEnergyLoss(false, p, pMasses.mass[particle], mat, path);
}

std::pair<double, double>
ionizationEnergyLossMpv(double p, double m, const Material& mat, double path)
{
  return ionizationEnergyLoss(false, p, m, mat, path);
}

std::pair<double, double>
radiationEnergyLoss(double p, const Material& mat, ParticleType particle)
{
  double sigma = 0.;
  if (!(mat)) {
    return std::pair<double, double>(0., 0.);
  }

  // preparation of kinetic constants
  double m     = particleMasses.mass[particle];
  double me    = particleMasses.mass[electron];
  double mfrac = me / m;
  double E     = std::sqrt(p * p + m * m);

  // Bremsstrahlung - Bethe-Heitler
  // See also ATL-SOFT-PUB-2008-003 equation (6)
  double Radiation = -E * mfrac * mfrac;
  // sigma is rms of steep exponential part of radiation
  sigma = -Radiation;

  // Add e+e- pair production and photonuclear effect for muons at energies
  // above 8 GeV
  //    Radiation gives mean Eloss including the long tail from 'catastrophic'
  //    Eloss
  //    sigma the rms of steep exponential
  /// @todo Units?
  // See also ATL-SOFT-PUB-2008-003 equation (7)(8)
  if ((particle == muon) && (E > 8000.)) {
    if (E < 1.e6) {
      Radiation += 0.5345 - 6.803e-5 * E - 2.278e-11 * E * E
          + 9.899e-18 * E * E * E;  // E below 1 TeV
      sigma += (0.1828 - 3.966e-3 * std::sqrt(E) + 2.151e-5 * E);  // idem
    } else {
      Radiation += 2.986 - 9.253e-5 * E;           // E above 1 TeV
      sigma += 17.73 + 2.409e-5 * (E - 1000000.);  // idem
    }
  }

  sigma = sigma / mat.X0();

  return std::pair<double, double>(Radiation / mat.X0(), sigma);
}

double
sigmaMS(double dInX0, double p, double beta)
{
  if (dInX0 == 0. || p == 0. || beta == 0.) {
    return 0.;
  }

  // Highland formula - projected sigma_s
  // ATL-SOFT-PUB-2008-003 equation (15)
  double sig_ms = 13.6 * units::_MeV * std::sqrt(dInX0) / (beta * p)
      * (1. + 0.038 * std::log(dInX0 / (beta * beta)));
  return sig_ms;
}
}
