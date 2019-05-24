// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/Material.hpp"
#include "Acts/Propagator/detail/Constants.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {
namespace detail {

/// The ionization energy loss along a given path length.
///
/// The mean energy loss should be used for reconstruction,
/// the most probable energy loss is more suited for fast simulation
///
/// The energy loss is calculated using the following paper
/// http://pdg.lbl.gov/2018/reviews/rpp2018-rev-passage-particles-matter.pdf
///
/// Formula 33.5 is used to calculate the mean energy loss [reco mode]
/// Formula 33.11 is used to calculate the most probable energy loss [sim
/// mode]
/// Sigma is calculated using the Landau width (FHWM) times the conversion
/// factor 1. / (2. * &radic(2. * log2))
///
struct IonisationLoss {
  /// @brief call operator for the Ionisation class
  ///
  /// @tparam material_t Type of the material class
  /// @param [in] m      Mass of the particle
  /// @param [in] lbeta  Beta factor of the particle
  /// @param [in] lgamma Gamma factor of the particle
  /// @param [in] mat    The traversed material
  /// @param [in] path   The path length (optional)
  /// @param [in] mean   Toggle between mean (true) and mode (false)
  /// @return A std::pair. The first entry is the mean energy loss due to
  /// ionization along a given path length or its mode. The second entry is
  /// the sigma of the distribution.
  template <typename material_t>
  std::pair<double, double> dEds(double m, double lbeta, double lgamma,
                                 const material_t& mat, double path = 1.,
                                 bool mean = true) const {
    // the return value
    double dE = 0.;

    // Ionization - Bethe-Bloch
    // See ATL-SOFT-PUB-2008-003 equation (4)
    // 16 eV * Z**0.9
    double I = constants::eionisation * std::pow(mat.Z(), 0.9);

    // See (1) table 33.1
    // K/A*Z = 0.5 * 30.7075MeV/(g/mm2) * Z/A * rho[g/mm3]
    double kaz = 0.5 * constants::ka_BetheBloch * mat.zOverAtimesRho();
    double eta2 = lbeta * lgamma;
    eta2 *= eta2;
    // density effect, only valid for high energies
    // (lgamma > 10 -> p > 1 GeV for muons)
    double delta = 0.;
    if (lbeta * lgamma > 10.) {
      // See (1) table 33.1
      double eplasma = constants::eplasma * sqrt(1000. * mat.zOverAtimesRho());
      // See (1) formula 33.6
      delta = 2. * std::log(eplasma / I) + std::log(eta2) - 1.;
    }

    // divide by lbeta^2 for non-electrons
    kaz /= lbeta * lbeta;
    double kazL = kaz * path;

    // The landau width (FWHM) is 4.*kazL
    // The factor is the conversion factor from FWHM to sigma for
    // gaussian curve: 1. / (2. * sqrt(2. * log(2.))).
    double sigma = kazL * constants::landau2gauss;
    if (mean) {
      // Calculate the fraction to the electron mass
      double mfrac = constants::me / m;
      // Calculate the mean value for reconstruction
      // See ATL-SOFT-PUB-2008-003 equation (2)
      double tMax = 2. * eta2 * constants::me /
                    (1. + 2. * lgamma * mfrac + mfrac * mfrac);
      // See ATL-SOFT-PUB-2008-003 equation (1)
      // or:
      // http://pdg.lbl.gov/2018/reviews/rpp2018-rev-passage-particles-matter.pdf
      // PDG formula 33.5
      dE = -kazL * 2.0 *
           (0.5 * std::log(2. * constants::me * eta2 * tMax / (I * I)) -
            (lbeta * lbeta) - delta * 0.5);
    } else {
      // Calculate the most probably value for simulation
      //
      // the landau sigmaL is path length dependent
      //    PDG formula 33.11 for MOP value from
      //    http://pdg.lbl.gov/2018/reviews/rpp2018-rev-passage-particles-matter.pdf
      //
      dE = kazL * (std::log(2. * m * eta2 / I) + std::log(kazL / I) + 0.2 -
                   (lbeta * lbeta) - delta);
    }
    // return the energy loss and stragling
    return std::make_pair(dE, sigma);
  }

  /// @brief Evaluation of the energy loss dEds by ionisation derived by q/p
  /// (=d(dE/ds)/d(q/p)). The underlying equations for dE/ds are given by
  /// Formula 33.5 (mean energy loss) and Formula 33.11 (most probable energy
  /// loss).
  ///
  /// @tparam material_t   Type of the material
  /// @param [in] energy   Energy of the particle
  /// @param [in] qop      Charge over momentum of the particle
  /// @param [in] mass     Mass of the particle
  /// @param [in] material Material that is penetrated
  /// @param [in] mean     Toggle between mean (true) and mode (false)
  /// @return The evaluated derivative
  template <typename material_t>
  double dqop(const double energy, const double qOverP, const double mass,
              const material_t& material, const bool mean = true) const {
    // Fast exit if material is invalid
    if (material.Z() == 0 || material.zOverAtimesRho() == 0) {
      return 0.;
    }

    // Constants for readability
    const double m2 = mass * mass;
    const double I = constants::eionisation * std::pow(material.Z(), 0.9);
    const double gamma = energy / mass;
    const double beta = std::abs(1 / (energy * qOverP));
    const double beta2 = beta * beta;
    const double kaz =
        0.5 * constants::ka_BetheBloch * material.zOverAtimesRho();
    // Storage of the result
    double betheBlochDerivative = 0.;
    if (mean) {
      // Constants related to the mean evaluation
      const double me = constants::me;
      const double me2 = me * me;
      const double qop3 = qOverP * qOverP * qOverP;
      const double qop4 = qop3 * qOverP;
      const double m4 = m2 * m2;
      const double I2 = I * I;

      // Parts of the derivative
      const double lnCore = 4. * me2 / (m4 * I2 * qop4) /
                            (1. + 2. * gamma * me / mass + me2 / m2);
      const double lnCoreDerivative =
          -4. * me2 / (m4 * I2) *
          std::pow(qop4 + 2. * gamma * qop4 * me / mass + qop4 * me2 / m2,
                   -2.) *
          (4. * qop3 + 8. * me * qop3 * gamma / mass -
           2. * me * qOverP / (m2 * mass * gamma) + 4. * qop3 * me2 / m2);

      // Combine parts
      const double lnDerivative = 2. * qOverP * m2 * std::log(lnCore) +
                                  lnCoreDerivative / (lnCore * beta2);
      betheBlochDerivative = -kaz * lnDerivative;
    } else {
      // Parts of the mode derivative
      const double lnCore = std::log(2. * mass * beta2 * gamma * gamma / I);
      const double lnCoreXi = std::log(kaz / (beta2 * I));
      const double lnCoreDerivative = -2. * mass * beta * gamma;
      const double lnCoreXiDerivative = 2. * beta2 * m2;

      // Result evaluation (without the density effect)
      betheBlochDerivative =
          kaz * (2. * m2 * qOverP * (lnCore + lnCoreXi + 0.2) +
                 (lnCoreDerivative + lnCoreXiDerivative) / beta2);
    }
    // Density effect, only valid for high energies (gamma > 10 -> p > 1GeV
    // for muons). The derivative includes the kinematic parts of the
    // pre-factor
    if (gamma > 10.) {
      const double delta =
          2. * std::log(constants::eplasma *
                        std::sqrt(1000. * material.zOverAtimesRho()) / I) +
          2. * std::log(beta * gamma) - 1.;
      if (mean) {
        const double deltaDerivative =
            -2. / (qOverP * beta2) + 2. * delta * qOverP * m2;
        betheBlochDerivative += kaz * deltaDerivative;
      } else {
        const double deltaDerivative =
            -kaz * 2. / (qOverP * beta2) + kaz * 2. * m2 * qOverP * delta;
        betheBlochDerivative -= deltaDerivative;
      }
    }
    return betheBlochDerivative;
  }
};

/// @brief Multiple scattering as function of dInX0
///
/// It supports MIP and electron scattering and return
/// the space angle which has to be transformed into actual
/// scattering components for the azimuthal and polar angles,
/// respectively (for reconstruction).
struct HighlandScattering {
  /// @brief call operator for the HighlandScattering formula
  ///
  /// @param p is the momentum
  /// @param lbeta is the Lorentz beta parameter
  /// @param dInX0 is the passed thickness expressed in units of X0
  double operator()(double p, double lbeta, double dInX0,
                    bool electron = false) const {
    if (dInX0 == 0. || p == 0. || lbeta == 0.) {
      return 0.;
    }
    // Highland formula - projected sigma_s
    // ATL-SOFT-PUB-2008-003 equation (15)
    if (!electron) {
      double sigmaSpace = constants::main_RutherfordScott * std::sqrt(dInX0) /
                          (lbeta * p) *
                          (1. + 0.038 * std::log(dInX0 / (lbeta * lbeta)));
      // return the space scattering angle
      return sigmaSpace;
    }
    // Electron specification multiple scattering effects
    // Source: Highland NIM 129 (1975)  p497-499
    // (Highland extension to the Rossi-Greisen formulation)
    // @note: The formula can be extended by replacing
    // the momentum^2 term with pi * pf
    double sigma2 = constants::main_RossiGreisen / (lbeta * p);
    sigma2 *= (sigma2 * dInX0);
    // logarithmic term
    double factor = 1. + constants::log_RossiGreisen * std::log10(10. * dInX0);
    factor *= factor;
    sigma2 *= factor;
    return std::sqrt(sigma2);
  }
};

/// @brief Structure for the energy loss of particles due to radiation in
/// dense material. It combines the effect of bremsstrahlung with direct e+e-
/// pair production and photonuclear interaction. The last two effects are
/// just included for muons.
struct RadiationLoss {
  /// @brief Main call operator for the energy loss. The following equations
  /// are provided by ATL-SOFT-PUB-2008-003.
  ///
  /// @tparam material_t Type of the material
  /// @param [in] E      Energy of the particle
  /// @param [in] m      Mass of the particle
  /// @param [in] mat    Material that is penetrated
  /// @param [in] pdg    PDG code of the particle
  /// @param [in] path   Path length of the particle through the material
  /// @return Radiation energy loss
  template <typename material_t>
  double dEds(double E, double m, const material_t& mat, int pdg,
              double path = 1.) const {
    using namespace Acts::UnitLiterals;

    // Easy exit
    if (mat.X0() == 0.) {
      return 0.;
    }

    // Calculate the bremsstrahlung energy loss (eq. 6)
    const double meOverm = constants::me / m;
    double energyLoss = -E * (meOverm * meOverm);

    // Calculate the energy loss due to direct e+e- pair production and
    // photonuclear interaction (eq. 7, 8) if the particle is a muon
    if ((pdg == 13 || pdg == -13) && E > 8_GeV) {
      if (E < 1_TeV) {
        energyLoss +=
            0.5345 - 6.803e-5 * E - 2.278e-11 * E * E + 9.899e-18 * E * E * E;
      } else {
        energyLoss += 2.986 - 9.253e-5 * E;
      }
    }
    return energyLoss * path / mat.X0();
  }

  /// @brief Evaluation of the energy loss dEds by radiation, direct e+e- pair
  /// production and photonuclear interaction derived by q/p
  /// (=d(dE/ds)/d(q/p)).
  ///
  /// @tparam material_t   Type of the material
  /// @param [in] mass     Mass of the particle
  /// @param [in] material Material that is penetrated
  /// @param [in] qop      Charge over momentum of the particle
  /// @param [in] energy   Energy of the particle
  /// @param [in] pdg      PDG code of the particle
  /// @return The evaluated derivative
  template <typename material_t>
  double dqop(const double mass, const material_t& material, const double qop,
              const double energy, const int pdg) const {
    using namespace Acts::UnitLiterals;

    // Fast exit if material is invalid
    if (material.X0() == 0.) {
      return 0.;
    }

    const double invqop3X0 = 1. / (qop * qop * qop * material.X0());
    double muonExpansion = 0.;
    if ((pdg == 13 || pdg == -13) && energy > 8_GeV) {
      if (energy < 1_TeV) {
        muonExpansion = 6.803e-5 * invqop3X0 / energy +
                        2. * 2.278e-11 * invqop3X0 -
                        3. * 9.899e-18 * invqop3X0 * energy;
      } else {
        muonExpansion = 9.253e-5 * invqop3X0 / energy;
      }
    }
    return constants::me * constants::me * invqop3X0 / (mass * mass * energy) +
           muonExpansion;
  }
};
}  // namespace detail
}  // namespace Acts
