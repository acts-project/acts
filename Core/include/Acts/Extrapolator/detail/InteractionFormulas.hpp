// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Constants.hpp"

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
  struct IonisationLoss
  {

    /// @brief call operator for the Ionisation class
    ///
    /// @tparam material_t Type of the material class
    ///
    /// @param [in] m The masses of the different particles
    /// @param [in] lbeta Beta factor of the particle
    /// @param [in] lgamma Gamma factor of the particle
    /// @param [in] mat The material that is traversed
    /// @param [in] path The path length (optional)
    /// @param [in] mean Toggle between calculation of mean (true) and mode
    /// (false)
    /// @param [in] siUnits Toggle for the unit system between SI and natural
    /// units
    /// @note The only conversion needed is for @p kazL if siUnits is true.
    ///
    /// @return A std::pair. The first entry is the mean energy loss due to
    /// ionization along a given path length. The second entry is the sigma of
    /// the
    /// distribution.
    template <typename material_t>
    std::pair<double, double>
    operator()(double            m,
               double            lbeta,
               double            lgamma,
               const material_t& mat,
               double            path    = 1.,
               bool              mean    = true,
               bool              siUnits = false) const
    {

      // the return value
      double dE = 0.;

      // Ionization - Bethe-Bloch
      // See ATL-SOFT-PUB-2008-003 equation (4)
      // 16 eV * Z**0.9
      double I = constants::eionisation * std::pow(mat.Z(), 0.9);

      // See (1) table 33.1
      // K/A*Z = 0.5 * 30.7075MeV/(g/mm2) * Z/A * rho[g/mm3]
      // Convert energy if SI units are used
      double kaz;
      if (siUnits) {
        kaz = 0.5 * units::Nat2SI<units::ENERGY>(30.7075 * units::_MeV)
            * units::_mm * units::_mm / units::_g * mat.zOverAtimesRho();
      } else
        kaz = 0.5 * constants::ka_BetheBloch * mat.zOverAtimesRho();

      double eta2 = lbeta * lgamma;
      eta2 *= eta2;
      // density effect, only valid for high energies
      // (lgamma > 10 -> p > 1 GeV for muons)
      double delta = 0.;
      if (lbeta * lgamma > 10.) {
        // See (1) table 33.1
        double eplasma
            = constants::eplasma * sqrt(1000. * mat.zOverAtimesRho());
        // See (1) formula 33.6
        delta = 2. * std::log(eplasma / I) + std::log(eta2) - 1.;
      }

      // divide by lbeta^2 for non-electrons
      kaz /= lbeta * lbeta;
      double kazL = kaz * path;

      // Multiply the charge for consistency
      if (siUnits) {
        kazL *= units::_e * units::_e;
      }

      // The landau width (FWHM) is 4.*kazL
      // The factor is the conversion factor from FWHM to sigma for
      // gaussian curve: 1. / (2. * sqrt(2. * log(2.))).
      double sigma = 2. * kazL * 1. / (std::sqrt(2. * std::log(2.)));
      if (mean) {
        // Calculate the fraction to the electron mass
        double mfrac = siUnits ? (constants::me / units::SI2Nat<units::MASS>(m))
                               : (constants::me / m);
        // Calculate the mean value for reconstruction
        // See ATL-SOFT-PUB-2008-003 equation (2)
        double tMax = 2. * eta2 * constants::me
            / (1. + 2. * lgamma * mfrac + mfrac * mfrac);
        // See ATL-SOFT-PUB-2008-003 equation (1)
        // or:
        // http://pdg.lbl.gov/2018/reviews/rpp2018-rev-passage-particles-matter.pdf
        // PDG formula 33.5
        dE = -kazL * 2.0
            * (0.5 * std::log(2. * constants::me * eta2 * tMax / (I * I))
               - (lbeta * lbeta)
               - delta * 0.5);
      } else {
        // Calculate the most probably value for simulation
        //
        // the landau sigmaL is path length dependent
        //    PDG formula 33.11 for MOP value from
        //    http://pdg.lbl.gov/2018/reviews/rpp2018-rev-passage-particles-matter.pdf
        //
        if (siUnits) {
          dE = kazL * (std::log(2. * units::SI2Nat<units::MASS>(m) * eta2 / I)
                       + std::log(units::SI2Nat<units::ENERGY>(kazL) / I)
                       + 0.2
                       - (lbeta * lbeta)
                       - delta);
        } else {
          dE = kazL * (std::log(2. * m * eta2 / I) + std::log(kazL / I) + 0.2
                       - (lbeta * lbeta)
                       - delta);
        }
      }
      // return the energy loss and stragling
      return std::make_pair(dE, sigma);
    }
  };

  /// @brield Multiple scattering as function of dInX0
  ///
  /// It supports MIP and electron scattering and return
  /// the space angle which has to be transformed into actual
  /// scattering components for the azimuthal and polar angles,
  /// respectively (for reconstruction).
  struct HighlandScattering
  {

    /// @brief call operator for the HighlandScattering formula
    ///
    /// @param p is the momentum
    /// @param lbeta is the Lorentz beta parameter
    /// @param dInX0 is the passed thickness expressed in units of X0
    double
    operator()(double p,
               double lbeta,
               double dInX0,
               bool   electron = false) const
    {
      if (dInX0 == 0. || p == 0. || lbeta == 0.) {
        return 0.;
      }
      // Highland formula - projected sigma_s
      // ATL-SOFT-PUB-2008-003 equation (15)
      if (!electron) {
        double sigmaSpace = constants::main_RutherfordScott * std::sqrt(dInX0)
            / (lbeta * p) * (1. + 0.038 * std::log(dInX0 / (lbeta * lbeta)));
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
      double factor
          = 1. + constants::log_RossiGreisen * std::log10(10. * dInX0);
      factor *= factor;
      sigma2 *= factor;
      return std::sqrt(sigma2);
    }
  };

  /// @brief Structure for the energy loss of muons in dense material. It
  /// combines the effect of energy loss by ionisation with bremsstrahlung,
  /// direct e+e- pair production and photonuclear interaction.
  struct MuonEnergyLoss
  {
    /// @brief Main call operator for the energy loss. The following equations
    /// are provided by ATL-SOFT-PUB-2008-003.
    ///
    /// @tparam material_t Type of the material
    /// @param [in] m Masst of the particle
    /// @param [in] lbeta Beta factor of the particle
    /// @param [in] lgamma Gamma factor of the particle
    /// @param [in] mat Material that is penetrated
    /// @param [in] path Path length of the particle through the material
    /// @param [in] mean Boolean flag for using the mean or the mode for
    /// ionisation losses
    template <typename material_t>
    std::pair<double, double>
    operator()(double            m,
               double            lbeta,
               double            lgamma,
               const material_t& mat,
               double            path = 1.,
               bool              mean = true)
    {
      // Calculate the ionisation energy loss
      IonisationLoss iLoss;
      auto           energyLoss = iLoss(m, lbeta, lgamma, mat, path, mean);

      // Calculate the bremsstrahlung energy loss (eq. 6)
      const double p       = m * lgamma * lbeta;
      const double E       = std::sqrt(p * p + m * m);
      const double meOverm = constants::me / m;
      energyLoss.first += E / mat.X0 * (meOverm * meOverm) * path;

      // Calculate the energy loss due to direct e+e- pair production and
      // photonuclear interaction (eq. 7, 8)
      if (E > 8. * units::_GeV) {
        const double invX0 = 1. / mat.X0;
        if (E > 1. * units::_TeV) {
          energyLoss.first += (0.5345 - 6.803e-5 * E - 2.278e-11 * E * E
                               + 9.899e-18 * E * E * E)
              * invX0 * path;
        } else {
          energyLoss.first += (2.986 - 9.253e-5 * E) * invX0 * path;
        }
      }

      // Calculate the width of the energy loss (eq. 12 - 14)
      const double a      = 121 + 3.9e-3 * E + 5.3e-9 * E * E;
      const double Eloss2 = a * path * mat.rho() * constants::ka_BetheBloch
          * 0.5 / (lbeta * lbeta);
      energyLoss.second = std::sqrt(Eloss2) / (lbeta * p * p);
      return energyLoss;
    }
  };

}  // namespace detail
}  // namespace Acts
