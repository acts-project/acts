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
    /// ionization along a given path length or its mode. The second entry is
    /// the sigma of the distribution.
    template <typename material_t>
    std::pair<double, double>
    dEds(double            m,
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
      } else {
        kaz = 0.5 * constants::ka_BetheBloch * mat.zOverAtimesRho();
      }
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
      double sigma = kazL * constants::landau2gauss;
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

    /// @brief Evaluation of the energy loss dEds by ionisation derived by q/p
    /// (=d(dE/ds)/d(q/p)). The underlying equations for dE/ds are given by
    /// Formula 33.5 (mean energy loss) and Formula 33.11 (most probable energy
    /// loss).
    ///
    /// @tparam material_t Type of the material
    /// @param [in] energy Energy of the particle
    /// @param [in] qop Charge over momentum of the particle
    /// @param [in] mass Mass of the particle
    /// @param [in] material Material that is penetrated
    /// @param [in] mean Boolean flag if mean or mode of dEds should be derived
    /// @param [in] siUnits Boolean flag if SI or natural units should be used
    ///
    /// @note The strategy is to evaluate everything in natural units as long as
    /// they are involved in terms that lead to canceling out of units. So, only
    /// the terms that define the resulting unit are treated in the individual
    /// unit system. Therewith the amount of calculations of constants and
    /// conversions gets reduced.
    /// @return The evaluated derivative
    template <typename material_t>
    double
    dqop(const double      energy,
         const double      qop,
         const double      mass,
         const material_t& material,
         const bool        mean    = true,
         const bool        siUnits = false) const
    {
      // Fast exit if material is invalid
      if (material.Z() == 0 || material.zOverAtimesRho() == 0) {
        return 0.;
      }

      // Constants for readability
      const double qop1
          = siUnits ? qop / units::SI2Nat<units::MOMENTUM>(1.) : qop;
      const double m  = siUnits ? units::SI2Nat<units::MASS>(mass) : mass;
      const double m2 = m * m;
      const double I  = constants::eionisation * std::pow(material.Z(), 0.9);
      const double E  = siUnits ? units::SI2Nat<units::ENERGY>(energy) : energy;
      const double gamma = E / m;
      const double beta  = std::abs(1 / (E * qop1));
      const double beta2 = beta * beta;
      const double kaz   = siUnits
          ? 0.5 * units::Nat2SI<units::ENERGY>(30.7075 * units::_MeV)
              * units::_mm2 * units::_e2 / units::_g * material.zOverAtimesRho()
          : 0.5 * constants::ka_BetheBloch * material.zOverAtimesRho();
      // Storage of the result
      double betheBlochDerivative = 0.;
      if (mean) {
        // Constants related to the mean evaluation
        const double me   = constants::me;
        const double me2  = me * me;
        const double qop3 = qop1 * qop1 * qop1;
        const double qop4 = qop3 * qop1;
        const double m4   = m2 * m2;
        const double I2   = I * I;

        // Parts of the derivative
        const double lnCore = 4. * me2 / (m4 * I2 * qop4)
            / (1. + 2. * gamma * me / m + me2 / m2);
        const double lnCoreDerivative = -4. * me2 / (m4 * I2)
            * std::pow(qop4 + 2. * gamma * qop4 * me / m + qop4 * me2 / m2, -2.)
            * (4. * qop3 + 8. * me * qop3 * gamma / m
               - 2. * me * qop1 / (m2 * m * gamma)
               + 4. * qop3 * me2 / m2);

        // Combine parts
        const double lnDerivative = 2. * qop1 * m2 * std::log(lnCore)
            + lnCoreDerivative / (lnCore * beta2);
        betheBlochDerivative = -kaz * lnDerivative;
      } else {
        // Parts of the mode derivative
        const double lnCore   = std::log(2. * m * beta2 * gamma * gamma / I);
        const double lnCoreXi = std::log(kaz / (beta2 * I));
        const double lnCoreDerivative   = -2. * m * beta * gamma;
        const double lnCoreXiDerivative = 2. * beta2 * m2;

        // Result evaluation (without the density effect)
        betheBlochDerivative
            = kaz * (2. * m2 * qop1 * (lnCore + lnCoreXi + 0.2)
                     + (lnCoreDerivative + lnCoreXiDerivative) / beta2);
      }
      // Density effect, only valid for high energies (gamma > 10 -> p > 1GeV
      // for muons). The derivative includes the kinematic parts of the
      // pre-factor
      if (gamma > 10.) {
        const double delta
            = 2. * std::log(constants::eplasma
                            * std::sqrt(1000. * material.zOverAtimesRho())
                            / I)
            + 2. * std::log(beta * gamma) - 1.;
        if (mean) {
          const double deltaDerivative
              = -2. / (qop1 * beta2) + 2. * delta * qop1 * m2;
          betheBlochDerivative += kaz * deltaDerivative;
        } else {
          const double deltaDerivative
              = -kaz * 2. / (qop1 * beta2) + kaz * 2. * m2 * qop1 * delta;
          betheBlochDerivative -= deltaDerivative;
        }
      }
      // Convert to right unit system and return
      if (siUnits) {
        return betheBlochDerivative * units::Nat2SI<units::MOMENTUM>(1.);
      } else {
        return betheBlochDerivative;
      }
    }
  };

  /// @brief Multiple scattering as function of dInX0
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

  /// @brief Structure for the energy loss of particles due to radiation in
  /// dense material. It combines the effect of bremsstrahlung with direct e+e-
  /// pair production and photonuclear interaction. The last two effects are
  /// just included for muons.
  struct RadiationLoss
  {
    /// @brief Main call operator for the energy loss. The following equations
    /// are provided by ATL-SOFT-PUB-2008-003.
    ///
    /// @tparam material_t Type of the material
    /// @param [in] E Energy of the particle
    /// @param [in] m Mass of the particle
    /// @param [in] mat Material that is penetrated
    /// @param [in] pdg PDG code of the particle
    /// @param [in] path Path length of the particle through the material
    /// @param [in] siUnits Boolean flag if SI or natural units should be used
    /// @return Radiation energy loss
    template <typename material_t>
    double
    dEds(double            E,
         double            m,
         const material_t& mat,
         int               pdg,
         double            path    = 1.,
         bool              siUnits = false) const
    {
      // Easy exit
      if (mat.X0() == 0.) {
        return 0.;
      }

      double       energyLoss;
      const double meOverm
          = constants::me / (siUnits ? units::SI2Nat<units::MASS>(m) : m);

      // Converting energy if needed
      if (siUnits) {
        E = units::SI2Nat<units::ENERGY>(E);
      }

      // Calculate the bremsstrahlung energy loss (eq. 6)
      energyLoss = -E * (meOverm * meOverm);

      // Calculate the energy loss due to direct e+e- pair production and
      // photonuclear interaction (eq. 7, 8) if the particle is a muon
      if ((pdg == 13 || pdg == -13) && E > 8. * units::_GeV) {
        if (E < 1. * units::_TeV) {
          energyLoss += 0.5345 - 6.803e-5 * E - 2.278e-11 * E * E
              + 9.899e-18 * E * E * E;
        } else {
          energyLoss += 2.986 - 9.253e-5 * E;
        }
      }

      // Return energy loss
      if (siUnits) {
        return units::Nat2SI<units::ENERGY>(energyLoss) * path / mat.X0();
      }
      return energyLoss * path / mat.X0();
    }

    /// @brief Evaluation of the energy loss dEds by radiation, direct e+e- pair
    /// production and photonuclear interaction derived by q/p
    /// (=d(dE/ds)/d(q/p)).
    ///
    /// @tparam material_t Type of the material
    /// @param [in] mass Mass of the particle
    /// @param [in] material Material that is penetrated
    /// @param [in] qop Charge over momentum of the particle
    /// @param [in] energy Energy of the particle
    /// @param [in] pdg PDG code of the particle
    /// @param [in] siUnits Boolean flag if SI or natural units should be used
    /// @return The evaluated derivative
    template <typename material_t>
    double
    dqop(const double      mass,
         const material_t& material,
         const double      qop,
         const double      energy,
         const int         pdg,
         const bool        siUnits = false) const
    {
      // Fast exit if material is invalid
      if (material.X0() == 0.) {
        return 0.;
      }

      const double invqop3X0 = 1. / (qop * qop * qop * material.X0());

      double muonExpansion = 0.;
      if ((pdg == 13 || pdg == -13)
          && (siUnits ? units::SI2Nat<units::ENERGY>(energy) : energy)
              > 8. * units::_GeV) {
        if ((siUnits ? units::SI2Nat<units::ENERGY>(energy) : energy)
            < 1. * units::_TeV) {
          muonExpansion = 6.803e-5 * invqop3X0 / energy
              + 2. * 2.278e-11 * invqop3X0
              - 3. * 9.899e-18 * invqop3X0 * energy;
        } else {
          muonExpansion = 9.253e-5 * invqop3X0 / energy;
        }
      }
      if (siUnits) {
        // Just rescale the mass to natural units, qop & energy are already
        // given in the right system
        const double scaling = units::SI2Nat<units::MASS>(1.);
        return constants::me * constants::me * invqop3X0
            / (mass * mass * energy * scaling * scaling)
            + muonExpansion;
      } else {
        return constants::me * constants::me * invqop3X0
            / (mass * mass * energy)
            + muonExpansion;
      }
    }
  };
}  // namespace detail
}  // namespace Acts
