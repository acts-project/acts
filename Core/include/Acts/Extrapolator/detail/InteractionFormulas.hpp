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
    /// @param[in] p The value of the momentum
    /// @param[in] m The masses of the different particles
    /// @param[in] mat The material that is traversed
    /// @param[in] path The path length (optional)
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
               double            path = 1.,
               bool              mean = true) const
    {

      // the return value
      double dE = 0.;

      // Ionization - Bethe-Bloch
      // See ATL-SOFT-PUB-2008-003 equation (4)
      // 16 eV * Z**0.9
      double I = constants::eionisation * std::pow(mat.Z(), 0.9);

      // See (1) table 33.1
      // K/A*Z = 0.5 * 30.7075MeV/(g/mm2) * Z/A * rho[g/mm3]
      double kaz  = 0.5 * constants::ka_BetheBloch * mat.zOverAtimesRho();
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

      // The landau width (FWHM) is 4.*kazL
      // The factor is the conversion factor from FWHM to sigma for
      // gaussian curve: 1. / (2. * sqrt(2. * log(2.))).
      double sigma = 2. * kazL * 1. / (std::sqrt(2. * std::log(2.)));
      if (mean) {
        // calculate the fraction to the electron mass
        double mfrac = constants::me / m;
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
        dE = kazL * (std::log(2. * m * eta2 / I) + std::log(kazL / I) + 0.2
                     - (lbeta * lbeta)
                     - delta);
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
      if (dInX0 == 0. || p == 0. || lbeta == 0.) return 0.;
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

}  // namespace detail
}  // namespace Acts
