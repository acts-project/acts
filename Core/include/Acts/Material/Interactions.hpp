// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <utility>

#include "Acts/Material/Material.hpp"
#include "Acts/Material/detail/Constants.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

/// Compute the mean energy loss due to ionisation.
///
/// @param material  Properties of the traversed material
/// @param thickness Thickness of the traversed material
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param q         Particle charge
/// @return Energy loss distribution mean and width
///
/// This computes the expected mean energy loss -dE(x) through a material of
/// the given properties and thickness x, i.e. it computes
///
///     -dE(x) = -dE/dx * x
///
std::pair<float, float> computeIonisationLossMean(const Material& material,
                                                  float thickness, float m,
                                                  float qOverP,
                                                  float q = UnitConstants::e);

/// Derivative of the mean ionisation energy loss with respect to q/p.
///
/// @see computeIonisationLossMean for parameters description
float deriveIonisationLossMeanQOverP(const Material& material, float thickness,
                                     float m, float qOverP,
                                     float q = UnitConstants::e);

/// Compute the most propable energy loss due to ionisation.
///
/// @see computeEnergyLossIonisationMean for parameters description
/// @return Energy loss distribution most probable value and width
///
/// This computes the most probable energy loss dE through a material of the
/// given properties and thickness as described by the mode of the
/// Landau-Vavilov-Bichsel distribution.
std::pair<float, float> computeIonisationLossMode(const Material& material,
                                                  float thickness, float m,
                                                  float qOverP,
                                                  float q = UnitConstants::e);

/// Derivative of the most probable ionisation energy loss with respect to q/p.
///
/// @see computeIonisationLossMode for parameters description
float deriveIonisationLossModeQOverP(const Material& material, float thickness,
                                     float m, float qOverP,
                                     float q = UnitConstants::e);

namespace detail {

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
