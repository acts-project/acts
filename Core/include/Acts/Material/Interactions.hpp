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
/// @see computeIonisationLossMean for parameters description
float deriveIonisationLossModeQOverP(const Material& material, float thickness,
                                     float m, float qOverP,
                                     float q = UnitConstants::e);

/// Compute the expected energy loss due to radiation effects.
///
/// @param material  Properties of the traversed material
/// @param thickness Thickness of the traversed material
/// @param pdg       Particle type PDG identifier
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param q         Particle charge
///
/// Bremsstrahlung is always included. Direct e+e- pair production and
/// photo-nuclear interactions only for muons.
float computeRadiationLoss(const Material& material, float thickness, int pdg,
                           float m, float qOverP, float q = UnitConstants::e);
/// Derivative of the expected radiation energy loss with respect to q/p.
///
/// @see computeRadiationLoss for parameters description
float deriveRadiationLossQOverP(const Material& material, float thickness,
                                int pdg, float m, float qOverP,
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

}  // namespace detail
}  // namespace Acts
