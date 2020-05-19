// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <utility>

#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

/// Compute the mean energy loss due to ionisation and excitation.
///
/// @param slab      The traversed material and its properties
/// @param pdg       Particle type PDG identifier
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param q         Particle charge
///
/// This computes the mean energy loss -dE(x) through a material with
/// the given properties, i.e. it computes
///
///     -dE(x) = -dE/dx * x
///
/// where -dE/dx is given by the Bethe formula. The computations are valid
/// for intermediate particle energies.
float computeEnergyLossBethe(const MaterialProperties& slab, int pdg, float m,
                             float qOverP, float q = UnitConstants::e);
/// Derivative of the Bethe energy loss with respect to q/p.
///
/// @see computeEnergyLossBethe for parameters description
float deriveEnergyLossBetheQOverP(const MaterialProperties& slab, int pdg,
                                  float m, float qOverP,
                                  float q = UnitConstants::e);

/// Compute the most propable energy loss due to ionisation and excitation.
///
/// @see computeEnergyLossBethe for parameters description
///
/// This computes the most probable energy loss -dE(x) through a material of
/// the given properties and thickness as described by the mode of the
/// Landau-Vavilov-Bichsel distribution. The computations are valid
/// for intermediate particle energies.
float computeEnergyLossLandau(const MaterialProperties& slab, int pdg, float m,
                              float qOverP, float q = UnitConstants::e);
/// Derivative of the most probable ionisation energy loss with respect to q/p.
///
/// @see computeEnergyLossBethe for parameters description
float deriveEnergyLossLandauQOverP(const MaterialProperties& slab, int pdg,
                                   float m, float qOverP,
                                   float q = UnitConstants::e);

/// Compute the Gaussian-equivalent sigma for the ionisation loss fluctuations.
///
/// @see computeEnergyLossBethe for parameters description
///
/// This is the sigma paramter of a Gaussian distribution with the same
/// full-width-half-maximum as the Landau-Vavilov-Bichsel distribution. The
/// computations are valid for intermediate particle energies.
float computeEnergyLossLandauSigma(const MaterialProperties& slab, int pdg,
                                   float m, float qOverP,
                                   float q = UnitConstants::e);
/// Compute q/p Gaussian-equivalent sigma due to ionisation loss fluctuations.
///
/// @see computeEnergyLossBethe for parameters description
float computeEnergyLossLandauSigmaQOverP(const MaterialProperties& slab,
                                         int pdg, float m, float qOverP,
                                         float q = UnitConstants::e);

/// Compute the mean energy loss due to radiative effects at high energies.
///
/// @param slab      The traversed material and its properties
/// @param pdg       Particle type PDG identifier
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param q         Particle charge
///
/// This computes the mean energy loss -dE(x) using an approximative formula.
/// Bremsstrahlung is always included; direct e+e- pair production and
/// photo-nuclear interactions only for muons.
float computeEnergyLossRadiative(const MaterialProperties& slab, int pdg,
                                 float m, float qOverP,
                                 float q = UnitConstants::e);
/// Derivative of the mean radiative energy loss with respect to q/p.
///
/// @see computeEnergyLossRadiative for parameters description
float deriveEnergyLossRadiativeQOverP(const MaterialProperties& slab, int pdg,
                                      float m, float qOverP,
                                      float q = UnitConstants::e);

/// Compute the combined mean energy loss.
///
/// @param slab      The traversed material and its properties
/// @param pdg       Particle type PDG identifier
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param q         Particle charge
///
/// This computes the combined mean energy loss -dE(x) including ionisation and
/// radiative effects. The computations are valid over a wide range of particle
/// energies.
float computeEnergyLossMean(const MaterialProperties& slab, int pdg, float m,
                            float qOverP, float q = UnitConstants::e);
/// Derivative of the combined mean energy loss with respect to q/p.
///
/// @see computeEnergyLossMean for parameters description.
float deriveEnergyLossMeanQOverP(const MaterialProperties& slab, int pdg,
                                 float m, float qOverP,
                                 float q = UnitConstants::e);

/// Compute the combined most probably energy loss.
///
/// @see computeEnergyLossMean for parameters description.
float computeEnergyLossMode(const MaterialProperties& slab, int pdg, float m,
                            float qOverP, float q = UnitConstants::e);
/// Derivative of the combined most probable energy loss with respect to q/p.
///
/// @see computeEnergyLossMean for parameters description.
float deriveEnergyLossModeQOverP(const MaterialProperties& slab, int pdg,
                                 float m, float qOverP,
                                 float q = UnitConstants::e);

/// Compute the core width of the projected planar scattering distribution.
///
/// @param slab      The traversed material and its properties
/// @param pdg       Particle type PDG identifier
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param q         Particle charge
float computeMultipleScatteringTheta0(const MaterialProperties& slab, int pdg,
                                      float m, float qOverP,
                                      float q = UnitConstants::e);

}  // namespace Acts
