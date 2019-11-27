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
#include "Acts/Utilities/Units.hpp"

namespace Acts {

/// Compute the mean energy loss due to ionisation and excitation.
///
/// @param material  Properties of the traversed material
/// @param thickness Thickness of the traversed material
/// @param pdg       Particle type PDG identifier
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param q         Particle charge
///
/// This computes the expected mean energy loss -dE(x) through a material of
/// the given properties and thickness x, i.e. it computes
///
///     -dE(x) = -dE/dx * x
///
/// where -dE/dx is given by the Bethe formula.
float computeIonisationLossMean(const Material& material, float thickness,
                                int pdg, float m, float qOverP,
                                float q = UnitConstants::e);
/// Derivative of the mean ionisation energy loss with respect to q/p.
///
/// @see computeIonisationLossMean for parameters description
float deriveIonisationLossMeanQOverP(const Material& material, float thickness,
                                     int pdg, float m, float qOverP,
                                     float q = UnitConstants::e);

/// Compute the most propable energy loss due to ionisation and excitation.
///
/// @see computeEnergyLossIonisationMean for parameters description
/// @return Energy loss distribution most probable value and width
///
/// This computes the most probable energy loss -dE(x) through a material of
/// the given properties and thickness as described by the mode of the
/// Landau-Vavilov-Bichsel distribution.
float computeIonisationLossMode(const Material& material, float thickness,
                                int pdg, float m, float qOverP,
                                float q = UnitConstants::e);
/// Derivative of the most probable ionisation energy loss with respect to q/p.
///
/// @see computeIonisationLossMean for parameters description
float deriveIonisationLossModeQOverP(const Material& material, float thickness,
                                     int pdg, float m, float qOverP,
                                     float q = UnitConstants::e);

/// Compute the Gaussian-equivalent sigma for the ionisation loss fluctuations.
///
/// @see computeIonisationLossMean for parameters description
float computeIonisationLossSigma(const Material& material, float thickness,
                                 int pdg, float m, float qOverP,
                                 float q = UnitConstants::e);
/// Compute q/p Gaussian-equivalent sigma due to ionisation loss fluctuations.
///
/// @see computeIonisationLossMean for parameters description
float computeIonisationLossSigmaQOverP(const Material& material,
                                       float thickness, int pdg, float m,
                                       float qOverP,
                                       float q = UnitConstants::e);

/// Compute the mean energy loss due to radiative effects at high energies.
///
/// @param material  Properties of the traversed material
/// @param thickness Thickness of the traversed material
/// @param pdg       Particle type PDG identifier
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param q         Particle charge
///
/// This computes the mean energy loss -dE(x) using an approximative formula.
/// Bremsstrahlung is always included; direct e+e- pair production and
/// photo-nuclear interactions only for muons.
float computeEnergyLossRadiative(const Material& material, float thickness,
                                 int pdg, float m, float qOverP,
                                 float q = UnitConstants::e);
/// Derivative of the mean radiative energy loss with respect to q/p.
///
/// @see computeEnergyLossRadiative for parameters description
float deriveEnergyLossRadiativeQOverP(const Material& material, float thickness,
                                      int pdg, float m, float qOverP,
                                      float q = UnitConstants::e);

/// Compute the combined mean energy loss.
///
/// @param material  Properties of the traversed material
/// @param thickness Thickness of the traversed material
/// @param pdg       Particle type PDG identifier
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param q         Particle charge
///
/// This computes the combined mean energy loss -dE(x) including ionisation and
/// radiative effects.
float computeEnergyLossMean(const Material& material, float thickness, int pdg,
                            float m, float qOverP, float q = UnitConstants::e);
/// Derivative of the combined mean energy loss with respect to q/p.
///
/// @see computeEnergyLossMean for parameters description.
float deriveEnergyLossMeanQOverP(const Material& material, float thickness,
                                 int pdg, float m, float qOverP,
                                 float q = UnitConstants::e);

/// Compute the combined most probably energy loss.
///
/// @see computeEnergyLossMean for parameters description.
float computeEnergyLossMode(const Material& material, float thickness, int pdg,
                            float m, float qOverP, float q = UnitConstants::e);
/// Derivative of the combined most probable energy loss with respect to q/p.
///
/// @see computeEnergyLossMean for parameters description.
float deriveEnergyLossModeQOverP(const Material& material, float thickness,
                                 int pdg, float m, float qOverP,
                                 float q = UnitConstants::e);

/// Compute the core width of the projected planar scattering distribution.
///
/// @param material  Properties of the traversed material
/// @param thickness Thickness of the traversed material
/// @param pdg       Particle type PDG identifier
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param q         Particle charge
float computeMultipleScatteringTheta0(const Material& material, float thickness,
                                      int pdg, float m, float qOverP,
                                      float q = UnitConstants::e);

}  // namespace Acts
