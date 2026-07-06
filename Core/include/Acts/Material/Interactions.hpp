// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Material/MaterialSlab.hpp"

namespace Acts {

/// Compute the mean energy loss due to ionisation and excitation.
///
/// @param slab      The traversed material and its properties
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param absQ      Absolute particle charge
///
/// This computes the mean ionisation energy loss @f$-dE(x)@f$ of a particle
/// traversing the material slab,
///
/// @f[
///   -dE(x) = -\frac{dE}{dx}\, x,
/// @f]
///
/// where @f$-dE/dx@f$ is given by the Bethe formula
/// @cite ParticleDataGroup:2018ovx (eq. 33.5), including the density-effect
/// correction. The result is the magnitude of the loss (always @f$\geq 0@f$);
/// The formula is valid for intermediate energies,
/// roughly @f$0.1 \lesssim \beta\gamma \lesssim 1000@f$.
///
/// @return Mean ionisation energy loss through the slab in native energy units
/// @see @ref computeEnergyLossLandau for the most probable value,
///   @ref computeEnergyLossRadiative for radiative losses, and
///   @ref computeEnergyLossMean for the sum of both
float computeEnergyLossBethe(const MaterialSlab& slab, float m, float qOverP,
                             float absQ);
/// Derivative of the Bethe energy loss with respect to q/p.
///
/// @copydoc computeEnergyLossBethe
/// @return Derivative of the mean ionisation energy loss with respect to q/p
float deriveEnergyLossBetheQOverP(const MaterialSlab& slab, float m,
                                  float qOverP, float absQ);

/// Compute the most probable energy loss due to ionisation and excitation.
///
/// @param slab      The traversed material and its properties
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param absQ      Absolute particle charge
///
/// This computes the most probable ionisation energy loss @f$-dE(x)@f$ through
/// the material slab, i.e. the mode of the Landau-Vavilov-Bichsel distribution
/// @cite ParticleDataGroup:2018ovx (eq. 33.12), including the density-effect
/// correction. Unlike @ref computeEnergyLossBethe (which returns the mean), this is
/// the most probable value, which for thin slabs is noticeably smaller than the
/// mean because of the long tail of the distribution. The formula is valid for
/// intermediate energies, roughly @f$0.1 \lesssim \beta\gamma \lesssim 1000@f$.
///
/// @return Most probable ionisation energy loss through the slab in native
///   energy units
float computeEnergyLossLandau(const MaterialSlab& slab, float m, float qOverP,
                              float absQ);
/// Derivative of the most probable ionisation energy loss with respect to q/p.
///
/// @copydoc computeEnergyLossLandau
/// @return Derivative of the most probable ionisation energy loss with respect
///   to q/p
float deriveEnergyLossLandauQOverP(const MaterialSlab& slab, float m,
                                   float qOverP, float absQ);

/// Compute the Gaussian-equivalent sigma for the ionisation loss fluctuations.
///
/// @see @ref computeEnergyLossBethe for parameters description
///
/// This is the sigma parameter of a Gaussian distribution with the same
/// full-width-half-maximum as the Landau-Vavilov-Bichsel distribution. The
/// computations are valid for intermediate particle energies.
/// @param slab The traversed material and its properties
/// @param m Particle mass
/// @param qOverP Particle charge divided by absolute momentum
/// @param absQ Absolute particle charge
/// @return Gaussian-equivalent sigma for energy loss fluctuations
float computeEnergyLossLandauSigma(const MaterialSlab& slab, float m,
                                   float qOverP, float absQ);

/// Compute the full with half maximum of landau energy loss distribution
///
/// @param slab The traversed material and its properties
/// @param m Particle mass
/// @param qOverP Particle charge divided by absolute momentum
/// @param absQ Absolute particle charge
/// @return Full width half maximum of the Landau distribution
float computeEnergyLossLandauFwhm(const MaterialSlab& slab, float m,
                                  float qOverP, float absQ);

/// Compute the Gaussian-equivalent sigma of q/p due to ionisation fluctuations.
///
/// @param slab   The traversed material and its properties
/// @param m      Particle mass
/// @param qOverP Particle charge divided by absolute momentum
/// @param absQ   Absolute particle charge
///
/// This propagates the energy-loss straggling (the Gaussian-equivalent sigma
/// @f$\sigma_E@f$ from @ref computeEnergyLossLandauSigma) into a standard deviation
/// on @f$q/p@f$ using the Jacobian @f$d(q/p)/dE@f$,
///
/// @f[
///   \sigma_{q/p} = \left|\frac{d(q/p)}{dE}\right|\, \sigma_E .
/// @f]
///
/// This is the quantity used as the @f$q/p@f$ process noise in the Kalman
/// fitters.
///
/// @return Gaussian-equivalent standard deviation of q/p
float computeEnergyLossLandauSigmaQOverP(const MaterialSlab& slab, float m,
                                         float qOverP, float absQ);

/// Compute the mean energy loss due to radiative effects at high energies.
///
/// @param slab      The traversed material and its properties
/// @param absPdg    Absolute particle type PDG identifier
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param absQ      Absolute particle charge
///
/// This computes the mean radiative energy loss @f$-dE(x)@f$ using the
/// approximation of @cite Lund:2008ad. Bremsstrahlung, scaling with
/// @f$(m_e/m)^2@f$, is always included; direct @f$e^+e^-@f$ pair production and
/// photo-nuclear interactions are added only for muons above @f$8\,\mathrm{GeV}@f$.
/// Like @ref computeEnergyLossBethe the result is the magnitude of the loss (always
/// @f$\geq 0@f$).
///
/// @return Mean radiative energy loss through the slab in native energy units
float computeEnergyLossRadiative(const MaterialSlab& slab, PdgParticle absPdg,
                                 float m, float qOverP, float absQ);
/// Derivative of the mean radiative energy loss with respect to q/p.
///
/// @copydoc computeEnergyLossRadiative
/// @return Derivative of radiative energy loss with respect to q/p
float deriveEnergyLossRadiativeQOverP(const MaterialSlab& slab,
                                      PdgParticle absPdg, float m, float qOverP,
                                      float absQ);

/// Compute the combined mean energy loss.
///
/// @param slab      The traversed material and its properties
/// @param absPdg    Absolute particle type PDG identifier
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param absQ      Absolute particle charge
///
/// This computes the combined mean energy loss -dE(x) including ionisation and
/// radiative effects. The computations are valid over a wide range of particle
/// energies.
/// @return Combined mean energy loss through the material slab
float computeEnergyLossMean(const MaterialSlab& slab, PdgParticle absPdg,
                            float m, float qOverP, float absQ);
/// Derivative of the combined mean energy loss with respect to q/p.
///
/// @copydoc computeEnergyLossMean
/// @return Derivative of combined mean energy loss with respect to q/p
float deriveEnergyLossMeanQOverP(const MaterialSlab& slab, PdgParticle absPdg,
                                 float m, float qOverP, float absQ);

/// Compute the combined most probably energy loss.
///
/// @copydoc computeEnergyLossMean
/// @return Combined most probable energy loss through the material slab
float computeEnergyLossMode(const MaterialSlab& slab, PdgParticle absPdg,
                            float m, float qOverP, float absQ);
/// Derivative of the combined most probable energy loss with respect to q/p.
///
/// @copydoc computeEnergyLossMean
/// @return Derivative of combined most probable energy loss with respect to q/p
float deriveEnergyLossModeQOverP(const MaterialSlab& slab, PdgParticle absPdg,
                                 float m, float qOverP, float absQ);

/// Compute the core width of the projected planar scattering distribution.
///
/// @param slab      The traversed material and its properties
/// @param absPdg    Absolute particle type PDG identifier
/// @param m         Particle mass
/// @param qOverP    Particle charge divided by absolute momentum
/// @param absQ      Absolute particle charge
///
/// The returned @f$\theta_0@f$ is the standard deviation of the central
/// (Gaussian) part of the multiple-Coulomb-scattering angle, projected onto a
/// plane. For all particles except electrons and positrons it is evaluated with
/// the Highland formula @cite Highland:1975pq in the parametrisation of
/// @cite ParticleDataGroup:2018ovx (eq. 33.15); for electrons and positrons the
/// Rossi-Greisen form is used instead.
///
/// @note This is the projected (single-plane) width; the width of the polar
///   space angle is larger by a factor @f$\sqrt{2}@f$.
/// @return Projected scattering angle standard deviation @f$\theta_0@f$ in radians
/// @see @ref approximateHighlandScattering for a charge- and momentum-independent
///   approximation
float computeMultipleScatteringTheta0(const MaterialSlab& slab,
                                      PdgParticle absPdg, float m, float qOverP,
                                      float absQ);

/// Approximate the core width of the projected planar scattering distribution
/// with Highland's formula.
///
/// In contrast to @ref computeMultipleScatteringTheta0, this ignores the particle
/// charge and velocity (assuming a singly-charged, ultra-relativistic particle
/// with @f$q^2/\beta^2 = 1@f$) and does not divide by the momentum. It therefore
/// returns @f$\theta_0 \cdot p@f$ rather than @f$\theta_0@f$ itself, which is
/// convenient when the momentum is not yet known.
///
/// @param xOverX0  The thickness of the material in radiation lengths
/// @return The projected scattering angle scaled by momentum,
///   @f$\theta_0 \cdot p@f$, in native units (radians times momentum)
float approximateHighlandScattering(float xOverX0);

}  // namespace Acts
