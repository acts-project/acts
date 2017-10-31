// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_EXTRAPOLATIONUTILS_MATERIALINTERACTION_H
#define ACTS_EXTRAPOLATIONUTILS_MATERIALINTERACTION_H 1

#include "ACTS/EventData/ParticleDefinitions.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Utilities/Definitions.hpp"

/// Collection of parametrizations used for
/// energy loss and scattering

namespace Acts {

/// @enum Acts::InteractionType
/// Depending on this type, the energy loss will be calculated differently
enum InteractionType {
  /// Reconstruction - the mean value of the energy loss will be calcuated
  reco = 0,
  /// Simulation - the most probable value of the energy loss will be calculated
  sim = 1,
};

/// The total energy loss consists of the energy loss due to ionization and the
/// energy loss due to radiation

/// The ionization energy loss along a given path length
///
///  Calculate the mean ionization that is pathlength dependent
///
/// The energy loss is calculated using the following paper
/// http://http://pdg.lbl.gov/2014/reviews/rpp2014-rev-passage-particles-matter.pdf
/// Formula 32.5 is used to calculate the mean energy loss [reco mode]
/// Formula 32.11 is used to calculate the mean energy loss [sim mode]
/// Sigma is calculated using the Landau width (FHWM) times the conversion
/// factor 1. / (2. * &radic(2. * log2))
///
/// @param[in] interactionType The type of interaction can either be reco
/// (Reconstruction) - mean energy loss is needed or sim (Simulation) - most
/// probable energy loss needed
/// @param[out] sigma Returns the sigma of the distribution (since it does not
/// exist for the landau distribution the 'Landau Width' is returned)
/// @param[in] p The value of the momentum
/// @param[in] mat The material
/// @param[in] particle The particle type
/// @param[in] particleMasses The masses of the different particles
/// @param[in] path The path length (optional)
/// @return The energy loss due to ionization along a given path length
double
ionizationEnergyLoss(InteractionType       interactionType,
                     double&               sigma,
                     double                p,
                     const Material&       mat,
                     ParticleType          particle,
                     const ParticleMasses& particleMasses,
                     double                path = 1.);

/// @todo To be validated
/// Radiation energy loss along a given path length
///
/// @param[in] p The value of particle momentum
/// @param[in] mat The material
/// @param[in] particleMasses The masses of the different particles
/// @return dEdx from radiation and associate sigma (straggling)
std::pair<double, double>
radiationEnergyLoss(double                p,
                    const Material&       mat,
                    ParticleType          particle,
                    const ParticleMasses& particleMasses);

/// multiple scattering as function of dInX0
///
///
double
sigmaMS(double dInX0, double p, double beta);

}  // end of namespace

#endif
