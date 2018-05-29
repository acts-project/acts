// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/EventData/ParticleDefinitions.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/Definitions.hpp"

/// Collection of parametrizations used for
/// energy loss and scattering

namespace Acts {

/// The total energy loss consists of the energy loss due to ionization and the
/// energy loss due to radiation.

/// The mean ionization energy loss along a given path length. The mean energy
/// loss should be used for reconstruction.
/// The energy loss is calculated using the following paper
/// http://http://pdg.lbl.gov/2014/reviews/rpp2014-rev-passage-particles-matter.pdf
/// Formula 32.5 is used to calculate the mean energy loss [reco mode]
/// Sigma is calculated using the Landau width (FHWM) times the conversion
/// factor 1. / (2. * &radic(2. * log2))
///
/// @param[in] p The value of the momentum
/// @param[in] mat The material
/// @param[in] particle The particle type
/// @param[in] particleMasses The masses of the different particles
/// @param[in] path The path length (optional)
/// @return A std::pair. The first entry is the mean energy loss due to
/// ionization along a given path length. The second entry is the sigma of the
/// distribution.
std::pair<double, double>
ionizationEnergyLoss_mean(double                p,
                          const Material&       mat,
                          ParticleType          particle,
                          const ParticleMasses& particleMasses,
                          double                path = 1.);

/// The most probable ionization energy loss along a given path length. The
/// meost probable energy loss should be used for simulation (fatras).
/// The energy loss is calculated using the following paper
/// http://http://pdg.lbl.gov/2014/reviews/rpp2014-rev-passage-particles-matter.pdf
/// Formula 32.11 is used to calculate the most probable energy loss [sim mode]
/// Sigma is calculated using the Landau width (FHWM) times the conversion
/// factor 1. / (2. * &radic(2. * log2))
///
/// @param[in] p The value of the momentum
/// @param[in] mat The material
/// @param[in] particle The particle type
/// @param[in] particleMasses The masses of the different particles
/// @param[in] path The path length (optional)
/// @return A std::pair. The first entry is the most probable energy loss due to
/// ionization along a given path length. The second entry is the sigma of the
/// distribution.
std::pair<double, double>
ionizationEnergyLoss_mop(double                p,
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