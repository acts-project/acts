// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Physics/ElectroMagnetic/BetheHeitler.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

ActsFatras::Particle ActsFatras::BetheHeitler::bremPhoton(
    const Particle &particle, Scalar gammaE, Scalar rndPsi, Scalar rndTheta1,
    Scalar rndTheta2, Scalar rndTheta3) const {
  // ------------------------------------------------------
  // simple approach
  // (a) simulate theta uniform within the opening angle of the relativistic
  // Hertz dipole
  //      theta_max = 1/gamma
  // (b)Following the Geant4 approximation from L. Urban -> encapsulate that
  // later
  //      the azimutal angle

  Scalar psi = 2. * M_PI * rndPsi;

  // the start of the equation
  Scalar theta = 0.;
  if (uniformHertzDipoleAngle) {
    // the simplest simulation
    theta = particle.mass() / particle.energy() * rndTheta1;
  } else {
    // ----->
    theta = particle.mass() / particle.energy();
    // follow
    constexpr Scalar a = 0.625;  // 5/8
    Scalar u = -log(rndTheta2 * rndTheta3) / a;
    theta *= (rndTheta1 < 0.25) ? u : u / 3.;  // 9./(9.+27) = 0.25
  }

  Vector3 particleDirection = particle.direction();
  Vector3 photonDirection = particleDirection;

  // construct the combined rotation to the scattered direction
  Acts::RotationMatrix3 rotation(
      // rotation of the scattering deflector axis relative to the reference
      Acts::AngleAxis3(psi, particleDirection) *
      // rotation by the scattering angle around the deflector axis
      Acts::AngleAxis3(theta, Acts::makeCurvilinearUnitU(particleDirection)));
  photonDirection.applyOnTheLeft(rotation);

  Particle photon(particle.particleId().makeDescendant(0),
                  Acts::PdgParticle::eGamma);
  photon.setProcess(ActsFatras::ProcessType::eBremsstrahlung)
      .setPosition4(particle.fourPosition())
      .setDirection(photonDirection)
      .setAbsoluteMomentum(gammaE);
  return photon;
}
