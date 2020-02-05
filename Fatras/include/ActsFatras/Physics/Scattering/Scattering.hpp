// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <random>
#include <vector>

#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsFatras/EventData/Particle.hpp"

namespace ActsFatras {

/// Simulate (multiple) scattering using a configurable scattering model.
///
/// @tparam scattering_model_t Model implementation to draw a scattering angle.
template <typename scattering_model_t>
struct Scattering {
  /// The flag to include scattering or not
  bool scattering = true;
  /// Include the log term
  bool parametric = false;
  /// The scattering formula
  scattering_model_t angle;

  /// Simulate scattering and update the particle parameters.
  ///
  /// @param[in]     generator is the random number generator
  /// @param[in]     slab      defines the passed material
  /// @param[in,out] particle  is the particle being updated
  /// @return Empty secondaries containers.
  ///
  /// @tparam generator_t is a RandomNumberEngine
  template <typename generator_t>
  std::vector<Particle> operator()(generator_t &generator,
                                   const Acts::MaterialProperties &slab,
                                   Particle &particle) const {
    // Do nothing if the flag is set to false
    if (not scattering) {
      return {};
    }

    // 3D scattering angle
    double angle3D = angle(generator, slab, particle);

    // parametric scattering
    if (parametric) {
      // the initial values
      double theta = Acts::VectorHelpers::theta(particle.direction());
      double phi = Acts::VectorHelpers::phi(particle.direction());
      double sinTheta = (sin(theta) * sin(theta) > 10e-10) ? sin(theta) : 1.;

      // sample them in an independent way
      const double projectionFactor = 1. / std::sqrt(2.);
      double deltaTheta = projectionFactor * angle3D;
      double numDetlaPhi = 0.;  //?? @THIS IS WRONG HERE !
      double deltaPhi = projectionFactor * numDetlaPhi / sinTheta;

      // @todo: use bound parameter
      // (i) phi
      phi += deltaPhi;
      if (phi >= M_PI)
        phi -= M_PI;
      else if (phi < -M_PI)
        phi += M_PI;
      // (ii) theta
      theta += deltaTheta;
      if (theta > M_PI)
        theta -= M_PI;
      else if (theta < 0.)
        theta += M_PI;

      double sphi = std::sin(phi);
      double cphi = std::cos(phi);
      double stheta = std::sin(theta);
      double ctheta = std::cos(theta);

      // assign the new values
      particle.setDirection(
          Acts::Vector3D(cphi * stheta, sphi * stheta, ctheta));
    } else {
      /// uniform distribution
      std::uniform_real_distribution<double> uniformDist(0., 1.);

      // Create a random uniform distribution between in the intervall [0,1]
      double psi = 2. * M_PI * uniformDist(generator);

      // more complex but "more true"
      Acts::Vector3D pDirection(particle.direction());
      double x = -pDirection.y();
      double y = pDirection.x();
      double z = 0.;

      // if it runs along the z axis - no good ==> take the x axis
      if (pDirection.z() * pDirection.z() > 0.999999) {
        x = 1.;
        y = 0.;
      }
      // deflector direction
      Acts::Vector3D deflector(x, y, z);
      // rotate the new direction for scattering using theta and  psi
      Acts::RotationMatrix3D rotation;
      rotation = Acts::AngleAxis3D(psi, pDirection) *
                 Acts::AngleAxis3D(angle3D, deflector);
      // rotate and set a new direction to the cache
      particle.setDirection(rotation * pDirection);
    }
    // scattering always returns an empty list
    // - it is a non-distructive process
    return {};
  }
};

}  // namespace ActsFatras
