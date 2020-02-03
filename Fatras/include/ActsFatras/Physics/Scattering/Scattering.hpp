// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsFatras/Kernel/detail/RandomNumberDistributions.hpp"

namespace ActsFatras {

/// @class Scatterering
///
/// This is the (multiple) scattering plugin to the
/// Physics list. It needs a scattering formula in order
/// to provide the the scattering angle in 3D space.
///
/// There's two options to apply the scattering
/// - a parametric action that relates phi and theta (default: off)
/// - an actuall out of direction scattering applying two random numbers
template <typename formula_t>
struct Scattering {
  /// The flag to include scattering or not
  bool scattering = true;

  /// Include the log term
  bool parametric = false;
  double projectionFactor = 1. / std::sqrt(2.);

  /// The scattering formula
  formula_t angle;

  /// This is the scattering call operator
  template <typename generator_t, typename detector_t, typename particle_t>
  std::vector<particle_t> operator()(generator_t &gen, const detector_t &det,
                                     particle_t &in) const {
    // Do nothing if the flag is set to false
    if (not scattering) {
      return {};
    }

    // 3D scattering angle
    double angle3D = angle(gen, det, in);

    // parametric scattering
    if (parametric) {
      // the initial values
      double theta = Acts::VectorHelpers::theta(in.momentum());
      double phi = Acts::VectorHelpers::phi(in.momentum());
      double sinTheta = (sin(theta) * sin(theta) > 10e-10) ? sin(theta) : 1.;

      // sample them in an independent way
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
      in.scatter(in.p() * Acts::Vector3D(cphi * stheta, sphi * stheta, ctheta));
    } else {
      /// uniform distribution
      UniformDist uniformDist = UniformDist(0., 1.);

      // Create a random uniform distribution between in the intervall [0,1]
      double psi = 2. * M_PI * uniformDist(gen);

      // more complex but "more true"
      Acts::Vector3D pDirection(in.momentum().normalized());
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
      in.scatter(in.p() * rotation * pDirection.normalized());
    }
    // scattering always returns an empty list
    // - it is a non-distructive process
    return {};
  }
};

}  // namespace ActsFatras
