// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <random>

#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsFatras/EventData/Particle.hpp"

namespace ActsFatras {
namespace detail {

/// Simulate (multiple) scattering using a configurable scattering model.
///
/// @tparam scattering_model_t Model implementation to draw a scattering angle.
template <typename scattering_model_t>
struct Scattering {
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
  std::array<Particle, 0> operator()(generator_t &generator,
                                     const Acts::MaterialProperties &slab,
                                     Particle &particle) const {
    // the scattered direction can be computed by rotating the initial
    // direction around a vector orthogonal to the initial direction, i.e. the
    // scattering deflector, by the scattering angle. there are an infinite
    // number of vectors orthogonal to the initial direction. the deflector is
    // rotated by some angle relative to some fixpoint.
    //
    // thus two random angles are required: the random deflector orientation
    // angle drawn uniformly from the [-pi,pi) range and the scattering angle
    // drawn from the specific scattering model distribution.

    // draw the random orientation angle
    const auto psi =
        std::uniform_real_distribution<double>(-M_PI, M_PI)(generator);
    // draw the scattering angle
    const auto theta = angle(generator, slab, particle);

    Acts::Vector3D direction = particle.unitDirection();
    // construct the combined rotation to the scattered direction
    Acts::RotationMatrix3D rotation(
        // rotation of the scattering deflector axis relative to the reference
        Acts::AngleAxis3D(psi, direction) *
        // rotation by the scattering angle around the deflector axis
        Acts::AngleAxis3D(theta, Acts::makeCurvilinearUnitU(direction)));
    direction.applyOnTheLeft(rotation);
    particle.setDirection(direction);

    // scattering is non-destructive and produces no secondaries
    return {};
  }
};

}  // namespace detail
}  // namespace ActsFatras
