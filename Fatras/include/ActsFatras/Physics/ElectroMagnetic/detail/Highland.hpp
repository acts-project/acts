// This file is part of the Acts project.
//
// Copyright (C) 2018-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/Interactions.hpp"

#include <random>

namespace ActsFatras::detail {

/// Generate scattering angles using the Highland/PDG parametrization.
///
/// The angles are drawn from a single normal distribution with a width
/// given by the Highland/PDG formula.
struct Highland {
  /// Generate a single scattering angle.
  ///
  /// @param[in]     generator is the random number generator
  /// @param[in]     slab      defines the passed material
  /// @param[in,out] particle  is the particle being scattered
  /// @return the scattering angle
  ///
  /// @tparam generator_t is a RandomNumberEngine
  template <typename generator_t>
  double operator()(generator_t &generator, const Acts::MaterialSlab &slab,
                    Particle &particle) const {
    // compute the planar scattering angle
    const auto theta0 = Acts::computeMultipleScatteringTheta0(
        slab, particle.absolutePdg(), particle.mass(), particle.qOverP(),
        particle.absoluteCharge());
    return std::normal_distribution<double>(0.0, theta0)(generator);
  }
};

}  // namespace ActsFatras::detail
