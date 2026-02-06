// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialSlab.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Kernel/InteractionList.hpp"

#include <limits>

namespace ActsFatras {

/// A point like simulation process based on a physics model plus selectors.
///
/// @tparam physics_t is the physics model type
/// @tparam input_particle_selector_t is the input particle selector
/// @tparam output_particle_selector_t is the output particle selector
/// @tparam child_particle_selector_t is the child particle selector
///
/// The physics model type **must** provide a call operator with the following
/// two member functions
///
///     // generate X0/L0 limits
///     template <typename generator_t>
///     std::pair<Scalar, Scalar>
///     generatePathLimits(
///         generator& rng,
///         const Particle& particle) const
///
///     // run the process simulation
///     template <typename generator_t>
///     bool
///     run(
///         generator_t& rng,
///         Particle& particle,
///         std::vector<Particle>& generatedParticles) const
///
/// The return type of generatePathLimits() can be any `Container` with
/// `Particle` elements.
///
/// The input selector defines whether the process is applied while the
/// output selector defines a break condition, i.e. whether to continue
/// simulating the particle propagation. The child selector is used to
/// filter the generated child particles.
///
/// @note The output and child particle selectors are identical unless the
///       child particle selector is explicitly specified.
template <detail::PointLikeProcessConcept physics_t,
          typename input_particle_selector_t,
          typename output_particle_selector_t,
          typename child_particle_selector_t = output_particle_selector_t>
struct PointLikeProcess {
  /// The physics interactions implementation.
  physics_t physics;
  /// Input selection: if this process applies to this particle.
  input_particle_selector_t selectInputParticle;
  /// Output selection: if the particle is still valid after the interaction.
  output_particle_selector_t selectOutputParticle;
  /// Child selection: if a generated child particle should be kept.
  child_particle_selector_t selectChildParticle;

  /// Generate path limits for the process
  /// @param generator Random number generator
  /// @param particle Input particle
  /// @return Pair of minimum and maximum path limits
  template <class generator_t>
  std::pair<double, double> generatePathLimits(generator_t& generator,
                                               const Particle& particle) const {
    return physics.generatePathLimits(generator, particle);
  }

  /// Run the process
  /// @param rng Random number generator
  /// @param particle Input particle
  /// @param generatedParticles Output vector of generated particles
  /// @return True if process was applied
  template <class generator_t>
  bool run(generator_t& rng, Particle& particle,
           std::vector<Particle>& generatedParticles) const {
    if (!selectInputParticle(particle)) {
      return false;
    }

    std::vector<Particle> children;
    physics.run(rng, particle, children);

    std::copy_if(std::begin(children), std::end(children),
                 std::back_inserter(generatedParticles), selectChildParticle);

    return !selectOutputParticle(particle);
  }
};

}  // namespace ActsFatras
