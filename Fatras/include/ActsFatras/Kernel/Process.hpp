// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialProperties.hpp"
#include "ActsFatras/EventData/Particle.hpp"

namespace ActsFatras {

/// No-op particle selector that selects all particles.
struct EveryParticle {
  constexpr bool operator()(const Particle &) const { return true; }
};

// No-op input selector that selects all configurations.
struct EveryInput {
  constexpr bool operator()(const Particle &,
                            const Acts::MaterialProperties &) const {
    return true;
  }
};

/// Enable usage of a particle selector as an input selector.
template <typename ParticleSelector>
struct AsInputSelector : public ParticleSelector {
  bool operator()(const Particle &particle,
                  const Acts::MaterialProperties &) const {
    return ParticleSelector::operator()(particle);
  }
};

/// A simulation process based on a physics interaction plus selectors.
///
/// @tparam physics_t is the physics interaction type
/// @tparam input_selector_t is the input material + particle selector
/// @tparam output_particle_selector_t is the output particle selector
/// @tparam child_particle_selector_t is the child particle selector
///
/// The input selector defines whether the interaction is applied while the
/// output selector defines a break condition, i.e. whether to continue
/// simulating the particle propagation. The child selector is used to
/// filter the generated child particles.
///
/// @note The output and child particle selectors are identical unless the
///       child particle selector is explicitely specified.
template <typename physics_t, typename input_selector_t = EveryInput,
          typename output_particle_selector_t = EveryParticle,
          typename child_particle_selector_t = output_particle_selector_t>
struct Process {
  /// The physics interactions implementation.
  physics_t physics;
  /// Input selection: if this process applies to material + particle.
  input_selector_t selectInput;
  /// Output selection: if the particle is still valid after the interaction.
  output_particle_selector_t selectOutputParticle;
  /// Child selection: if a generated child particle should be kept.
  child_particle_selector_t selectChildParticle;

  /// Execute the physics process considering the configured selectors.
  ///
  /// @param[in]     generator is the random number generator
  /// @param[in]     slab      is the passed material
  /// @param[in,out] particle  is the particle being updated
  /// @param[out]    generated is the container of generated particles
  /// @return Break condition, i.e. whether this process stops the propagation
  ///
  /// @tparam generator_t must be a RandomNumberEngine
  template <typename generator_t>
  bool operator()(generator_t &generator, const Acts::MaterialProperties &slab,
                  Particle &particle, std::vector<Particle> &generated) const {
    // not selecting this process is not a break condition
    if (not selectInput(particle, slab)) {
      return false;
    }
    // modify particle according to the physics process
    auto children = physics(generator, slab, particle);
    // move selected child particles to the output container
    std::copy_if(children.begin(), children.end(),
                 std::back_inserter(generated), selectChildParticle);
    // break condition is defined by whether the output particle is still valid
    // or not e.g. because it has fallen below a momentum threshold.
    return not selectOutputParticle(particle);
  }
};

}  // namespace ActsFatras
