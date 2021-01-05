// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Kernel/SimulationResult.hpp"

#include <cassert>
#include <limits>

namespace ActsFatras {
namespace detail {

/// Fatras simulator plugin for the Acts propagator.
///
/// This plugin must be added to the action list of the propagator and is the
/// equivalent to the `MaterialInteractor` for the reconstruction. This
/// implements surface-based simulation of particle interactions with matter
/// using a configurable physics lists as well as some parts of the decay
/// simulation. The physics lists is called for every surface with valid
/// material.
///
/// @tparam generator_t is a random number generator
/// @tparam physics_list_t is a simulation physics lists
/// @tparam hit_surface_selector_t is a selector of sensitive hit surfaces
template <typename generator_t, typename physics_list_t,
          typename hit_surface_selector_t>
struct Interactor {
  using result_type = SimulationResult;

  /// Abort if the particle was killed during a previous interaction.
  struct ParticleNotAlive {
    // This references the Interactor to automatically access its result type.
    using action_type = Interactor;

    template <typename propagator_state_t, typename stepper_t>
    constexpr bool operator()(propagator_state_t &, const stepper_t &,
                              const result_type &result) const {
      return (result.particleStatus != SimulationParticleStatus::eAlive);
    }
  };

  /// Random number generator used for the simulation.
  generator_t *generator = nullptr;
  /// Physics list detailing the simulated interactions and processes.
  physics_list_t physics;
  /// Selector for surfaces that should generate hits.
  hit_surface_selector_t selectHitSurface;
  /// Initial particle state.
  Particle initialParticle;
  Particle::Scalar properTimeLimit =
      std::numeric_limits<Particle::Scalar>::infinity();

  /// Simulate the interaction with a single surface.
  ///
  /// @tparam propagator_state_t is propagator state
  /// @tparam stepper_t is the stepper instance
  ///
  /// @param state is the mutable propagator state object
  /// @param stepper is the propagation stepper object
  /// @param result is the mutable result/cache object
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t &state, stepper_t &stepper,
                  result_type &result) const {
    assert(generator and "The generator pointer must be valid");

    // compute the change in proper time first. this needs the particle state
    // from the previous step for reference and thus must occur before the state
    // is updated. this also means, that the first simulation step needs to be
    // treated differently from subsequent ones as there is no previous state.
    // for the decay, only the elapsed proper time during simulation is
    // relevant; the absolute time does matter.
    //
    // a particle can loose energy and thus its gamma factor is not a constant
    // of motion. since the stepper provides only the lab time, we need to
    // compute the change in proper time for each step separately. this assumes
    // that the gamma factor is constant over one stepper step.
    if (not std::isfinite(result.properTime)) {
      // compute the first step w/ respect to the initial particle
      auto deltaLabTime = stepper.time(state.stepping) - initialParticle.time();
      // proper-time = time / gamma = (1/gamma) * time
      //       beta² = p²/E²
      //       gamma = 1 / sqrt(1 - beta²) = sqrt(m² + p²) / m
      //     1/gamma = m / sqrt(m² + p²) = m / E
      auto gammaInv = initialParticle.mass() / initialParticle.energy();
      // first step resets the proper time
      result.properTime = gammaInv * deltaLabTime;
    } else {
      // compute the current step w/ respect to the previous state
      auto deltaLabTime = stepper.time(state.stepping) - result.particle.time();
      auto gammaInv = result.particle.mass() / result.particle.energy();
      // all other steps accumulate proper time
      result.properTime += gammaInv * deltaLabTime;
    }

    // computing the proper time above always requires the particle state from
    // the previous stepper state. thus, it must be updated unconditionally
    // regardless of whether any other simulation actions are taken below.
    result.particle = makeParticle(stepper, state.stepping, result);

    // decay check
    // TODO limit the stepsize when close to the lifetime limit to avoid
    //   overstepping and decaying the particle systematically too late
    if (properTimeLimit < result.properTime) {
      // result.particle was already updated
      result.particleStatus = SimulationParticleStatus::eDecayed;
      return;
    }

    // If we are on target, everything should have been done
    if (state.navigation.targetReached) {
      return;
    }
    // If we are not on a surface, there is nothing further for us to do
    if (not state.navigation.currentSurface) {
      return;
    }
    const Acts::Surface &surface = *state.navigation.currentSurface;
    // we want to keep the particle state before and after the interaction.
    // since the particle is modified in-place we need a copy.
    const Particle &before = result.particle;
    Particle after = result.particle;

    // interactions only make sense if there is material to interact with.
    if (surface.surfaceMaterial()) {
      // TODO - is this the right thing to do when globalToLocal fails
      // it should in principle never happen, so probably it would be best
      // to change to a model using transform() directly
      auto lpResult = surface.globalToLocal(state.geoContext, before.position(),
                                            before.unitDirection());
      if (lpResult.ok()) {
        Acts::Vector2 local = lpResult.value();
        Acts::MaterialSlab slab =
            surface.surfaceMaterial()->materialSlab(local);

        // again: no valid material -> no interaction
        if (slab) {
          // adapt material for non-zero incidence
          auto normal = surface.normal(state.geoContext, local);
          // dot-product(unit normal, direction) = cos(incidence angle)
          // particle direction is normalized, not sure about surface normal
          auto cosIncidenceInv =
              normal.norm() / normal.dot(before.unitDirection());
          slab.scaleThickness(cosIncidenceInv);
          // physics list returns true if the particle was killed
          if (physics(*generator, slab, after, result.generatedParticles)) {
            result.particleStatus = SimulationParticleStatus::eInteracted;
          }
          // add the accumulated material; assumes the full material was passsed
          // event if the particle was killed.
          result.pathInX0 += slab.thicknessInX0();
          result.pathInL0 += slab.thicknessInL0();
        }
      }
    }

    // store results of this interaction step, including potential hits
    result.particle = after;
    if (selectHitSurface(surface)) {
      result.hits.emplace_back(
          surface.geometryId(), before.particleId(),
          // the interaction could potentially modify the particle position
          Hit::Scalar(0.5) * (before.fourPosition() + after.fourPosition()),
          before.fourMomentum(), after.fourMomentum(), result.hits.size());
    }

    // continue the propagation with the modified parameters
    stepper.update(state.stepping, after.position(), after.unitDirection(),
                   after.absoluteMomentum(), after.time());
  }

  /// Pure observer interface. Does not apply to the Fatras simulator.
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t &, stepper_t &) const {}

  /// Construct the particle state from the stepper state.
  template <typename stepper_t>
  Particle makeParticle(const stepper_t &stepper,
                        const typename stepper_t::State &state,
                        const result_type &result) const {
    // avoid having a clumsy `initialized` flag by reconstructing the particle
    // state directly from the propagation state; using only the identity
    // parameters from the initial particle state.
    return Particle(initialParticle)
        .setPosition4(stepper.position(state), stepper.time(state))
        .setDirection(stepper.direction(state))
        .setAbsoluteMomentum(stepper.momentum(state));
  }
};

}  // namespace detail
}  // namespace ActsFatras
