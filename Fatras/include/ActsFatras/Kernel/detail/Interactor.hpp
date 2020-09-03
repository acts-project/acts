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

namespace ActsFatras {
namespace detail {

/// Fatras interactor plugin for the Acts propagator.
///
/// This plugin must be added to the action list of the propagator and is the
/// equivalent to the `MaterialInteractor` for the reconstruction. This
/// implements surface-based simulation of particle interactions with matter
/// using a configurable physics lists. The physics lists is called for
/// every surface with valid material.
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
      return not result.isAlive;
    }
  };

  /// Random number generator used for the simulation.
  generator_t *generator = nullptr;
  /// Physics list detailing the simulated interactions and processes.
  physics_list_t physics;
  /// Selector for surfaces that should generate hits.
  hit_surface_selector_t selectHitSurface;
  /// Initial particle state.
  Particle particle;

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

    // If we are on target, everything should have been done
    if (state.navigation.targetReached) {
      return;
    }
    // If we are not on a surface, there is nothing for us to do
    if (not state.navigation.currentSurface) {
      return;
    }
    const Acts::Surface &surface = *state.navigation.currentSurface;

    // avoid having a clumsy `initialized` flag by reconstructing the particle
    // state directly from the propagation state; using only the identity
    // parameters from the initial particle state.
    const auto before =
        Particle(particle)
            // include passed material from the initial particle state
            .setMaterialPassed(particle.pathInX0() + result.pathInX0,
                               particle.pathInL0() + result.pathInL0)
            .setPosition4(stepper.position(state.stepping),
                          stepper.time(state.stepping))
            .setDirection(stepper.direction(state.stepping))
            .setAbsMomentum(stepper.momentum(state.stepping));
    // we want to keep the particle state before and after the interaction.
    // since the particle is modified in-place we need a copy.
    auto after = before;

    // interactions only make sense if there is material to interact with.
    if (surface.surfaceMaterial()) {
      // TODO - is this the right thing to do when globalToLocal fails
      // it should in principle never happen, so probably it would be best
      // to change to a model using transform() directly
      auto lpResult = surface.globalToLocal(state.geoContext, before.position(),
                                            before.unitDirection());
      if (lpResult.ok()) {
        Acts::Vector2D local(0., 0.);

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
          // physics list returns if the particle was killed.
          result.isAlive =
              not physics(*generator, slab, after, result.generatedParticles);
          // add the accumulated material; assumes the full material was passsed
          // event if the particle was killed.
          result.pathInX0 += slab.thicknessInX0();
          result.pathInL0 += slab.thicknessInL0();
          // WARNING this overwrites changes that the physics interactions
          //         might have performed with regard to the passed material.
          //         ensures consistent material counting by making the one
          //         component that by construction will see all material
          //         contributions (this Interactor) responsible.
          // TODO review this for supporting multiple interactions within the
          // same material slab
          after.setMaterialPassed(before.pathInX0() + slab.thicknessInX0(),
                                  before.pathInL0() + slab.thicknessInL0());
        }
      }
    }

    // store results of this interaction step, including potential hits
    result.particle = after;
    if (selectHitSurface(surface)) {
      result.hits.emplace_back(
          surface.geometryId(), before.particleId(),
          // the interaction could potentially modify the particle position
          Hit::Scalar(0.5) * (before.position4() + after.position4()),
          before.momentum4(), after.momentum4(), result.hits.size());
    }

    // continue the propagation with the modified parameters
    stepper.update(state.stepping, after.position(), after.unitDirection(),
                   after.absMomentum(), after.time());
  }

  /// Pure observer interface. Does not apply to the Fatras simulator.
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t &, stepper_t &) const {}
};

}  // namespace detail
}  // namespace ActsFatras
