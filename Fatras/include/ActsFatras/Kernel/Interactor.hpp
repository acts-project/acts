// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <climits>
#include <cmath>
#include <sstream>
#include <utility>

#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsFatras/Kernel/PhysicsList.hpp"

namespace ActsFatras {

struct VoidSelector {
  bool operator()(const Acts::Surface &) const { return false; }
};

/// The Fatras Interactor
///
/// This is the Fatras plugin to the ACTS Propagator, it replaces
/// the MaterialInteractor of the reconstruction
///
/// @tparam generator_t Type of the random generator
/// @tparam particle_t is Type of the particle
/// @tparam hit_t Type of the simulated hit
/// @tparam hit_creator_t Type of the hit creator (does thruth association)
/// @tparam sensitive_selector_t The Selector type to identify sensitive
/// surfaces
/// @tparam physics_list_t Type of Extendable physics list that is called
/// @tparam decay_list_t Type of Extendable decay list that is called
///
/// The physics list plays a central role in this DetectorInteractor
/// it is called on each process that is defined at compile time
/// if a process triggers an abort, this will be forwarded to
/// the propagation cache.
template <typename generator_t, typename particle_t, typename hit_t,
          typename hit_creator_t, typename sensitive_selector_t = VoidSelector,
          typename physics_list_t = PhysicsList<>>
struct Interactor {
  using PhysicsList_t = physics_list_t;

  /// The random generator to be spawnper event
  generator_t *generator = nullptr;

  /// The slector for sensitive surfaces
  sensitive_selector_t sensitiveSelector;

  /// The physics list provided for this call
  physics_list_t physicsList;

  /// Simple result struct to be returned
  particle_t initialParticle;

  /// The hit creator helper class
  hit_creator_t hitCreator;

  /// It mainly acts as an internal state cache which is
  /// created for every propagation/extrapolation step
  struct this_result {
    /// result initialization
    bool initialized = false;

    /// The current particle - updated along the way
    particle_t particle;

    /// The outgoing particles due to physics processes
    std::vector<particle_t> outgoing;

    /// The simulated hits created along the way
    std::vector<hit_t> simulatedHits;
  };

  typedef this_result result_type;

  /// Interaction with detector material for the ActionList of the Propagator
  ///
  /// It checks if the cache has a current surface, in which case the action
  /// is performed according to the physics list content.
  ///
  /// Eventual particles produced in electromagnetic or hadronic interactions
  /// are stored in the result struct and can thus be retrieved by the caller
  ///
  /// @tparam  propagator_state_t is the type of Propagtor state
  /// @tparam  stepper_t the type of the Stepper for the access to the state
  ///
  /// @param state is the mutable propagator state object
  /// @param stepper is the propagation stepper object
  /// @param result is the mutable result cache object
  ///
  /// return value is void as it is a standard actor in the
  /// propagation
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t &state, stepper_t &stepper,
                  result_type &result) const {
    // If we are on target, everything should have been done
    if (state.navigation.targetReached)
      return;

    // Initialize the result, the state is thread local
    if (!result.initialized) {
      // set the initial particle parameters
      result.particle = initialParticle;
      result.initialized = true;
    }
    // get position and momentum presetp
    auto position = stepper.position(state.stepping);
    auto direction = stepper.direction(state.stepping);
    auto p = stepper.momentum(state.stepping);

    // set the stepping position to the particle
    result.particle.update(position, p * direction, 0., 0.,
                           stepper.time(state.stepping));

    // Check if the current surrface a senstive one
    bool sensitive = state.navigation.currentSurface
                         ? sensitiveSelector(*state.navigation.currentSurface)
                         : false;
    double depositedEnergy = 0.;

    // a current surface has been assigned by the navigator
    if (state.navigation.currentSurface &&
        state.navigation.currentSurface->surfaceMaterial()) {
      // get the surface material and the corresponding material properties
      auto sMaterial = state.navigation.currentSurface->surfaceMaterial();
      const Acts::MaterialProperties &mProperties =
          sMaterial->materialProperties(position);
      bool breakIndicator = false;
      if (mProperties) {
        // run the Fatras physics list - only when there's material
        breakIndicator = physicsList(*generator, mProperties, result.particle,
                                     result.outgoing);
      }
    }
    // Update the stepper cache with the current particle parameters
    position = result.particle.position();
    direction = result.particle.momentum().normalized();
    stepper.update(state.stepping, position, direction, result.particle.p(),
                   result.particle.time());
    // create the hit on a senstive element
    if (sensitive) {
      // create and fill the hit
      double htime =
          stepper.time(state.stepping);  //!< todo calculate from delta time
      hit_t simHit =
          hitCreator(*state.navigation.currentSurface, position, direction,
                     depositedEnergy, htime, result.particle);
      result.simulatedHits.push_back(std::move(simHit));
    }
  }

  /// Pure observer interface
  /// This does not apply to the Fatras simulator
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t &, stepper_t &) const {}
};

}  // namespace ActsFatras
