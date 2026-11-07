// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/DefaultActor.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Kernel/SingleParticleSimulationResult.hpp"
#include "ActsFatras/Kernel/detail/SimulationActor.hpp"

#include <algorithm>
#include <cassert>
#include <memory>

namespace ActsFatras {

/// Single particle simulation with fixed propagator, interactions, and decay.
///
/// @tparam generator_t random number generator
/// @tparam interactions_t interaction list
/// @tparam hit_surface_selector_t selector for hit surfaces
/// @tparam decay_t decay module
template <typename propagator_t, typename interactions_t,
          typename hit_surface_selector_t, typename decay_t>
struct SingleParticleSimulation {
  /// How and within which geometry to propagate the particle.
  propagator_t propagator;
  /// Absolute maximum step size
  double maxStepSize = std::numeric_limits<double>::max();
  /// Absolute maximum path length
  double pathLimit = std::numeric_limits<double>::max();
  /// Decay module.
  decay_t decay;
  /// Interaction list containing the simulated interactions.
  interactions_t interactions;
  /// Selector for surfaces that should generate hits.
  hit_surface_selector_t selectHitSurface;
  /// Logger for debug output.
  std::unique_ptr<const Acts::Logger> logger;

  /// Alternatively construct the simulator with an external logger.
  /// @param propagator_ Propagator to use for particle simulation
  /// @param _logger Logger instance for debug output
  SingleParticleSimulation(propagator_t &&propagator_,
                           std::unique_ptr<const Acts::Logger> _logger)
      : propagator(propagator_), logger(std::move(_logger)) {}

  /// Simulate a single particle without secondaries.
  ///
  /// @tparam generator_t is the type of the random number generator
  ///
  /// @param geoCtx is the geometry context to access surface geometries
  /// @param magCtx is the magnetic field context to access field values
  /// @param generator is the random number generator
  /// @param particle is the initial particle state
  /// @returns Simulated particle state, hits, and generated particles.
  template <typename generator_t>
  Acts::Result<SingleParticleSimulationResult> simulate(
      const Acts::GeometryContext &geoCtx,
      const Acts::MagneticFieldContext &magCtx, generator_t &generator,
      const Particle &particle) const {
    // propagator-related additional types
    using Actor = detail::SimulationActor<generator_t, decay_t, interactions_t,
                                          hit_surface_selector_t>;
    using Result = typename Actor::state_type;
    using ActorList = Acts::ActorList<Acts::DefaultActor, Actor>;
    using PropagatorOptions =
        typename propagator_t::template Options<ActorList>;

    // Construct per-call options.
    PropagatorOptions options(geoCtx, magCtx);
    options.stepping.maxStepSize = maxStepSize;
    options.pathLimit = pathLimit;
    // setup the interactor as part of the propagator options
    auto &actor = options.actor.template get<Actor>();
    actor.generator = &generator;
    actor.decay = decay;
    actor.interactions = interactions;
    actor.selectHitSurface = selectHitSurface;
    actor.initialParticle = particle;

    auto state = propagator.makeState(options);
    auto result =
        particle.hasReferenceSurface()
            ? propagator.propagate(state,
                                   particle.boundParameters(geoCtx).value())
            : propagator.propagate(state, particle.curvilinearParameters());
    if (!result.ok()) {
      return result.error();
    }
    return std::move(state.actorState.template get<Result>());
  }
};

}  // namespace ActsFatras
