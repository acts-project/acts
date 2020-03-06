// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>
#include <memory>

#include "Acts/EventData/ChargePolicy.hpp"
#include "Acts/EventData/SingleCurvilinearTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Kernel/Interactor.hpp"

namespace ActsFatras {

/// Single particle simulator with a fixed propagator and physics list.
///
/// @tparam propagator_t is the type of the underlying propagator
/// @tparam physics_list_t is the type of the simulated physics list
/// @tparam hit_surface_selector_t is the type that selects hit surfaces
template <typename propagator_t, typename physics_list_t,
          typename hit_surface_selector_t>
struct ParticleSimulator {
  /// How and within which geometry to propagate the particle.
  propagator_t propagator;
  /// What should be simulated. Will be copied to the per-call interactor.
  physics_list_t physics;
  /// Where hits are registiered. Will be copied to the per-call interactor.
  hit_surface_selector_t selectHitSurface;
  /// Local logger for debug output.
  std::shared_ptr<const Acts::Logger> localLogger = nullptr;

  /// Construct the simulator with the underlying propagator.
  ParticleSimulator(propagator_t &&propagator_, Acts::Logging::Level lvl)
      : propagator(propagator_),
        localLogger(Acts::getDefaultLogger("Simulator", lvl)) {}

  /// Provide access to the local logger instance, e.g. for logging macros.
  const Acts::Logger &logger() const { return *localLogger; }

  /// Simulate a single particle without secondaries.
  ///
  /// @param geoCtx is the geometry context to access surface geometries
  /// @param magCtx is the magnetic field context to access field values
  /// @param generator is the random number generator
  /// @param particle is the initial particle state
  /// @returns the result of the corresponding Interactor propagator action.
  ///
  /// @tparam generator_t is the type of the random number generator
  template <typename generator_t>
  Acts::Result<typename Interactor<generator_t, physics_list_t,
                                   hit_surface_selector_t>::result_type>
  simulate(const Acts::GeometryContext &geoCtx,
           const Acts::MagneticFieldContext &magCtx, generator_t &generator,
           const Particle &particle) const {
    assert(localLogger and "Missing local logger");

    // propagator-related additional types
    using Interact =
        Interactor<generator_t, physics_list_t, hit_surface_selector_t>;
    using InteractResult = typename Interact::result_type;
    using Actions = Acts::ActionList<Interact, Acts::detail::DebugOutputActor>;
    using Abort = Acts::AbortList<typename Interact::ParticleNotAlive,
                                  Acts::detail::EndOfWorldReached>;
    using PropagatorOptions = Acts::PropagatorOptions<Actions, Abort>;

    // Construct per-call options.
    PropagatorOptions options(geoCtx, magCtx);
    options.absPdgCode = particle.pdg();
    options.mass = particle.mass();
    options.debug = localLogger->doPrint(Acts::Logging::Level::DEBUG);
    // setup the interactor as part of the propagator options
    auto &interactor = options.actionList.template get<Interact>();
    interactor.generator = &generator;
    interactor.physics = physics;
    interactor.selectHitSurface = selectHitSurface;
    interactor.particle = particle;

    // run with a start parameter type depending on the particle charge.
    // TODO make track parameters consistently constructible regardless
    //      of the charge policy and template the class on the parameter.
    if (particle.charge() != 0) {
      Acts::SingleCurvilinearTrackParameters<Acts::ChargedPolicy> start(
          std::nullopt, particle.position(),
          particle.absMomentum() * particle.unitDirection(), particle.charge(),
          particle.time());
      auto result = propagator.propagate(start, options);
      if (result.ok()) {
        return result.value().template get<InteractResult>();
      } else {
        return result.error();
      }
    } else {
      Acts::SingleCurvilinearTrackParameters<Acts::NeutralPolicy> start(
          std::nullopt, particle.position(),
          particle.absMomentum() * particle.unitDirection(), particle.time());
      auto result = propagator.propagate(start, options);
      if (result.ok()) {
        return result.value().template get<InteractResult>();
      } else {
        return result.error();
      }
    }
  }
};

/// Fatras multi-particle simulator.
///
/// @tparam charged_selector_t Callable selector type for charged particles
/// @tparam charged_simulator_t Single particle simulator for charged particles
/// @tparam neutral_selector_t Callable selector type for neutral particles
/// @tparam neutral_simulator_t Single particle simulator for neutral particles
template <typename charged_selector_t, typename charged_simulator_t,
          typename neutral_selector_t, typename neutral_simulator_t>
struct Simulator {
  charged_selector_t selectCharged;
  neutral_selector_t selectNeutral;
  charged_simulator_t charged;
  neutral_simulator_t neutral;

  /// Construct from the single charged/neutral particle simulators.
  Simulator(charged_simulator_t &&charged_, neutral_simulator_t &&neutral_)
      : charged(std::move(charged_)), neutral(std::move(neutral_)) {}

  /// Simulate multiple particles and generated secondaries.
  ///
  /// @param geoCtx is the geometry context to access surface geometries
  /// @param magCtx is the magnetic field context to access field values
  /// @param inputParticles contains all particles that should be simulated
  /// @param simulatedParticlesInitial contains initial particle states
  /// @param simulatedParticlesFinal contains final particle states
  /// @param hits contains all generated hits
  ///
  /// This takes all input particles and simulates those passing the selection
  /// using the appropriate simulator. All selected particle states including
  /// additional ones generated from interactions are stored in separate output
  /// containers; both the initial state at the production vertex and the final
  /// state after propagation are stored. Hits generated from selected input and
  /// generated particles are stored in the hit container.
  ///
  /// @warning Particle-hit association is based on particle ids generated
  ///          during the simulation. This requires that all input particles
  ///          **must** have generation and sub-particle number set to zero.
  ///
  /// @tparam generator_t is the type of the random number generator
  /// @tparam input_particles_t is a Container for particles
  /// @tparam output_particles_t is a SequenceContainer for particles
  /// @tparam hits_t is a SequenceContainer for hits
  template <typename generator_t, typename input_particles_t,
            typename output_particles_t, typename hits_t>
  Acts::Result<void> simulate(const Acts::GeometryContext &geoCtx,
                              const Acts::MagneticFieldContext &magCtx,
                              generator_t &generator,
                              const input_particles_t &inputParticles,
                              output_particles_t &simulatedParticlesInitial,
                              output_particles_t &simulatedParticlesFinal,
                              hits_t &hits) const {
    assert(
        (simulatedParticlesInitial.size() == simulatedParticlesFinal.size()) and
        "Inconsistent initial sizes of the simulated particle containers");

    for (const Particle &inputParticle : inputParticles) {
      if (not selectParticle(inputParticle)) {
        continue;
      }
      // required to allow correct particle id numbering for secondaries later
      if ((inputParticle.particleId().generation() != 0u) or
          (inputParticle.particleId().subParticle() != 0u)) {
        // TODO add meaningfull error code
        return std::error_code(-1, std::generic_category());
      }

      // Do a *depth-first* simulation of the particle and its secondaries,
      // i.e. we simulate all secondaries, tertiaries, ... before simulating
      // the next primary particle. Use the end of the output container as
      // a queue to store particles that should be simulated.
      //
      // WARNING the initial particle output container will be modified during
      //         iteration as new secondaries are added directly to it. to avoid
      //         issues, access must always occur via indices.
      auto iinitial = simulatedParticlesInitial.size();
      simulatedParticlesInitial.push_back(inputParticle);
      for (; iinitial < simulatedParticlesInitial.size(); ++iinitial) {
        // only simulatable particles are pushed to the container
        // no additional check is necessary here.
        const auto &initialParticle = simulatedParticlesInitial[iinitial];

        if (selectCharged(initialParticle)) {
          auto result =
              charged.simulate(geoCtx, magCtx, generator, initialParticle);
          if (not result.ok()) {
            // do not keep unsimulated/ failed particles in the output
            simulatedParticlesInitial.resize(iinitial);
            return result.error();
          }
          copyOutputs(result.value(), simulatedParticlesInitial,
                      simulatedParticlesFinal, hits);
        } else if (selectNeutral(initialParticle)) {
          auto result =
              neutral.simulate(geoCtx, magCtx, generator, initialParticle);
          if (not result.ok()) {
            // do not keep unsimulated/ failed particles in the output
            simulatedParticlesInitial.resize(iinitial);
            return result.error();
          }
          copyOutputs(result.value(), simulatedParticlesInitial,
                      simulatedParticlesFinal, hits);
        }

        // since physics processes are independent, there can be particle id
        // collisions within the generated secondaries. they can be resolved by
        // renumbering within each sub-particle generation. this must happen
        // before the particle is simulated since the particle id is used to
        // associate generated hits back to the particle.
        renumberTailParticleIds(simulatedParticlesInitial, iinitial);
      }
    }

    return Acts::Result<void>::success();
  }

 private:
  // Select if the particle should be simulated at all.
  //
  // This also enforces mutual-exclusivity of the two charge selections. If both
  // charge selections evaluate true, they are probably not setup correctly and
  // not simulating them at all is a reasonable fall-back.
  bool selectParticle(const Particle &particle) const {
    const bool isValidCharged = selectCharged(particle);
    const bool isValidNeutral = selectNeutral(particle);
    assert(not(isValidCharged and isValidNeutral) and
           "Charge selectors are not mutually exclusive");
    return isValidCharged xor isValidNeutral;
  }

  /// Copy Interactor results to output containers.
  ///
  /// @tparam interactor_result_t is an Interactor result struct
  /// @tparam particles_t is a SequenceContainer for particles
  /// @tparam hits_t is a SequenceContainer for hits
  template <typename interactor_result_t, typename particles_t, typename hits_t>
  void copyOutputs(const interactor_result_t &result,
                   particles_t &particlesInitial, particles_t &particlesFinal,
                   hits_t &hits) const {
    // initial particle state was already pushed to the container before
    // store final particle state at the end of the simulation
    particlesFinal.push_back(result.particle);
    // move generated secondaries that should be simulated to the output
    std::copy_if(
        result.generatedParticles.begin(), result.generatedParticles.end(),
        std::back_inserter(particlesInitial),
        [this](const Particle &particle) { return selectParticle(particle); });
    std::copy(result.hits.begin(), result.hits.end(), std::back_inserter(hits));
  }

  /// Renumber particle ids in the tail of the container.
  ///
  /// Ensures particle ids are unique by modifying the sub-particle number
  /// within each generation.
  ///
  /// @param particles particle container in which particles are renumbered
  /// @param lastValid index of the last particle with a valid particle id
  ///
  /// @tparam particles_t is a SequenceContainer for particles
  ///
  /// @note This function assumes that all particles in the tail have the same
  ///       vertex numbers (primary/secondary) and particle number and are
  ///       ordered according to their generation number.
  ///
  template <typename particles_t>
  static void renumberTailParticleIds(particles_t &particles,
                                      std::size_t lastValid) {
    // iterate over adjacent pairs; potentially modify the second element.
    // assume e.g. a primary particle 2 with generation=subparticle=0 that
    // generates two secondaries during simulation. we have the following
    // ids (decoded as vertex|particle|generation|subparticle)
    //
    //     v|2|0|0, v|2|1|0, v|2|1|0
    //
    // where v represents the vertex numbers. this will be renumbered to
    //
    //     v|2|0|0, v|2|1|0, v|2|1|1
    //
    // if each secondary then generates two tertiaries we could have e.g.
    //
    //     v|2|0|0, v|2|1|0, v|2|1|1, v|2|2|0, v|2|2|1, v|2|2|0, v|2|2|1
    //
    // which would then be renumbered to
    //
    //     v|2|0|0, v|2|1|0, v|2|1|1, v|2|2|0, v|2|2|1, v|2|2|2, v|2|2|3
    //
    for (auto j = lastValid; (j + 1u) < particles.size(); ++j) {
      const auto prevId = particles[j].particleId();
      auto currId = particles[j + 1u].particleId();
      // NOTE primary/secondary vertex and particle are assumed to be equal
      // only renumber within one generation
      if (prevId.generation() != currId.generation()) {
        continue;
      }
      // ensure the sub-particle is strictly monotonic within a generation
      if (prevId.subParticle() < currId.subParticle()) {
        continue;
      }
      // sub-particle numbering must be non-zero
      currId.setSubParticle(prevId.subParticle() + 1u);
      particles[j + 1u] = particles[j + 1u].withParticleId(currId);
    }
  }
};

}  // namespace ActsFatras
