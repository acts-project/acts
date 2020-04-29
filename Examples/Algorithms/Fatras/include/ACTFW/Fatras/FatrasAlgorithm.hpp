// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <string>

#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

namespace FW {

/// Fast track simulation using the Acts propagation and navigation.
///
/// @tparam simulator_t the Fatras simulation kernel type
template <typename simulator_t>
class FatrasAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// The particles input collection.
    std::string inputParticles;
    /// The simulated particles initial state collection.
    std::string outputParticlesInitial;
    /// The simulated particles final state collection.
    std::string outputParticlesFinal;
    /// The simulated hits output collection.
    std::string outputHits;
    /// The simulator kernel.
    simulator_t simulator;
    /// Random number service.
    std::shared_ptr<const RandomNumbers> randomNumbers;

    /// Construct the algorithm config with the simulator kernel.
    Config(simulator_t&& simulator_) : simulator(std::move(simulator_)) {}
  };

  /// Construct the algorithm from a config.
  ///
  /// @param cfg is the configuration struct
  /// @param lvl is the logging level
  FatrasAlgorithm(Config cfg, Acts::Logging::Level lvl)
      : FW::BareAlgorithm("FatrasAlgorithm", lvl), m_cfg(std::move(cfg)) {
    ACTS_DEBUG("hits on sensitive surfaces: "
               << m_cfg.simulator.charged.selectHitSurface.sensitive);
    ACTS_DEBUG("hits on material surfaces: "
               << m_cfg.simulator.charged.selectHitSurface.material);
    ACTS_DEBUG("hits on passive surfaces: "
               << m_cfg.simulator.charged.selectHitSurface.passive);
  }

  /// Run the simulation for a single event.
  ///
  /// @param ctx the algorithm context containing all event information
  FW::ProcessCode execute(const AlgorithmContext& ctx) const final override {
    // read input containers
    const auto& inputParticles =
        ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
    // prepare output containers
    SimParticleContainer::sequence_type particlesInitialUnordered;
    SimParticleContainer::sequence_type particlesFinalUnordered;
    SimHitContainer::sequence_type hitsUnordered;
    // reserve appropriate resources
    constexpr auto meanHitsPerParticle = 16u;
    particlesInitialUnordered.reserve(inputParticles.size());
    particlesFinalUnordered.reserve(inputParticles.size());
    hitsUnordered.reserve(meanHitsPerParticle * inputParticles.size());

    // run the simulation w/ a local random generator
    auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);
    auto ret = m_cfg.simulator.simulate(
        ctx.geoContext, ctx.magFieldContext, rng, inputParticles,
        particlesInitialUnordered, particlesFinalUnordered, hitsUnordered);
    // fatal error leads to panic
    if (not ret.ok()) {
      ACTS_FATAL("event " << ctx.eventNumber << " simulation failed with error "
                          << ret.error());
      return ProcessCode::ABORT;
    }
    // failed particles are just logged. assumes that failed particles are due
    // to edge-cases representing a tiny fraction of the event; not due to a
    // fundamental issue.
    for (const auto& failed : ret.value()) {
      ACTS_ERROR("event " << ctx.eventNumber << " particle " << failed.particle
                          << " failed to simulate with error " << failed.error
                          << ": " << failed.error.message());
    }
    // TODO is there a point where too many failed particles or failed particles
    //      of a particular type (e.g. from hard interaction or any primary
    //      particle) should also lead to a panic?

    ACTS_DEBUG(inputParticles.size() << " input particles");
    ACTS_DEBUG(particlesInitialUnordered.size()
               << " simulated particles (initial state)");
    ACTS_DEBUG(particlesFinalUnordered.size()
               << " simulated particles (final state)");
    ACTS_DEBUG(hitsUnordered.size() << " hits");

    // restore ordering for output containers
    SimParticleContainer particlesInitial;
    SimParticleContainer particlesFinal;
    SimHitContainer hits;
    particlesInitial.adopt_sequence(std::move(particlesInitialUnordered));
    particlesFinal.adopt_sequence(std::move(particlesFinalUnordered));
    hits.adopt_sequence(std::move(hitsUnordered));

    // store ordered output containers
    ctx.eventStore.add(m_cfg.outputParticlesInitial,
                       std::move(particlesInitial));
    ctx.eventStore.add(m_cfg.outputParticlesFinal, std::move(particlesFinal));
    ctx.eventStore.add(m_cfg.outputHits, std::move(hits));

    return FW::ProcessCode::SUCCESS;
  }

 private:
  Config m_cfg;
};

}  // namespace FW
