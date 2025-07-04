// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Generators/EventGenerator.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Util.hpp"

#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <span>
#include <stdexcept>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/Print.h>

using namespace Acts::UnitLiterals;

namespace ActsExamples {

EventGenerator::EventGenerator(const Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("EventGenerator", lvl)) {
  if (m_cfg.generators.empty()) {
    throw std::invalid_argument("No generators are configured");
  }

  for (const auto& generator : m_cfg.generators) {
    if (generator.multiplicity == nullptr) {
      throw std::invalid_argument("Missing multiplicity generator");
    }
    if (generator.vertex == nullptr) {
      throw std::invalid_argument("Missing vertex generator");
    }
    if (generator.particles == nullptr) {
      throw std::invalid_argument("Missing particles generator");
    }
  }

  if (!m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers service");
  }

  m_outputEvent.maybeInitialize(m_cfg.outputEvent);
}

std::string EventGenerator::name() const {
  return "EventGenerator";
}

std::pair<std::size_t, std::size_t> EventGenerator::availableEvents() const {
  return {0u, std::numeric_limits<std::size_t>::max()};
}

ProcessCode EventGenerator::read(const AlgorithmContext& ctx) {
  ACTS_VERBOSE("EventGenerator::read");
  std::vector<SimParticle> particlesUnordered;
  std::vector<SimVertex> verticesUnordered;

  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);

  std::vector<std::shared_ptr<HepMC3::GenEvent>> genEvents;

  std::size_t nPrimaryVertices = 0;
  ACTS_VERBOSE("Using " << m_cfg.generators.size() << " generators");
  {
    Acts::AveragingScopedTimer genTimer("Generating primary vertices", logger(),
                                        Acts::Logging::DEBUG);

    for (std::size_t iGenerate = 0; iGenerate < m_cfg.generators.size();
         ++iGenerate) {
      auto& generate = m_cfg.generators[iGenerate];

      // generate the primary vertices from this generator
      std::size_t multiplicity = (*generate.multiplicity)(rng);
      ACTS_VERBOSE("Generating " << multiplicity << " primary vertices");
      for (std::size_t n = multiplicity; 0 < n; --n) {
        std::optional sample{genTimer.sample()};
        nPrimaryVertices += 1;

        // generate primary vertex position
        auto vertexPosition = (*generate.vertex)(rng);
        ACTS_VERBOSE("Generate vertex at " << vertexPosition.transpose());

        // generate particles associated to this vertex
        auto genEvent = (*generate.particles)(rng);
        if (genEvent->length_unit() != HepMC3::Units::MM) {
          throw std::runtime_error("EventGenerator: length unit is not mm");
        }
        if (genEvent->momentum_unit() != HepMC3::Units::GEV) {
          throw std::runtime_error("EventGenerator: momentum unit is not GeV");
        }

        ACTS_VERBOSE("Shifting event to " << vertexPosition.transpose());
        // Our internal time unit is ctau, so is HepMC3's, make sure we convert
        // to mm
        HepMC3::FourVector vtxPosHepMC(vertexPosition[Acts::eFreePos0] / 1_mm,
                                       vertexPosition[Acts::eFreePos1] / 1_mm,
                                       vertexPosition[Acts::eFreePos2] / 1_mm,
                                       vertexPosition[Acts::eFreeTime] / 1_mm);
        genEvent->shift_position_to(vtxPosHepMC);
        genEvent->add_attribute("acts",
                                std::make_shared<HepMC3::BoolAttribute>(true));

        if (m_cfg.printListing) {
          ACTS_VERBOSE("Generated event:\n"
                       << [&]() {
                            std::stringstream ss;
                            HepMC3::Print::listing(ss, *genEvent);
                            return ss.str();
                          }());
        }
        genEvents.push_back(genEvent);

        sample.reset();  // reset the gen timer
      }
    }
  }

  std::vector<const HepMC3::GenEvent*> eventPtrs;
  eventPtrs.reserve(genEvents.size());
  for (const auto& evt : genEvents) {
    eventPtrs.push_back(evt.get());
  }

  auto event = std::make_shared<HepMC3::GenEvent>();
  event->set_units(HepMC3::Units::GEV, HepMC3::Units::MM);

  HepMC3Util::mergeEvents(*event, eventPtrs, logger());
  event->set_event_number(static_cast<int>(ctx.eventNumber));

  ACTS_VERBOSE("Vertices size: " << event->vertices().size());
  if (m_cfg.printListing) {
    ACTS_VERBOSE("Final event:\n"
                 << [&]() {
                      std::stringstream ss;
                      HepMC3::Print::listing(ss, *event);
                      return ss.str();
                    }());
  }

  ACTS_DEBUG("event=" << ctx.eventNumber
                      << " n_primary_vertices=" << nPrimaryVertices
                      << " n_particles=" << event->particles().size());

  if (m_outputEvent.isInitialized()) {
    m_outputEvent(ctx, std::move(event));
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
