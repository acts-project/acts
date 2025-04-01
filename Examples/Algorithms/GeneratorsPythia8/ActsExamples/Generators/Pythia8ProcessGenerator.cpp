// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Generators/Pythia8ProcessGenerator.hpp"

#include "Acts/Utilities/MathHelpers.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <algorithm>
#include <iterator>
#include <ostream>
#include <random>
#include <utility>

#include <Pythia8/Pythia.h>

#if defined(_HAS_HEPMC3)
#include <HepMC3/WriterAscii.h>
#include <Pythia8Plugins/HepMC3.h>
#endif

namespace ActsExamples {

struct Pythia8RandomEngineWrapper : public Pythia8::RndmEngine {
  RandomEngine* rng{nullptr};

  struct {
    std::size_t numUniformRandomNumbers = 0;
    double first = std::numeric_limits<double>::quiet_NaN();
    double last = std::numeric_limits<double>::quiet_NaN();
  } statistics;

  Pythia8RandomEngineWrapper() = default;

  double flat() override {
    if (rng == nullptr) {
      throw std::runtime_error(
          "Pythia8RandomEngineWrapper: no random engine set");
    }

    double value = std::uniform_real_distribution<double>(0.0, 1.0)(*rng);
    if (statistics.numUniformRandomNumbers == 0) {
      statistics.first = value;
    }
    statistics.last = value;
    statistics.numUniformRandomNumbers++;
    return value;
  }

  void setRandomEngine(RandomEngine& rng_) { rng = &rng_; }
  void clearRandomEngine() { rng = nullptr; }
};

struct Pythia8GeneratorImpl {
#if defined(_HAS_HEPMC3)
  std::unique_ptr<HepMC3::Writer> m_hepMC3Writer;
  std::unique_ptr<HepMC3::Pythia8ToHepMC3> m_hepMC3Converter;
#endif
  std::shared_ptr<Pythia8RandomEngineWrapper> m_pythia8RndmEngine;
};

Pythia8Generator::Pythia8Generator(const Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("Pythia8Generator", lvl)),
      m_pythia8(std::make_unique<Pythia8::Pythia>("", false)) {
  // disable all output by default but allow re-enable via config
  m_pythia8->settings.flag("Print:quiet", true);
  for (const auto& setting : m_cfg.settings) {
    ACTS_VERBOSE("use Pythia8 setting '" << setting << "'");
    m_pythia8->readString(setting.c_str());
  }
  m_pythia8->settings.mode("Beams:idA", m_cfg.pdgBeam0);
  m_pythia8->settings.mode("Beams:idB", m_cfg.pdgBeam1);
  m_pythia8->settings.mode("Beams:frameType", 1);
  m_pythia8->settings.parm("Beams:eCM",
                           m_cfg.cmsEnergy / Acts::UnitConstants::GeV);

  m_impl = std::make_unique<Pythia8GeneratorImpl>();

  m_impl->m_pythia8RndmEngine = std::make_shared<Pythia8RandomEngineWrapper>();

#if PYTHIA_VERSION_INTEGER >= 8310
  m_pythia8->setRndmEnginePtr(m_impl->m_pythia8RndmEngine);
#else
  m_pythia8->setRndmEnginePtr(m_impl->m_pythia8RndmEngine.get());
#endif

  RandomEngine rng{m_cfg.initializationSeed};
  m_impl->m_pythia8RndmEngine->setRandomEngine(rng);
  m_pythia8->init();
  m_impl->m_pythia8RndmEngine->clearRandomEngine();

  if (m_cfg.writeHepMC3.has_value()) {
#if defined(_HAS_HEPMC3)
    ACTS_DEBUG("Initializing HepMC3 output to: " << m_cfg.writeHepMC3.value());
    m_impl->m_hepMC3Converter = std::make_unique<HepMC3::Pythia8ToHepMC3>();
    m_impl->m_hepMC3Writer =
        std::make_unique<HepMC3::WriterAscii>(m_cfg.writeHepMC3.value());
#else
    throw std::runtime_error("HepMC3 output requested but not available");
#endif
  }
}

// needed to allow unique_ptr of forward-declared Pythia class
Pythia8Generator::~Pythia8Generator() {
#if defined(_HAS_HEPMC3)
  if (m_impl->m_hepMC3Writer) {
    m_impl->m_hepMC3Writer->close();
  }
#endif

  ACTS_INFO("Pythia8Generator produced "
            << m_impl->m_pythia8RndmEngine->statistics.numUniformRandomNumbers
            << " uniform random numbers");
  ACTS_INFO("                 first = "
            << m_impl->m_pythia8RndmEngine->statistics.first);
  ACTS_INFO("                  last = "
            << m_impl->m_pythia8RndmEngine->statistics.last);
}

std::pair<SimVertexContainer, SimParticleContainer>
Pythia8Generator::operator()(RandomEngine& rng) {
  using namespace Acts::UnitLiterals;

  SimVertexContainer::sequence_type vertices;
  SimParticleContainer::sequence_type particles;

  // pythia8 is not thread safe and generation needs to be protected
  std::lock_guard<std::mutex> lock(m_pythia8Mutex);
  // use per-thread random engine also in pythia

  m_impl->m_pythia8RndmEngine->setRandomEngine(rng);

  {
    Acts::FpeMonitor mon{0};  // disable all FPEs while we're in Pythia8
    m_pythia8->next();
  }

  if (m_cfg.printShortEventListing) {
    m_pythia8->process.list();
  }
  if (m_cfg.printLongEventListing) {
    m_pythia8->event.list();
  }

#if defined(_HAS_HEPMC3)
  if (m_impl->m_hepMC3Converter && m_impl->m_hepMC3Writer) {
    auto hepmc_event = std::make_shared<HepMC3::GenEvent>();
    m_impl->m_hepMC3Converter->fill_next_event(*m_pythia8, hepmc_event.get());
    hepmc_event->set_units(HepMC3::Units::GEV, HepMC3::Units::MM);
    m_impl->m_hepMC3Writer->write_event(*hepmc_event);
  }
#endif

  // create the primary vertex
  vertices.emplace_back(SimVertexBarcode{0}, Acts::Vector4(0., 0., 0., 0.));

  // convert generated final state particles into internal format
  for (int ip = 0; ip < m_pythia8->event.size(); ++ip) {
    const auto& genParticle = m_pythia8->event[ip];

    // ignore beam particles
    if (genParticle.statusHepMC() == 4) {
      continue;
    }
    // only interested in final, visible particles
    if (!genParticle.isFinal()) {
      continue;
    }
    if (!genParticle.isVisible()) {
      continue;
    }

    // production vertex. Pythia8 time uses units mm/c, and we use c=1
    Acts::Vector4 pos4(genParticle.xProd() * 1_mm, genParticle.yProd() * 1_mm,
                       genParticle.zProd() * 1_mm, genParticle.tProd() * 1_mm);

    // define the particle identifier including possible secondary vertices

    SimBarcode particleId(0u);
    // ensure particle identifier component is non-zero
    particleId.setParticle(1u + particles.size());
    // only secondaries have a defined vertex position
    if (m_cfg.labelSecondaries && genParticle.hasVertex()) {
      // either add to existing secondary vertex if exists or create new one

      // check if an existing vertex is close enough
      auto it =
          std::ranges::find_if(vertices, [&pos4, this](const SimVertex& v) {
            return (pos4.head<3>() - v.position()).norm() <
                   m_cfg.spatialVertexThreshold;
          });

      if (it != vertices.end()) {
        particleId.setVertexSecondary(std::distance(vertices.begin(), it));
        it->outgoing.insert(particleId);
      } else {
        // no matching secondary vertex exists -> create new one
        particleId.setVertexSecondary(vertices.size());
        auto& vertex = vertices.emplace_back(
            static_cast<SimVertexBarcode>(particleId.vertexId()), pos4);
        vertex.outgoing.insert(particleId);
        ACTS_VERBOSE("created new secondary vertex " << pos4.transpose());
      }
    } else {
      auto& primaryVertex = vertices.front();
      primaryVertex.outgoing.insert(particleId);
    }

    // construct internal particle
    const auto pdg = static_cast<Acts::PdgParticle>(genParticle.id());
    const auto charge = genParticle.charge() * 1_e;
    const auto mass = genParticle.m0() * 1_GeV;
    SimParticleState particle(particleId, pdg, charge, mass);
    particle.setPosition4(pos4);
    // normalization/ units are not import for the direction
    particle.setDirection(genParticle.px(), genParticle.py(), genParticle.pz());
    particle.setAbsoluteMomentum(
        Acts::fastHypot(genParticle.px(), genParticle.py(), genParticle.pz()) *
        1_GeV);

    particles.push_back(SimParticle(particle, particle));
  }

  std::pair<SimVertexContainer, SimParticleContainer> out;
  out.first.insert(vertices.begin(), vertices.end());
  out.second.insert(particles.begin(), particles.end());

  m_impl->m_pythia8RndmEngine->clearRandomEngine();

  return out;
}

}  // namespace ActsExamples
