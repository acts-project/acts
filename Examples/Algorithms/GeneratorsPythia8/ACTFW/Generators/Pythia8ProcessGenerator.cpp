// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Generators/Pythia8ProcessGenerator.hpp"

#include <Pythia8/Pythia.h>
#include <algorithm>
#include <iterator>
#include <random>

namespace {
struct FrameworkRndmEngine : public Pythia8::RndmEngine {
  FW::RandomEngine& rng;

  FrameworkRndmEngine(FW::RandomEngine& rng_) : rng(rng_) {}
  double flat() {
    return std::uniform_real_distribution<double>(0.0, 1.0)(rng);
  }
};
}  // namespace

std::function<std::vector<FW::SimVertex>(FW::RandomEngine&)>
FW::Pythia8Generator::makeFunction(const FW::Pythia8Generator::Config& cfg,
                                   Acts::Logging::Level lvl) {
  auto gen = std::make_shared<Pythia8Generator>(cfg, lvl);
  return [=](RandomEngine& rng) { return (*gen)(rng); };
}

FW::Pythia8Generator::Pythia8Generator(const FW::Pythia8Generator::Config& cfg,
                                       Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("Pythia8Generator", lvl)),
      m_pythia8(std::make_unique<Pythia8::Pythia>("", false)) {
  // disable all output by default but allow reenable via config
  m_pythia8->settings.flag("Print:quiet", true);
  for (const auto& str : m_cfg.settings) {
    ACTS_VERBOSE("use Pythia8 setting '" << str << "'");
    m_pythia8->readString(str.c_str());
  }
  m_pythia8->settings.mode("Beams:idA", m_cfg.pdgBeam0);
  m_pythia8->settings.mode("Beams:idB", m_cfg.pdgBeam1);
  m_pythia8->settings.mode("Beams:frameType", 1);
  m_pythia8->settings.parm("Beams:eCM",
                           m_cfg.cmsEnergy / Acts::UnitConstants::GeV);
  m_pythia8->init();
}

// needed to allow unique_ptr of forward-declared Pythia class
FW::Pythia8Generator::~Pythia8Generator() {}

std::vector<FW::SimVertex> FW::Pythia8Generator::operator()(
    FW::RandomEngine& rng) {
  using namespace Acts::UnitLiterals;

  // first process vertex is the primary one at origin with time=0
  std::vector<SimVertex> vertices = {
      SimVertex(SimVertex::Vector4::Zero()),
  };

  // pythia8 is not thread safe and generation needs to be protected
  std::lock_guard<std::mutex> lock(m_pythia8Mutex);
  // use per-thread random engine also in pythia
  FrameworkRndmEngine rndmEngine(rng);
  m_pythia8->rndm.rndmEnginePtr(&rndmEngine);
  m_pythia8->next();

  // convert generated final state particles into internal format
  for (int ip = 0; ip < m_pythia8->event.size(); ++ip) {
    const auto& genParticle = m_pythia8->event[ip];

    // ignore beam particles
    if (genParticle.statusHepMC() == 4) {
      continue;
    }
    // only interested in final, visible particles
    if (not genParticle.isFinal()) {
      continue;
    }
    if (not genParticle.isVisible()) {
      continue;
    }

    // production vertex. Pythia8 time uses units mm/c, and we use c=1
    SimVertex::Vector4 pos4(
        genParticle.xProd() * 1_mm, genParticle.yProd() * 1_mm,
        genParticle.zProd() * 1_mm, genParticle.tProd() * 1_mm);
    // identify vertex
    std::vector<SimVertex>::iterator vertex;
    if (not genParticle.hasVertex()) {
      // w/o defined vertex, must belong to the first (primary) process vertex
      vertex = vertices.begin();
    } else {
      // either add to existing secondary vertex if exists or create new one
      // TODO can we do this w/o the manual search and position check?
      vertex = std::find_if(
          vertices.begin(), vertices.end(),
          [=](const SimVertex& v) { return (v.position4 == pos4); });
      if (vertex == vertices.end()) {
        // no matching secondary vertex exists -> create new one
        vertices.emplace_back(pos4);
        vertex = std::prev(vertices.end());

        ACTS_VERBOSE("created new secondary vertex " << pos4.transpose());
      }
    }

    // Ensure particle identifier components are defaulted to zero.
    ActsFatras::Barcode particleId(0u);
    // first vertex w/ distance=0 contains all direct particles
    particleId.setVertexSecondary(std::distance(vertices.begin(), vertex));
    // ensure particle identifier component is non-zero
    particleId.setParticle(1u + vertex->outgoing.size());
    // reuse PDG id from generator
    const auto pdg = static_cast<Acts::PdgParticle>(genParticle.id());

    // construct internal particle
    ActsFatras::Particle particle(particleId, pdg, genParticle.charge() * 1_e,
                                  genParticle.m0() * 1_GeV);
    particle.setPosition4(pos4);
    // normalization/ units are not import for the direction
    particle.setDirection(genParticle.px(), genParticle.py(), genParticle.pz());
    particle.setAbsMomentum(
        std::hypot(genParticle.px(), genParticle.py(), genParticle.pz()) *
        1_GeV);

    vertex->outgoing.push_back(std::move(particle));
  }
  return vertices;
}
