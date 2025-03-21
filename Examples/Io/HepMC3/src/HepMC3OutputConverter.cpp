// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3OutputConverter.hpp"

#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsExamples {

HepMC3OutputConverter::HepMC3OutputConverter(const Config& config,
                                             Acts::Logging::Level level)
    : IAlgorithm{"HepMC3OutputConverter", level}, m_cfg{config} {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles collection");
  }

  if (m_cfg.inputVertices.empty()) {
    throw std::invalid_argument("Missing input vertices collection");
  }

  if (m_cfg.outputEvents.empty()) {
    throw std::invalid_argument("Missing output event collection");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputVertices.initialize(m_cfg.inputVertices);
  m_outputEvents.initialize(m_cfg.outputEvents);
}

ProcessCode HepMC3OutputConverter::execute(const AlgorithmContext& ctx) const {
  std::vector<HepMC3::GenEvent> genEvents;
  auto& genEvent =
      genEvents.emplace_back(HepMC3::Units::GEV, HepMC3::Units::MM);

  const auto& inputParticles = m_inputParticles(ctx);
  const auto& inputVertices = m_inputVertices(ctx);
  ACTS_DEBUG("Converting " << inputParticles.size() << " and "
                           << inputVertices.size() << " vertices to HepMC3");

  std::unordered_map<SimBarcode, std::shared_ptr<HepMC3::GenParticle>>
      barcodeMap;
  for (const auto& inParticle : inputParticles) {
    const Vector4 momentum =
        inParticle.fourMomentum() / 1_GeV;  // Let's ensure we're using GeV
    const HepMC3::FourVector vec(momentum[0], momentum[1], momentum[2],
                                 momentum[3]);
    auto hepmc3Particle =
        std::make_shared<HepMC3::GenParticle>(vec, inParticle.pdg());
    hepmc3Particle->set_generated_mass(inParticle.mass());
    genEvent.add_particle(hepmc3Particle);
    barcodeMap.insert({inParticle.particleId(), hepmc3Particle});
  }

  for (const auto& inVertex : inputVertices) {
    const Vector4 position =
        inVertex.position4 / 1_mm;  // Let's ensure we're using mm
    const HepMC3::FourVector vec(position[0], position[1], position[2],
                                 position[3]);
    auto hepmc3Vertex = std::make_shared<HepMC3::GenVertex>(vec);
    genEvent.add_vertex(hepmc3Vertex);

    for (const auto& particleId : inVertex.incoming) {
      auto it = barcodeMap.find(particleId);
      if (it == barcodeMap.end()) {
        ACTS_ERROR("Particle with barcode " << particleId << " not found");
        continue;
      }
      hepmc3Vertex->add_particle_in(it->second);
    }

    for (const auto& particleId : inVertex.outgoing) {
      auto it = barcodeMap.find(particleId);
      if (it == barcodeMap.end()) {
        ACTS_ERROR("Particle with barcode " << particleId << " not found");
        continue;
      }
      hepmc3Vertex->add_particle_out(it->second);
    }
  }

  m_outputEvents(ctx, std::move(genEvents));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
