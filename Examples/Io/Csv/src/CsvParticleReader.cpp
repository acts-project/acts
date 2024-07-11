// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <map>
#include <iostream>

#include <dfe/dfe_io_dsv.hpp>

#include "CsvOutputData.hpp"

ActsExamples::CsvParticleReader::CsvParticleReader(
    const ActsExamples::CsvParticleReader::Config& config,
    Acts::Logging::Level level)
    : m_cfg(config),
      m_eventsRange(
          determineEventFilesRange(m_cfg.inputDir, m_cfg.inputStem + ".csv")),
      m_logger(Acts::getDefaultLogger("CsvParticleReader", level)) {
  if (m_cfg.inputStem.empty()) {
    throw std::invalid_argument("Missing input filename stem");
  }
  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument("Missing output collection");
  }

  m_outputParticles.initialize(m_cfg.outputParticles);
  m_outputVertices.initialize(m_cfg.outputVertices);

}

std::string ActsExamples::CsvParticleReader::CsvParticleReader::name() const {
  return "CsvParticleReader";
}

std::pair<std::size_t, std::size_t>
ActsExamples::CsvParticleReader::availableEvents() const {
  return m_eventsRange;
}

ActsExamples::ProcessCode ActsExamples::CsvParticleReader::read(
    const ActsExamples::AlgorithmContext& ctx) {
  SimParticleContainer::sequence_type unordered;
  SimVertexContainer::sequence_type vertices;

  std::map<ActsFatras::Barcode, SimVertex> foundVertices;

  auto path = perEventFilepath(m_cfg.inputDir, m_cfg.inputStem + ".csv",
                               ctx.eventNumber);
  // vt and m are an optional columns
  dfe::NamedTupleCsvReader<ParticleData> reader(path, {"vt", "m"});
  ParticleData data;
  
  while (reader.read(data)) {
    ActsFatras::Barcode particle_id(data.particle_id);
    ActsFatras::Barcode vertex_id = particle_id.vertexId();
  
    ActsFatras::Particle particle(particle_id,
                                  Acts::PdgParticle{data.particle_type},
                                  data.q * Acts::UnitConstants::e,
                                  data.m * Acts::UnitConstants::GeV);
    particle.setProcess(static_cast<ActsFatras::ProcessType>(data.process));
    particle.setPosition4(
        data.vx * Acts::UnitConstants::mm, data.vy * Acts::UnitConstants::mm,
        data.vz * Acts::UnitConstants::mm, data.vt * Acts::UnitConstants::mm);
    // Only used for direction; normalization/units do not matter
    particle.setDirection(data.px, data.py, data.pz);
    particle.setAbsoluteMomentum(std::hypot(data.px, data.py, data.pz) *
                                 Acts::UnitConstants::GeV);
    unordered.push_back(std::move(particle));

    if (foundVertices.find(vertex_id) != foundVertices.end()) {
      foundVertices[vertex_id].outgoing.insert(particle_id);
    }
    else
    {
      foundVertices[vertex_id] = vertices.emplace_back(vertex_id, SimVertex::Vector4(data.vx * Acts::UnitConstants::mm,
                                                      data.vy * Acts::UnitConstants::mm,
                                                      data.vz * Acts::UnitConstants::mm,
                                                      data.vt * Acts::UnitConstants::mm));
    }
  }

  // Write ordered particles container to the EventStore
  SimParticleContainer particles;
  particles.insert(unordered.begin(), unordered.end());
  SimVertexContainer vert;
  vert.insert(vertices.begin(), vertices.end());
  m_outputParticles(ctx, std::move(particles));
  m_outputVertices(ctx, std::move(vert));
  return ProcessCode::SUCCESS;
}
