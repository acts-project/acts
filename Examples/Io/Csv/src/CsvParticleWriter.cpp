// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Csv/CsvParticleWriter.hpp"

#include <Acts/Utilities/Units.hpp>
#include <dfe/dfe_io_dsv.hpp>
#include <map>
#include <stdexcept>

#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "TrackMlData.hpp"

FW::CsvParticleWriter::CsvParticleWriter(
    const FW::CsvParticleWriter::Config& cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.inputParticles, "CsvParticleWriter", lvl), m_cfg(cfg) {
  // inputParticles is already checked by base constructor
  if (m_cfg.outputStem.empty()) {
    throw std::invalid_argument("Missing ouput filename stem");
  }
}

FW::ProcessCode FW::CsvParticleWriter::writeT(
    const FW::AlgorithmContext& ctx, const SimParticleContainer& particles) {
  auto pathParticles = perEventFilepath(
      m_cfg.outputDir, m_cfg.outputStem + ".csv", ctx.eventNumber);
  dfe::NamedTupleCsvWriter<ParticleData> writer(pathParticles,
                                                m_cfg.outputPrecision);

  ParticleData data;
  for (const auto& particle : particles) {
    data.particle_id = particle.particleId().value();
    data.particle_type = particle.pdg();
    data.process = static_cast<decltype(data.process)>(particle.process());
    data.vx = particle.position().x() / Acts::UnitConstants::mm;
    data.vy = particle.position().y() / Acts::UnitConstants::mm;
    data.vz = particle.position().z() / Acts::UnitConstants::mm;
    data.vt = particle.time() / Acts::UnitConstants::ns;
    const auto p = particle.absMomentum() / Acts::UnitConstants::GeV;
    data.px = p * particle.unitDirection().x();
    data.py = p * particle.unitDirection().y();
    data.pz = p * particle.unitDirection().z();
    data.m = particle.mass() / Acts::UnitConstants::GeV;
    data.q = particle.charge() / Acts::UnitConstants::e;
    writer.append(data);
  }

  return ProcessCode::SUCCESS;
}
