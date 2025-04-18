// This file is part of the Acts project.
//
// Copyright (C) 2017-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootParticleReader.hpp"

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/Root/RootUtility.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>

#include <TChain.h>

namespace ActsExamples {

RootParticleReader::RootParticleReader(const RootParticleReader::Config& config,
                                       Acts::Logging::Level level)
    : IReader(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  m_inputChain = new TChain(m_cfg.treeName.c_str());

  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing input filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_outputParticles.initialize(m_cfg.outputParticles);

  // Set the branches
  m_inputChain->SetBranchAddress("event_id", &m_eventId);
  m_inputChain->SetBranchAddress("particle_id", &m_particleId);
  m_inputChain->SetBranchAddress("particle_type", &m_particleType);
  m_inputChain->SetBranchAddress("process", &m_process);
  m_inputChain->SetBranchAddress("vx", &m_vx);
  m_inputChain->SetBranchAddress("vy", &m_vy);
  m_inputChain->SetBranchAddress("vz", &m_vz);
  m_inputChain->SetBranchAddress("vt", &m_vt);
  m_inputChain->SetBranchAddress("p", &m_p);
  m_inputChain->SetBranchAddress("px", &m_px);
  m_inputChain->SetBranchAddress("py", &m_py);
  m_inputChain->SetBranchAddress("pz", &m_pz);
  m_inputChain->SetBranchAddress("m", &m_m);
  m_inputChain->SetBranchAddress("q", &m_q);
  m_inputChain->SetBranchAddress("eta", &m_eta);
  m_inputChain->SetBranchAddress("phi", &m_phi);
  m_inputChain->SetBranchAddress("pt", &m_pt);
  m_inputChain->SetBranchAddress("vertex_primary", &m_vertexPrimary);
  m_inputChain->SetBranchAddress("vertex_secondary", &m_vertexSecondary);
  m_inputChain->SetBranchAddress("particle", &m_particle);
  m_inputChain->SetBranchAddress("generation", &m_generation);
  m_inputChain->SetBranchAddress("sub_particle", &m_subParticle);

  auto path = m_cfg.filePath;

  // add file to the input chain
  m_inputChain->Add(path.c_str());
  ACTS_DEBUG("Adding File " << path << " to tree '" << m_cfg.treeName << "'.");

  m_events = m_inputChain->GetEntries();
  ACTS_DEBUG("The full chain has " << m_events << " entries.");

  // Sort the entry numbers of the events
  {
    m_entryNumbers.resize(m_events);
    m_inputChain->Draw("event_id", "", "goff");
    RootUtility::stableSort(m_inputChain->GetEntries(), m_inputChain->GetV1(),
                            m_entryNumbers.data(), false);
  }
}

std::pair<std::size_t, std::size_t> RootParticleReader::availableEvents()
    const {
  return {0u, m_events};
}

RootParticleReader::~RootParticleReader() {
  delete m_particleId;
  delete m_particleType;
  delete m_process;
  delete m_vx;
  delete m_vy;
  delete m_vz;
  delete m_vt;
  delete m_p;
  delete m_px;
  delete m_py;
  delete m_pz;
  delete m_m;
  delete m_q;
  delete m_eta;
  delete m_phi;
  delete m_pt;
  delete m_vertexPrimary;
  delete m_vertexSecondary;
  delete m_particle;
  delete m_generation;
  delete m_subParticle;
}

ProcessCode RootParticleReader::read(const AlgorithmContext& context) {
  ACTS_DEBUG("Trying to read recorded particles.");

  if (m_inputChain == nullptr || context.eventNumber >= m_events) {
    return ProcessCode::SUCCESS;
  }

  // lock the mutex
  std::lock_guard<std::mutex> lock(m_read_mutex);
  // now read

  // The particle collection to be filled
  SimParticleContainer particles;

  // Read the correct entry
  auto entry = m_entryNumbers.at(context.eventNumber);
  m_inputChain->GetEntry(entry);
  ACTS_DEBUG("Reading event: " << context.eventNumber
                               << " stored as entry: " << entry);

  unsigned int nParticles = m_particleId->size();

  for (unsigned int i = 0; i < nParticles; i++) {
    SimParticle p;

    p.setProcess(static_cast<ActsFatras::ProcessType>((*m_process)[i]));
    p.setPdg(static_cast<Acts::PdgParticle>((*m_particleType)[i]));
    p.setCharge((*m_q)[i] * Acts::UnitConstants::e);
    p.setMass((*m_m)[i] * Acts::UnitConstants::GeV);
    p.setParticleId((*m_particleId)[i]);
    p.setPosition4((*m_vx)[i] * Acts::UnitConstants::mm,
                   (*m_vy)[i] * Acts::UnitConstants::mm,
                   (*m_vz)[i] * Acts::UnitConstants::mm,
                   (*m_vt)[i] * Acts::UnitConstants::mm);
    // NOTE: direction is normalized inside `setDirection`
    p.setDirection((*m_px)[i], (*m_py)[i], (*m_pz)[i]);
    p.setAbsoluteMomentum((*m_p)[i] * Acts::UnitConstants::GeV);

    particles.insert(p);
  }

  ACTS_DEBUG("Read " << particles.size() << " particles for event "
                     << context.eventNumber);

  // Write the collections to the EventStore
  m_outputParticles(context, std::move(particles));

  // Return success flag
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
