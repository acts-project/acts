// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootVertexReader.hpp"

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>

#include <TChain.h>
#include <TMathBase.h>

namespace ActsExamples {

RootVertexReader::RootVertexReader(const RootVertexReader::Config& config,
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

  m_outputVertices.initialize(m_cfg.outputVertices);

  // Set the branches
  m_inputChain->SetBranchAddress("event_id", &m_eventId);
  m_inputChain->SetBranchAddress("vertex_id", &m_vertexId);
  m_inputChain->SetBranchAddress("process", &m_process);
  m_inputChain->SetBranchAddress("vx", &m_vx);
  m_inputChain->SetBranchAddress("vy", &m_vy);
  m_inputChain->SetBranchAddress("vz", &m_vz);
  m_inputChain->SetBranchAddress("vt", &m_vt);
  m_inputChain->SetBranchAddress("outgoing_particles", &m_outgoingParticles);
  m_inputChain->SetBranchAddress("vertex_primary", &m_vertexPrimary);
  m_inputChain->SetBranchAddress("vertex_secondary", &m_vertexSecondary);
  m_inputChain->SetBranchAddress("generation", &m_generation);

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
    // Sort to get the entry numbers of the ordered events
    TMath::Sort(m_inputChain->GetEntries(), m_inputChain->GetV1(),
                m_entryNumbers.data(), false);
  }
}

std::pair<std::size_t, std::size_t> RootVertexReader::availableEvents() const {
  return {0u, m_events};
}

RootVertexReader::~RootVertexReader() {
  delete m_vertexId;
  delete m_process;
  delete m_vx;
  delete m_vy;
  delete m_vz;
  delete m_vt;
  delete m_outgoingParticles;
  delete m_vertexPrimary;
  delete m_vertexSecondary;
  delete m_generation;
}

ProcessCode RootVertexReader::read(const AlgorithmContext& context) {
  ACTS_DEBUG("Trying to read recorded vertices.");

  if (m_inputChain == nullptr || context.eventNumber >= m_events) {
    return ProcessCode::SUCCESS;
  }

  // lock the mutex
  std::lock_guard<std::mutex> lock(m_read_mutex);
  // now read

  // The vertex collection to be filled
  SimVertexContainer vertices;

  // Read the correct entry
  auto entry = m_entryNumbers.at(context.eventNumber);
  m_inputChain->GetEntry(entry);
  ACTS_DEBUG("Reading event: " << context.eventNumber
                               << " stored as entry: " << entry);

  unsigned int nVertices = m_vertexId->size();

  for (unsigned int i = 0; i < nVertices; i++) {
    SimVertex v;

    v.id = (*m_vertexId)[i];
    v.process = static_cast<ActsFatras::ProcessType>((*m_process)[i]);
    v.position4 = Acts::Vector4((*m_vx)[i] * Acts::UnitConstants::mm,
                                (*m_vy)[i] * Acts::UnitConstants::mm,
                                (*m_vz)[i] * Acts::UnitConstants::mm,
                                (*m_vt)[i] * Acts::UnitConstants::mm);

    // TODO ingoing particles

    for (auto& id : (*m_outgoingParticles)[i]) {
      v.outgoing.insert(static_cast<std::uint64_t>(id));
    }

    vertices.insert(v);
  }

  ACTS_DEBUG("Read " << vertices.size() << " vertices for event "
                     << context.eventNumber);

  // Write the collections to the EventStore
  m_outputVertices(context, std::move(vertices));

  // Return success flag
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
