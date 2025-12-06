// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootVertexReader.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/Root/RootUtility.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

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
  m_inputChain = std::make_unique<TChain>(m_cfg.treeName.c_str());

  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing input filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_outputVertices.initialize(m_cfg.outputVertices);

  // Set the branches
  m_inputChain->SetBranchAddress("event_id", &m_eventId);
  m_inputChain->SetBranchAddress("process", &m_process.get());
  m_inputChain->SetBranchAddress("vx", &m_vx.get());
  m_inputChain->SetBranchAddress("vy", &m_vy.get());
  m_inputChain->SetBranchAddress("vz", &m_vz.get());
  m_inputChain->SetBranchAddress("vt", &m_vt.get());
  if (m_inputChain->GetBranch("incoming_particles") != nullptr) {
    m_hasCombinedIncoming = true;
    m_incomingParticles.allocate();
    m_inputChain->SetBranchAddress("incoming_particles",
                                   &m_incomingParticles.get());
  } else {
    m_hasCombinedIncoming = false;
    m_incomingParticlesVertexPrimary.allocate();
    m_incomingParticlesVertexSecondary.allocate();
    m_incomingParticlesParticle.allocate();
    m_incomingParticlesGeneration.allocate();
    m_incomingParticlesSubParticle.allocate();
    m_inputChain->SetBranchAddress("incoming_particles_vertex_primary",
                                   &m_incomingParticlesVertexPrimary.get());
    m_inputChain->SetBranchAddress("incoming_particles_vertex_secondary",
                                   &m_incomingParticlesVertexSecondary.get());
    m_inputChain->SetBranchAddress("incoming_particles_particle",
                                   &m_incomingParticlesParticle.get());
    m_inputChain->SetBranchAddress("incoming_particles_generation",
                                   &m_incomingParticlesGeneration.get());
    m_inputChain->SetBranchAddress("incoming_particles_sub_particle",
                                   &m_incomingParticlesSubParticle.get());
  }

  if (m_inputChain->GetBranch("outgoing_particles") != nullptr) {
    m_hasCombinedOutgoing = true;
    m_outgoingParticles.allocate();
    m_inputChain->SetBranchAddress("outgoing_particles",
                                   &m_outgoingParticles.get());
  } else {
    m_hasCombinedOutgoing = false;
    m_outgoingParticlesVertexPrimary.allocate();
    m_outgoingParticlesVertexSecondary.allocate();
    m_outgoingParticlesParticle.allocate();
    m_outgoingParticlesGeneration.allocate();
    m_outgoingParticlesSubParticle.allocate();
    m_inputChain->SetBranchAddress("outgoing_particles_vertex_primary",
                                   &m_outgoingParticlesVertexPrimary.get());
    m_inputChain->SetBranchAddress("outgoing_particles_vertex_secondary",
                                   &m_outgoingParticlesVertexSecondary.get());
    m_inputChain->SetBranchAddress("outgoing_particles_particle",
                                   &m_outgoingParticlesParticle.get());
    m_inputChain->SetBranchAddress("outgoing_particles_generation",
                                   &m_outgoingParticlesGeneration.get());
    m_inputChain->SetBranchAddress("outgoing_particles_sub_particle",
                                   &m_outgoingParticlesSubParticle.get());
  }
  m_inputChain->SetBranchAddress("vertex_primary", &m_vertexPrimary.get());
  m_inputChain->SetBranchAddress("vertex_secondary", &m_vertexSecondary.get());
  m_inputChain->SetBranchAddress("generation", &m_generation.get());

  auto path = m_cfg.filePath;

  // add file to the input chain
  m_inputChain->Add(path.c_str());
  ACTS_DEBUG("Adding File " << path << " to tree '" << m_cfg.treeName << "'.");

  m_events = m_inputChain->GetEntries();
  ACTS_DEBUG("The full chain has " << m_events << " entries.");

  // Sort the entry numbers of the events
  {
    // necessary to guarantee that m_inputChain->GetV1() is valid for the
    // entire range
    m_inputChain->SetEstimate(m_events + 1);

    m_entryNumbers.resize(m_events);
    m_inputChain->Draw("event_id", "", "goff");
    RootUtility::stableSort(m_inputChain->GetEntries(), m_inputChain->GetV1(),
                            m_entryNumbers.data(), false);
  }
}

std::pair<std::size_t, std::size_t> RootVertexReader::availableEvents() const {
  return {0u, m_events};
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

  unsigned int nVertices = m_vx->size();

  for (unsigned int i = 0; i < nVertices; i++) {
    SimVertex v;

    v.id = SimVertexBarcode()
               .withVertexPrimary((*m_vertexPrimary)[i])
               .withVertexSecondary((*m_vertexSecondary)[i])
               .withGeneration((*m_generation)[i]);

    v.process = static_cast<ActsFatras::ProcessType>((*m_process)[i]);
    v.position4 = Acts::Vector4((*m_vx)[i] * Acts::UnitConstants::mm,
                                (*m_vy)[i] * Acts::UnitConstants::mm,
                                (*m_vz)[i] * Acts::UnitConstants::mm,
                                (*m_vt)[i] * Acts::UnitConstants::mm);

    // incoming particles
    if (m_hasCombinedIncoming && m_incomingParticles.hasValue()) {
      for (auto& barcodeComponents : (*m_incomingParticles)[i]) {
        v.incoming.insert(SimBarcode().withData(barcodeComponents));
      }
    } else if (m_incomingParticlesVertexPrimary.hasValue()) {
      const auto& incomingPrimaries = (*m_incomingParticlesVertexPrimary)[i];
      const auto& incomingSecondaries =
          (*m_incomingParticlesVertexSecondary)[i];
      const auto& incomingParticles = (*m_incomingParticlesParticle)[i];
      const auto& incomingGenerations = (*m_incomingParticlesGeneration)[i];
      const auto& incomingSubParticles = (*m_incomingParticlesSubParticle)[i];
      for (std::size_t j = 0; j < incomingPrimaries.size(); ++j) {
        v.incoming.insert(
            SimBarcode()
                .withVertexPrimary(static_cast<SimBarcode::PrimaryVertexId>(
                    incomingPrimaries[j]))
                .withVertexSecondary(static_cast<SimBarcode::SecondaryVertexId>(
                    incomingSecondaries[j]))
                .withParticle(
                    static_cast<SimBarcode::ParticleId>(incomingParticles[j]))
                .withGeneration(static_cast<SimBarcode::GenerationId>(
                    incomingGenerations[j]))
                .withSubParticle(static_cast<SimBarcode::SubParticleId>(
                    incomingSubParticles[j])));
      }
    }

    // outgoing particles
    if (m_hasCombinedOutgoing && m_outgoingParticles.hasValue()) {
      for (auto& barcodeComponents : (*m_outgoingParticles)[i]) {
        v.outgoing.insert(SimBarcode().withData(barcodeComponents));
      }
    } else if (m_outgoingParticlesVertexPrimary.hasValue()) {
      const auto& outgoingPrimaries = (*m_outgoingParticlesVertexPrimary)[i];
      const auto& outgoingSecondaries =
          (*m_outgoingParticlesVertexSecondary)[i];
      const auto& outgoingParticles = (*m_outgoingParticlesParticle)[i];
      const auto& outgoingGenerations = (*m_outgoingParticlesGeneration)[i];
      const auto& outgoingSubParticles = (*m_outgoingParticlesSubParticle)[i];
      for (std::size_t j = 0; j < outgoingPrimaries.size(); ++j) {
        v.outgoing.insert(
            SimBarcode()
                .withVertexPrimary(static_cast<SimBarcode::PrimaryVertexId>(
                    outgoingPrimaries[j]))
                .withVertexSecondary(static_cast<SimBarcode::SecondaryVertexId>(
                    outgoingSecondaries[j]))
                .withParticle(
                    static_cast<SimBarcode::ParticleId>(outgoingParticles[j]))
                .withGeneration(static_cast<SimBarcode::GenerationId>(
                    outgoingGenerations[j]))
                .withSubParticle(static_cast<SimBarcode::SubParticleId>(
                    outgoingSubParticles[j])));
      }
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
