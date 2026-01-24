// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootVertexWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <cstdint>
#include <ios>
#include <ostream>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

namespace ActsExamples {

RootVertexWriter::RootVertexWriter(const RootVertexWriter::Config& cfg,
                                   Acts::Logging::Level lvl)
    : WriterT(cfg.inputVertices, "RootVertexWriter", lvl), m_cfg(cfg) {
  // inputParticles is already checked by base constructor
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing file path");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // open root file and create the tree
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  // setup the branches
  m_outputTree->Branch("event_id", &m_eventId);
  m_outputTree->Branch("process", &m_process);
  m_outputTree->Branch("vx", &m_vx);
  m_outputTree->Branch("vy", &m_vy);
  m_outputTree->Branch("vz", &m_vz);
  m_outputTree->Branch("vt", &m_vt);
  m_outputTree->Branch("incoming_particles_vertex_primary",
                       &m_incomingParticlesVertexPrimary);
  m_outputTree->Branch("incoming_particles_vertex_secondary",
                       &m_incomingParticlesVertexSecondary);
  m_outputTree->Branch("incoming_particles_particle",
                       &m_incomingParticlesParticle);
  m_outputTree->Branch("incoming_particles_generation",
                       &m_incomingParticlesGeneration);
  m_outputTree->Branch("incoming_particles_sub_particle",
                       &m_incomingParticlesSubParticle);
  m_outputTree->Branch("outgoing_particles_vertex_primary",
                       &m_outgoingParticlesVertexPrimary);
  m_outputTree->Branch("outgoing_particles_vertex_secondary",
                       &m_outgoingParticlesVertexSecondary);
  m_outputTree->Branch("outgoing_particles_particle",
                       &m_outgoingParticlesParticle);
  m_outputTree->Branch("outgoing_particles_generation",
                       &m_outgoingParticlesGeneration);
  m_outputTree->Branch("outgoing_particles_sub_particle",
                       &m_outgoingParticlesSubParticle);
  m_outputTree->Branch("vertex_primary", &m_vertexPrimary);
  m_outputTree->Branch("vertex_secondary", &m_vertexSecondary);
  m_outputTree->Branch("generation", &m_generation);
}

RootVertexWriter::~RootVertexWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode RootVertexWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  ACTS_INFO("Wrote vertices to tree '" << m_cfg.treeName << "' in '"
                                       << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}

ProcessCode RootVertexWriter::writeT(const AlgorithmContext& ctx,
                                     const SimVertexContainer& vertices) {
  // ensure exclusive access to tree/file while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  m_eventId = ctx.eventNumber;
  for (const auto& vertex : vertices) {
    m_process.push_back(static_cast<std::uint32_t>(vertex.process));
    // position
    m_vx.push_back(Acts::clampValue<float>(vertex.position4.x() /
                                           Acts::UnitConstants::mm));
    m_vy.push_back(Acts::clampValue<float>(vertex.position4.y() /
                                           Acts::UnitConstants::mm));
    m_vz.push_back(Acts::clampValue<float>(vertex.position4.z() /
                                           Acts::UnitConstants::mm));
    m_vt.push_back(Acts::clampValue<float>(vertex.position4.w() /
                                           Acts::UnitConstants::mm));

    // incoming particles
    std::vector<std::uint32_t> incoming_vertex_primary;
    std::vector<std::uint32_t> incoming_vertex_secondary;
    std::vector<std::uint32_t> incoming_particle_component;
    std::vector<std::uint32_t> incoming_generation;
    std::vector<std::uint32_t> incoming_sub_particle;
    for (const auto& particle : vertex.incoming) {
      incoming_vertex_primary.push_back(particle.vertexPrimary());
      incoming_vertex_secondary.push_back(particle.vertexSecondary());
      incoming_particle_component.push_back(particle.particle());
      incoming_generation.push_back(particle.generation());
      incoming_sub_particle.push_back(particle.subParticle());
    }
    m_incomingParticlesVertexPrimary.push_back(
        std::move(incoming_vertex_primary));
    m_incomingParticlesVertexSecondary.push_back(
        std::move(incoming_vertex_secondary));
    m_incomingParticlesParticle.push_back(
        std::move(incoming_particle_component));
    m_incomingParticlesGeneration.push_back(std::move(incoming_generation));
    m_incomingParticlesSubParticle.push_back(std::move(incoming_sub_particle));

    // outgoing particles
    std::vector<std::uint32_t> outgoing_vertex_primary;
    std::vector<std::uint32_t> outgoing_vertex_secondary;
    std::vector<std::uint32_t> outgoing_particle_component;
    std::vector<std::uint32_t> outgoing_generation;
    std::vector<std::uint32_t> outgoing_sub_particle;
    for (const auto& particle : vertex.outgoing) {
      outgoing_vertex_primary.push_back(particle.vertexPrimary());
      outgoing_vertex_secondary.push_back(particle.vertexSecondary());
      outgoing_particle_component.push_back(particle.particle());
      outgoing_generation.push_back(particle.generation());
      outgoing_sub_particle.push_back(particle.subParticle());
    }
    m_outgoingParticlesVertexPrimary.push_back(
        std::move(outgoing_vertex_primary));
    m_outgoingParticlesVertexSecondary.push_back(
        std::move(outgoing_vertex_secondary));
    m_outgoingParticlesParticle.push_back(
        std::move(outgoing_particle_component));
    m_outgoingParticlesGeneration.push_back(std::move(outgoing_generation));
    m_outgoingParticlesSubParticle.push_back(std::move(outgoing_sub_particle));

    // decoded barcode components
    m_vertexPrimary.push_back(vertex.vertexId().vertexPrimary());
    m_vertexSecondary.push_back(vertex.vertexId().vertexSecondary());
    m_generation.push_back(vertex.vertexId().generation());
  }

  m_outputTree->Fill();

  m_process.clear();
  m_vx.clear();
  m_vy.clear();
  m_vz.clear();
  m_vt.clear();
  m_incomingParticlesVertexPrimary.clear();
  m_incomingParticlesVertexSecondary.clear();
  m_incomingParticlesParticle.clear();
  m_incomingParticlesGeneration.clear();
  m_incomingParticlesSubParticle.clear();
  m_outgoingParticlesVertexPrimary.clear();
  m_outgoingParticlesVertexSecondary.clear();
  m_outgoingParticlesParticle.clear();
  m_outgoingParticlesGeneration.clear();
  m_outgoingParticlesSubParticle.clear();
  m_vertexPrimary.clear();
  m_vertexSecondary.clear();
  m_generation.clear();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
