// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootParticleReader.hpp"

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/Root/RootUtility.hpp"
#include "ActsFatras/EventData/GenerationProcess.hpp"
#include "ActsFatras/EventData/SimulationOutcome.hpp"

#include <iostream>
#include <stdexcept>

#include <TChain.h>

namespace ActsExamples {

RootParticleReader::RootParticleReader(const RootParticleReader::Config& config,
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

  m_outputParticles.initialize(m_cfg.outputParticles);

  auto path = m_cfg.filePath;

  // add file to the input chain. This has to happen before the branches are
  // set up, otherwise the branch layout of the file is not known yet and
  // missing branches cannot be detected.
  m_inputChain->Add(path.c_str());
  ACTS_DEBUG("Adding File " << path << " to tree '" << m_cfg.treeName << "'.");

  // Without a loaded tree the branch layout is unknown and cannot be
  // validated. Reading is a no-op in that case since the chain has no entries.
  const bool treeLoaded = m_inputChain->LoadTree(0) >= 0;

  // Set the branches. A required branch that is absent would otherwise leave
  // the branch pointer null and be dereferenced during reading.
  auto setBranch = [&](const char* name, auto address) {
    if (treeLoaded && m_inputChain->FindBranch(name) == nullptr) {
      throw std::runtime_error("Branch '" + std::string(name) +
                               "' not found in input file '" + path + "'");
    }
    m_inputChain->SetBranchAddress(name, address);
  };

  // Branches which are not present in files written before the HF truth
  // information was introduced. Those fall back to the SimParticle defaults.
  auto setOptionalBranch = [&](const char* name, auto address) {
    if (treeLoaded && m_inputChain->FindBranch(name) == nullptr) {
      ACTS_WARNING("Branch '" << name << "' not found in input file '" << path
                              << "', falling back to default values");
      return;
    }
    m_inputChain->SetBranchAddress(name, address);
  };

  setBranch("event_id", &m_eventId);
  setBranch("particle_hash", &m_particleHash.get());
  setBranch("particle_type", &m_particleType.get());
  setBranch("process", &m_process.get());
  setBranch("vx", &m_vx.get());
  setBranch("vy", &m_vy.get());
  setBranch("vz", &m_vz.get());
  setBranch("vt", &m_vt.get());
  setBranch("p", &m_p.get());
  setBranch("px", &m_px.get());
  setBranch("py", &m_py.get());
  setBranch("pz", &m_pz.get());
  setBranch("m", &m_m.get());
  setBranch("q", &m_q.get());
  setBranch("eta", &m_eta.get());
  setBranch("phi", &m_phi.get());
  setBranch("pt", &m_pt.get());
  setBranch("vertex_primary", &m_vertexPrimary.get());
  setBranch("vertex_secondary", &m_vertexSecondary.get());
  setBranch("particle", &m_particle.get());
  setBranch("generation", &m_generation.get());
  setBranch("sub_particle", &m_subParticle.get());
  setOptionalBranch("orig_part_idx", &m_origParticleIdx.get());
  setOptionalBranch("hf_origin", &m_hfOrigin.get());

  setBranch("e_loss", &m_eLoss.get());
  setBranch("total_x0", &m_pathInX0.get());
  setBranch("total_l0", &m_pathInL0.get());
  setBranch("number_of_hits", &m_numberOfHits.get());
  setBranch("outcome", &m_outcome.get());

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

std::pair<std::size_t, std::size_t> RootParticleReader::availableEvents()
    const {
  return {0u, m_events};
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

  unsigned int nParticles = m_particleType->size();

  for (unsigned int i = 0; i < nParticles; i++) {
    SimParticle p;

    p.setProcess(
        static_cast<ActsFatras::GenerationProcess>((*m_process).at(i)));
    p.setPdg(static_cast<Acts::PdgParticle>((*m_particleType).at(i)));
    p.setCharge((*m_q).at(i) * Acts::UnitConstants::e);
    p.setMass((*m_m).at(i) * Acts::UnitConstants::GeV);
    p.setParticleId(SimBarcode()
                        .withVertexPrimary((*m_vertexPrimary).at(i))
                        .withVertexSecondary((*m_vertexSecondary).at(i))
                        .withParticle((*m_particle).at(i))
                        .withGeneration((*m_generation).at(i))
                        .withSubParticle((*m_subParticle).at(i)));

    if (m_origParticleIdx.hasValue()) {
      p.setOrigParticleIdx((*m_origParticleIdx).at(i));
    }
    if (m_hfOrigin.hasValue()) {
      p.setHeavyFlavourOrigin(
          static_cast<HeavyFlavourOrigin>((*m_hfOrigin).at(i)));
    }

    SimParticleState& initialState = p.initialState();

    initialState.setPosition4((*m_vx).at(i) * Acts::UnitConstants::mm,
                              (*m_vy).at(i) * Acts::UnitConstants::mm,
                              (*m_vz).at(i) * Acts::UnitConstants::mm,
                              (*m_vt).at(i) * Acts::UnitConstants::mm);
    // NOTE: direction is normalized inside `setDirection`
    initialState.setDirection((*m_px).at(i), (*m_py).at(i), (*m_pz).at(i));
    initialState.setAbsoluteMomentum((*m_p).at(i) * Acts::UnitConstants::GeV);

    SimParticleState& finalState = p.finalState();

    // TODO eloss cannot be read since we need the final momentum
    finalState.setMaterialPassed((*m_pathInX0).at(i) * Acts::UnitConstants::mm,
                                 (*m_pathInL0).at(i) * Acts::UnitConstants::mm);
    finalState.setNumberOfHits((*m_numberOfHits).at(i));
    finalState.setOutcome(
        static_cast<ActsFatras::SimulationOutcome>((*m_outcome).at(i)));

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
