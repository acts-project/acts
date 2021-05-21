// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootSeedingPerformanceWriter.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <ios>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

ActsExamples::RootSeedingPerformanceWriter::RootSeedingPerformanceWriter(
    const ActsExamples::RootSeedingPerformanceWriter::Config& cfg,
    Acts::Logging::Level level)
    : WriterT(cfg.collection, "RootSeedingPerformanceWriter", level),
      m_cfg(cfg),
      m_outputFile(cfg.rootFile) {
  // An input collection name and tree name must be specified
  if (m_cfg.collection.empty()) {
    throw std::invalid_argument("Missing input collection");
  } else if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
    }
  }
  m_outputFile->cd();

  m_outputTree = new TTree(m_cfg.treeName.c_str(),
                           "TTree from RootSeedingPerformanceWriter");
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  // Set the branches
  m_outputTree->Branch("event_nr", &m_eventNr);
  m_outputTree->Branch("seeds", &m_seeds);
  for (unsigned int isp = 0; isp < 3; ++isp) {
    std::string spX = std::string("sp") + std::to_string(isp);
    // m_outputTree->Branch((spX+std::string("volume_id")).c_str(),
    // &m_volumeID[isp]);
    // m_outputTree->Branch((spX+std::string("layer_id")).c_str(),
    // &m_layerID[isp]);
    // m_outputTree->Branch((spX+std::string("sensitive_id")).c_str(),
    // &m_sensitiveID[isp]);
    m_outputTree->Branch((spX + std::string("g_x")).c_str(), &m_x[isp]);
    m_outputTree->Branch((spX + std::string("g_y")).c_str(), &m_y[isp]);
    m_outputTree->Branch((spX + std::string("g_z")).c_str(), &m_z[isp]);
  }

  // m_outputTree->Branch("seed_d0", &m_d0Rec);
  // m_outputTree->Branch("seed_z0", &m_z0Rec);
  // m_outputTree->Branch("seed_phi", &m_phiRec);
  // m_outputTree->Branch("seed_theta", &m_thetaRec);
  // m_outputTree->Branch("seed_qop", &m_qopRec);
}

ActsExamples::RootSeedingPerformanceWriter::~RootSeedingPerformanceWriter() {
  /// Close the file if it's yours
  if (m_cfg.rootFile == nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootSeedingPerformanceWriter::endRun() {
  // Write the tree
  m_outputFile->cd();
  m_outputTree->Write();
  ACTS_VERBOSE("Wrote particles to tree '" << m_cfg.treeName << "' in '"
                                           << m_cfg.filePath << "'");
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootSeedingPerformanceWriter::writeT(
    const AlgorithmContext& context, const SimSeedContainer& seeds) {
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = context.eventNumber;
  m_seeds = seeds.size();

  // Reserve the vectors
  for (unsigned int isp = 0; isp < 3; ++isp) {
    m_x[isp].reserve(seeds.size());
    m_y[isp].reserve(seeds.size());
    m_z[isp].reserve(seeds.size());
  }

  // Fill the vectors
  for (const auto& seed : seeds) {
    for (unsigned int isp = 0; isp < 3; ++isp) {
      const auto& sp = seed.sp()[isp];
      m_x[isp].push_back(sp->x());
      m_y[isp].push_back(sp->y());
      m_z[isp].push_back(sp->z());
    }
  }
  m_outputTree->Fill();

  // Clear the vectors
  for (unsigned int isp = 0; isp < 3; ++isp) {
    m_x[isp].clear();
    m_y[isp].clear();
    m_z[isp].clear();
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
