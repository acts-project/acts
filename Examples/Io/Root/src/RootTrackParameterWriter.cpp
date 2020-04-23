// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Root/RootTrackParameterWriter.hpp"

#include <Acts/Utilities/Helpers.hpp>
#include <TFile.h>
#include <TTree.h>
#include <ios>
#include <iostream>
#include <stdexcept>

FW::RootTrackParameterWriter::RootTrackParameterWriter(
    const FW::RootTrackParameterWriter::Config& cfg, Acts::Logging::Level level)
    : TrackParameterWriter(cfg.collection, "RootTrackParameterWriter", level),
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
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr)
    throw std::bad_alloc();
  else {
    // I/O parameters
    m_outputTree->Branch("event_nr", &m_eventNr);
    m_outputTree->Branch("d0", &m_d0);
    m_outputTree->Branch("z0", &m_z0);
    m_outputTree->Branch("phi", &m_phi);
    m_outputTree->Branch("theta", &m_theta);
    m_outputTree->Branch("qp", &m_qp);
    // MORE HERE
  }
}

FW::RootTrackParameterWriter::~RootTrackParameterWriter() {
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

FW::ProcessCode FW::RootTrackParameterWriter::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_outputTree->Write();
    ACTS_INFO("Wrote trackparameters to tree '" << m_cfg.treeName << "' in '"
                                                << m_cfg.filePath << "'");
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode FW::RootTrackParameterWriter::writeT(
    const FW::AlgorithmContext& ctx,
    const std::vector<BoundTrackParameters>& trackParams) {
  if (m_outputFile == nullptr)
    return ProcessCode::SUCCESS;

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = ctx.eventNumber;

  for (auto& params : trackParams) {
    m_d0 = params.parameters()[0];
    m_z0 = params.parameters()[1];
    m_phi = params.parameters()[2];
    m_theta = params.parameters()[3];
    m_qp = params.parameters()[4];

    m_outputTree->Fill();
  }

  return ProcessCode::SUCCESS;
}
