// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Root/RootVertexAndTracksWriter.hpp"

#include <Acts/Utilities/Helpers.hpp>
#include <TFile.h>
#include <TTree.h>
#include <ios>
#include <stdexcept>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

FW::RootVertexAndTracksWriter::RootVertexAndTracksWriter(
    const FW::RootVertexAndTracksWriter::Config& cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.collection, "RootVertexAndTracksWriter", lvl),
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
    m_outputTree->Branch("vx", &m_ptrVx);
    m_outputTree->Branch("vy", &m_ptrVy);
    m_outputTree->Branch("vz", &m_ptrVz);

    m_outputTree->Branch("d0", &m_ptrD0);
    m_outputTree->Branch("z0", &m_ptrZ0);
    m_outputTree->Branch("phi", &m_ptrPhi);
    m_outputTree->Branch("theta", &m_ptrTheta);
    m_outputTree->Branch("qp", &m_ptrQP);
    m_outputTree->Branch("time", &m_ptrTime);
    m_outputTree->Branch("vtxID", &m_ptrVtxID);

    m_outputTree->Branch("trkCov11", &m_ptrCov11);
    m_outputTree->Branch("trkCov12", &m_ptrCov12);
    m_outputTree->Branch("trkCov13", &m_ptrCov13);
    m_outputTree->Branch("trkCov14", &m_ptrCov14);
    m_outputTree->Branch("trkCov15", &m_ptrCov15);
    m_outputTree->Branch("trkCov16", &m_ptrCov16);

    m_outputTree->Branch("trkCov21", &m_ptrCov21);
    m_outputTree->Branch("trkCov22", &m_ptrCov22);
    m_outputTree->Branch("trkCov23", &m_ptrCov23);
    m_outputTree->Branch("trkCov24", &m_ptrCov24);
    m_outputTree->Branch("trkCov25", &m_ptrCov25);
    m_outputTree->Branch("trkCov26", &m_ptrCov26);

    m_outputTree->Branch("trkCov31", &m_ptrCov31);
    m_outputTree->Branch("trkCov32", &m_ptrCov32);
    m_outputTree->Branch("trkCov33", &m_ptrCov33);
    m_outputTree->Branch("trkCov34", &m_ptrCov34);
    m_outputTree->Branch("trkCov35", &m_ptrCov35);
    m_outputTree->Branch("trkCov36", &m_ptrCov36);

    m_outputTree->Branch("trkCov41", &m_ptrCov41);
    m_outputTree->Branch("trkCov42", &m_ptrCov42);
    m_outputTree->Branch("trkCov43", &m_ptrCov43);
    m_outputTree->Branch("trkCov44", &m_ptrCov44);
    m_outputTree->Branch("trkCov45", &m_ptrCov45);
    m_outputTree->Branch("trkCov46", &m_ptrCov46);

    m_outputTree->Branch("trkCov51", &m_ptrCov51);
    m_outputTree->Branch("trkCov52", &m_ptrCov52);
    m_outputTree->Branch("trkCov53", &m_ptrCov53);
    m_outputTree->Branch("trkCov54", &m_ptrCov54);
    m_outputTree->Branch("trkCov55", &m_ptrCov55);
    m_outputTree->Branch("trkCov56", &m_ptrCov56);

    m_outputTree->Branch("trkCov61", &m_ptrCov61);
    m_outputTree->Branch("trkCov62", &m_ptrCov62);
    m_outputTree->Branch("trkCov63", &m_ptrCov63);
    m_outputTree->Branch("trkCov64", &m_ptrCov64);
    m_outputTree->Branch("trkCov65", &m_ptrCov65);
    m_outputTree->Branch("trkCov66", &m_ptrCov66);
  }
}

FW::RootVertexAndTracksWriter::~RootVertexAndTracksWriter() {
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

FW::ProcessCode FW::RootVertexAndTracksWriter::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_outputTree->Write();
    ACTS_INFO("Wrote event to tree '" << m_cfg.treeName << "' in '"
                                      << m_cfg.filePath << "'");
  }
  return ProcessCode::SUCCESS;
}

void FW::RootVertexAndTracksWriter::ClearAll() {
  m_vx.clear();
  m_vy.clear();
  m_vz.clear();
  m_d0.clear();
  m_z0.clear();
  m_phi.clear();
  m_theta.clear();
  m_qp.clear();
  m_time.clear();
  m_vtxID.clear();

  m_cov11.clear();
  m_cov12.clear();
  m_cov13.clear();
  m_cov14.clear();
  m_cov15.clear();
  m_cov16.clear();

  m_cov21.clear();
  m_cov22.clear();
  m_cov23.clear();
  m_cov24.clear();
  m_cov25.clear();
  m_cov26.clear();

  m_cov31.clear();
  m_cov32.clear();
  m_cov33.clear();
  m_cov34.clear();
  m_cov35.clear();
  m_cov36.clear();

  m_cov41.clear();
  m_cov42.clear();
  m_cov43.clear();
  m_cov44.clear();
  m_cov45.clear();
  m_cov46.clear();

  m_cov51.clear();
  m_cov52.clear();
  m_cov53.clear();
  m_cov54.clear();
  m_cov55.clear();
  m_cov56.clear();

  m_cov61.clear();
  m_cov62.clear();
  m_cov63.clear();
  m_cov64.clear();
  m_cov65.clear();
  m_cov66.clear();
}

FW::ProcessCode FW::RootVertexAndTracksWriter::writeT(
    const AlgorithmContext& context,
    const std::vector<VertexAndTracks>& vertexAndTracksCollection) {
  if (m_outputFile == nullptr || vertexAndTracksCollection.empty()) {
    return ProcessCode::SUCCESS;
  }

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  ClearAll();

  // Get the event number
  m_eventNr = context.eventNumber;

  for (auto& vertexAndTracks : vertexAndTracksCollection) {
    // Collect the vertex information
    m_vx.push_back(vertexAndTracks.vertex.position().x());
    m_vy.push_back(vertexAndTracks.vertex.position().y());
    m_vz.push_back(vertexAndTracks.vertex.position().z());

    for (auto& track : vertexAndTracks.tracks) {
      // Collect the track information
      m_d0.push_back(track.parameters()[Acts::ParDef::eLOC_D0]);
      m_z0.push_back(track.parameters()[Acts::ParDef::eLOC_Z0]);
      m_phi.push_back(track.parameters()[Acts::ParDef::ePHI]);
      m_theta.push_back(track.parameters()[Acts::ParDef::eTHETA]);
      m_qp.push_back(track.parameters()[Acts::ParDef::eQOP]);
      m_time.push_back(track.parameters()[Acts::ParDef::eT]);
      // Current vertex index as vertex ID
      m_vtxID.push_back(m_vx.size() - 1);

      // Save track covariance
      Acts::BoundSymMatrix cov = *track.covariance();

      m_cov11.push_back(cov(0, 0));
      m_cov12.push_back(cov(0, 1));
      m_cov13.push_back(cov(0, 2));
      m_cov14.push_back(cov(0, 3));
      m_cov15.push_back(cov(0, 4));
      m_cov16.push_back(cov(0, 5));

      m_cov21.push_back(cov(1, 0));
      m_cov22.push_back(cov(1, 1));
      m_cov23.push_back(cov(1, 2));
      m_cov24.push_back(cov(1, 3));
      m_cov25.push_back(cov(1, 4));
      m_cov26.push_back(cov(1, 5));

      m_cov31.push_back(cov(2, 0));
      m_cov32.push_back(cov(2, 1));
      m_cov33.push_back(cov(2, 2));
      m_cov34.push_back(cov(2, 3));
      m_cov35.push_back(cov(2, 4));
      m_cov36.push_back(cov(2, 5));

      m_cov41.push_back(cov(3, 0));
      m_cov42.push_back(cov(3, 1));
      m_cov43.push_back(cov(3, 2));
      m_cov44.push_back(cov(3, 3));
      m_cov45.push_back(cov(3, 4));
      m_cov46.push_back(cov(3, 5));

      m_cov51.push_back(cov(4, 0));
      m_cov52.push_back(cov(4, 1));
      m_cov53.push_back(cov(4, 2));
      m_cov54.push_back(cov(4, 3));
      m_cov55.push_back(cov(4, 4));
      m_cov56.push_back(cov(4, 5));

      m_cov61.push_back(cov(5, 0));
      m_cov62.push_back(cov(5, 1));
      m_cov63.push_back(cov(5, 2));
      m_cov64.push_back(cov(5, 3));
      m_cov65.push_back(cov(5, 4));
      m_cov66.push_back(cov(5, 5));
    }
  }

  m_outputTree->Fill();

  return ProcessCode::SUCCESS;
}
