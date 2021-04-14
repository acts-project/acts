// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Root/RootRecVertexWriter.hpp"

#include <Acts/Utilities/Helpers.hpp>
#include <TFile.h>
#include <TTree.h>
#include <ios>
#include <stdexcept>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

FW::RootRecVertexWriter::RootRecVertexWriter(
    const FW::RootRecVertexWriter::Config& cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.collection, "RootRecVertexWriter", lvl),
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
    m_outputTree->Branch("fitquality", &m_ptr_vtx_fitquality);

    m_outputTree->Branch("d0", &m_ptrD0);
    m_outputTree->Branch("z0", &m_ptrZ0);
    m_outputTree->Branch("phi", &m_ptrPhi);
    m_outputTree->Branch("theta", &m_ptrTheta);
    m_outputTree->Branch("qp", &m_ptrQP);
    m_outputTree->Branch("time", &m_ptrTime);
    m_outputTree->Branch("vtxID", &m_ptrVtxID);

    m_outputTree->Branch("vtxCov11", &m_ptr_vtx_Cov11);
    m_outputTree->Branch("vtxCov12", &m_ptr_vtx_Cov12);
    m_outputTree->Branch("vtxCov13", &m_ptr_vtx_Cov13);
    m_outputTree->Branch("vtxCov14", &m_ptr_vtx_Cov14);
    m_outputTree->Branch("vtxCov21", &m_ptr_vtx_Cov21);
    m_outputTree->Branch("vtxCov22", &m_ptr_vtx_Cov22);
    m_outputTree->Branch("vtxCov23", &m_ptr_vtx_Cov23);
    m_outputTree->Branch("vtxCov24", &m_ptr_vtx_Cov24);
    m_outputTree->Branch("vtxCov31", &m_ptr_vtx_Cov31);
    m_outputTree->Branch("vtxCov32", &m_ptr_vtx_Cov32);
    m_outputTree->Branch("vtxCov33", &m_ptr_vtx_Cov33);
    m_outputTree->Branch("vtxCov34", &m_ptr_vtx_Cov34);
    m_outputTree->Branch("vtxCov41", &m_ptr_vtx_Cov41);
    m_outputTree->Branch("vtxCov42", &m_ptr_vtx_Cov42);
    m_outputTree->Branch("vtxCov43", &m_ptr_vtx_Cov43);
    m_outputTree->Branch("vtxCov44", &m_ptr_vtx_Cov44);

    m_outputTree->Branch("trkCov11", &m_ptr_trk_Cov11);
    m_outputTree->Branch("trkCov12", &m_ptr_trk_Cov12);
    m_outputTree->Branch("trkCov13", &m_ptr_trk_Cov13);
    m_outputTree->Branch("trkCov14", &m_ptr_trk_Cov14);
    m_outputTree->Branch("trkCov15", &m_ptr_trk_Cov15);
    m_outputTree->Branch("trkCov16", &m_ptr_trk_Cov16);

    m_outputTree->Branch("trkCov21", &m_ptr_trk_Cov21);
    m_outputTree->Branch("trkCov22", &m_ptr_trk_Cov22);
    m_outputTree->Branch("trkCov23", &m_ptr_trk_Cov23);
    m_outputTree->Branch("trkCov24", &m_ptr_trk_Cov24);
    m_outputTree->Branch("trkCov25", &m_ptr_trk_Cov25);
    m_outputTree->Branch("trkCov26", &m_ptr_trk_Cov26);

    m_outputTree->Branch("trkCov31", &m_ptr_trk_Cov31);
    m_outputTree->Branch("trkCov32", &m_ptr_trk_Cov32);
    m_outputTree->Branch("trkCov33", &m_ptr_trk_Cov33);
    m_outputTree->Branch("trkCov34", &m_ptr_trk_Cov34);
    m_outputTree->Branch("trkCov35", &m_ptr_trk_Cov35);
    m_outputTree->Branch("trkCov36", &m_ptr_trk_Cov36);

    m_outputTree->Branch("trkCov41", &m_ptr_trk_Cov41);
    m_outputTree->Branch("trkCov42", &m_ptr_trk_Cov42);
    m_outputTree->Branch("trkCov43", &m_ptr_trk_Cov43);
    m_outputTree->Branch("trkCov44", &m_ptr_trk_Cov44);
    m_outputTree->Branch("trkCov45", &m_ptr_trk_Cov45);
    m_outputTree->Branch("trkCov46", &m_ptr_trk_Cov46);

    m_outputTree->Branch("trkCov51", &m_ptr_trk_Cov51);
    m_outputTree->Branch("trkCov52", &m_ptr_trk_Cov52);
    m_outputTree->Branch("trkCov53", &m_ptr_trk_Cov53);
    m_outputTree->Branch("trkCov54", &m_ptr_trk_Cov54);
    m_outputTree->Branch("trkCov55", &m_ptr_trk_Cov55);
    m_outputTree->Branch("trkCov56", &m_ptr_trk_Cov56);

    m_outputTree->Branch("trkCov61", &m_ptr_trk_Cov61);
    m_outputTree->Branch("trkCov62", &m_ptr_trk_Cov62);
    m_outputTree->Branch("trkCov63", &m_ptr_trk_Cov63);
    m_outputTree->Branch("trkCov64", &m_ptr_trk_Cov64);
    m_outputTree->Branch("trkCov65", &m_ptr_trk_Cov65);
    m_outputTree->Branch("trkCov66", &m_ptr_trk_Cov66);
  }
}

FW::RootRecVertexWriter::~RootRecVertexWriter() {
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

FW::ProcessCode FW::RootRecVertexWriter::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_outputTree->Write();
    ACTS_INFO("Wrote event to tree '" << m_cfg.treeName << "' in '"
                                      << m_cfg.filePath << "'");
  }
  return ProcessCode::SUCCESS;
}

void FW::RootRecVertexWriter::ClearAll() {
  m_vx.clear();
  m_vy.clear();
  m_vz.clear();
  m_vtx_fitquality.clear();
  m_d0.clear();
  m_z0.clear();
  m_phi.clear();
  m_theta.clear();
  m_qp.clear();
  m_time.clear();
  m_vtxID.clear();

  m_vtx_cov11.clear();
  m_vtx_cov12.clear();
  m_vtx_cov13.clear();
  m_vtx_cov14.clear();
  m_vtx_cov21.clear();
  m_vtx_cov22.clear();
  m_vtx_cov23.clear();
  m_vtx_cov24.clear();
  m_vtx_cov31.clear();
  m_vtx_cov32.clear();
  m_vtx_cov33.clear();
  m_vtx_cov34.clear();
  m_vtx_cov41.clear();
  m_vtx_cov42.clear();
  m_vtx_cov43.clear();
  m_vtx_cov44.clear();

  m_trk_cov11.clear();
  m_trk_cov12.clear();
  m_trk_cov13.clear();
  m_trk_cov14.clear();
  m_trk_cov15.clear();
  m_trk_cov16.clear();
  m_trk_cov21.clear();
  m_trk_cov22.clear();
  m_trk_cov23.clear();
  m_trk_cov24.clear();
  m_trk_cov25.clear();
  m_trk_cov26.clear();
  m_trk_cov31.clear();
  m_trk_cov32.clear();
  m_trk_cov33.clear();
  m_trk_cov34.clear();
  m_trk_cov35.clear();
  m_trk_cov36.clear();
  m_trk_cov41.clear();
  m_trk_cov42.clear();
  m_trk_cov43.clear();
  m_trk_cov44.clear();
  m_trk_cov45.clear();
  m_trk_cov46.clear();
  m_trk_cov51.clear();
  m_trk_cov52.clear();
  m_trk_cov53.clear();
  m_trk_cov54.clear();
  m_trk_cov55.clear();
  m_trk_cov66.clear();
  m_trk_cov61.clear();
  m_trk_cov62.clear();
  m_trk_cov63.clear();
  m_trk_cov64.clear();
  m_trk_cov65.clear();
  m_trk_cov66.clear();
}

FW::ProcessCode FW::RootRecVertexWriter::writeT(
    const AlgorithmContext& context,
    const std::vector<Acts::Vertex<Acts::BoundParameters>>& vertexCollection) {
  // if (m_outputFile == nullptr || vertexCollection.empty()) {
  //   return ProcessCode::SUCCESS;
  // }

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  ClearAll();

  // Get the event number
  m_eventNr = context.eventNumber;

  ACTS_INFO("Found " << vertexCollection.size() << " vertices in event.");

  unsigned int count = 0;
  for (const auto& vtx : vertexCollection) {
    ACTS_INFO("\t" << ++count << ". vertex at "
                   << "(" << vtx.position().x() << "," << vtx.position().y()
                   << "," << vtx.position().z() << ") with "
                   << vtx.tracks().size() << " tracks.");

    // Collect the vertex information
    m_vx.push_back(vtx.position().x());
    m_vy.push_back(vtx.position().y());
    m_vz.push_back(vtx.position().z());
    m_vtx_fitquality.push_back(vtx.fitQuality());
    // std::cout << vtx.position().z() << "\t" << vtx.fitQuality().first<<
    // std::endl;
    std::cout << "fit size:" << m_vtx_fitquality.size() << std::endl;

    Acts::SpacePointSymMatrix vtx_cov = vtx.fullCovariance();
    m_vtx_cov11.push_back(vtx_cov(0, 0));
    m_vtx_cov12.push_back(vtx_cov(0, 1));
    m_vtx_cov13.push_back(vtx_cov(0, 2));
    m_vtx_cov14.push_back(vtx_cov(0, 3));

    m_vtx_cov21.push_back(vtx_cov(1, 0));
    m_vtx_cov22.push_back(vtx_cov(1, 1));
    m_vtx_cov23.push_back(vtx_cov(1, 2));
    m_vtx_cov24.push_back(vtx_cov(1, 3));

    m_vtx_cov31.push_back(vtx_cov(2, 0));
    m_vtx_cov32.push_back(vtx_cov(2, 1));
    m_vtx_cov33.push_back(vtx_cov(2, 2));
    m_vtx_cov34.push_back(vtx_cov(2, 3));

    m_vtx_cov41.push_back(vtx_cov(3, 0));
    m_vtx_cov42.push_back(vtx_cov(3, 1));
    m_vtx_cov43.push_back(vtx_cov(3, 2));
    m_vtx_cov44.push_back(vtx_cov(3, 3));

    // TODO: Add more branches
    for (auto& track : vtx.tracks()) {
      // Collect the track information
      // Store Fitted perigee tracks
      m_d0.push_back(track.fittedParams.parameters()[Acts::ParDef::eLOC_D0]);
      m_z0.push_back(track.fittedParams.parameters()[Acts::ParDef::eLOC_Z0]);
      m_phi.push_back(track.fittedParams.parameters()[Acts::ParDef::ePHI]);
      m_theta.push_back(track.fittedParams.parameters()[Acts::ParDef::eTHETA]);
      m_qp.push_back(track.fittedParams.parameters()[Acts::ParDef::eQOP]);
      m_time.push_back(track.fittedParams.parameters()[Acts::ParDef::eT]);
      // Current vertex index as vertex ID
      m_vtxID.push_back(vtx.tracks().size() - 1);

      // // Save track covariance
      auto cov = track.fittedParams.covariance();

      m_trk_cov11.push_back((*cov)(0, 0));
      m_trk_cov12.push_back((*cov)(0, 1));
      m_trk_cov13.push_back((*cov)(0, 2));
      m_trk_cov14.push_back((*cov)(0, 3));
      m_trk_cov15.push_back((*cov)(0, 4));
      m_trk_cov16.push_back((*cov)(0, 5));

      m_trk_cov21.push_back((*cov)(1, 0));
      m_trk_cov22.push_back((*cov)(1, 1));
      m_trk_cov23.push_back((*cov)(1, 2));
      m_trk_cov24.push_back((*cov)(1, 3));
      m_trk_cov25.push_back((*cov)(1, 4));
      m_trk_cov26.push_back((*cov)(1, 5));

      m_trk_cov31.push_back((*cov)(2, 0));
      m_trk_cov32.push_back((*cov)(2, 1));
      m_trk_cov33.push_back((*cov)(2, 2));
      m_trk_cov34.push_back((*cov)(2, 3));
      m_trk_cov35.push_back((*cov)(2, 4));
      m_trk_cov36.push_back((*cov)(2, 5));

      m_trk_cov41.push_back((*cov)(3, 0));
      m_trk_cov42.push_back((*cov)(3, 1));
      m_trk_cov43.push_back((*cov)(3, 2));
      m_trk_cov44.push_back((*cov)(3, 3));
      m_trk_cov45.push_back((*cov)(3, 4));
      m_trk_cov46.push_back((*cov)(3, 5));

      m_trk_cov51.push_back((*cov)(4, 0));
      m_trk_cov52.push_back((*cov)(4, 1));
      m_trk_cov53.push_back((*cov)(4, 2));
      m_trk_cov54.push_back((*cov)(4, 3));
      m_trk_cov55.push_back((*cov)(4, 4));
      m_trk_cov56.push_back((*cov)(4, 5));

      m_trk_cov61.push_back((*cov)(5, 0));
      m_trk_cov62.push_back((*cov)(5, 1));
      m_trk_cov63.push_back((*cov)(5, 2));
      m_trk_cov64.push_back((*cov)(5, 3));
      m_trk_cov65.push_back((*cov)(5, 4));
      m_trk_cov66.push_back((*cov)(5, 5));
    }
  }
  std::cout << "fit size:" << m_vtx_fitquality.size() << std::endl;
  std::cout << "vtx cov size:" << m_vtx_cov11.size() << std::endl;

  m_outputTree->Fill();

  return ProcessCode::SUCCESS;
}
