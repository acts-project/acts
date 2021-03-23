// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Root/RootVertexAndTracksReader.hpp"

#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/TruthTracking/VertexAndTracks.hpp"

FW::RootVertexAndTracksReader::RootVertexAndTracksReader(
    FW::RootVertexAndTracksReader::Config cfg, Acts::Logging::Level lvl)
    : m_cfg(std::move(cfg)),
      m_events(0),
      m_inputChain(nullptr),
      m_logger(Acts::getDefaultLogger("RootVertexAndTracksReader", lvl)) {
  m_inputChain = new TChain(m_cfg.treeName.c_str());

  m_inputChain->SetBranchAddress("event_nr", &m_eventNr);
  m_inputChain->SetBranchAddress("vx", &m_ptrVx);
  m_inputChain->SetBranchAddress("vy", &m_ptrVy);
  m_inputChain->SetBranchAddress("vz", &m_ptrVz);

  m_inputChain->SetBranchAddress("d0", &m_ptrD0);
  m_inputChain->SetBranchAddress("z0", &m_ptrZ0);
  m_inputChain->SetBranchAddress("phi", &m_ptrPhi);
  m_inputChain->SetBranchAddress("theta", &m_ptrTheta);
  m_inputChain->SetBranchAddress("qp", &m_ptrQP);
  m_inputChain->SetBranchAddress("time", &m_ptrTime);
  m_inputChain->SetBranchAddress("vtxID", &m_ptrVtxID);
  // m_inputChain->SetBranchAddress("trkCov", &m_ptrTrkCov);

  m_inputChain->SetBranchAddress("trkCov11", &m_ptrTrkCov11);
  m_inputChain->SetBranchAddress("trkCov12", &m_ptrTrkCov12);
  m_inputChain->SetBranchAddress("trkCov13", &m_ptrTrkCov13);
  m_inputChain->SetBranchAddress("trkCov14", &m_ptrTrkCov14);
  m_inputChain->SetBranchAddress("trkCov15", &m_ptrTrkCov15);
  m_inputChain->SetBranchAddress("trkCov16", &m_ptrTrkCov16);

  m_inputChain->SetBranchAddress("trkCov21", &m_ptrTrkCov21);
  m_inputChain->SetBranchAddress("trkCov22", &m_ptrTrkCov22);
  m_inputChain->SetBranchAddress("trkCov23", &m_ptrTrkCov23);
  m_inputChain->SetBranchAddress("trkCov24", &m_ptrTrkCov24);
  m_inputChain->SetBranchAddress("trkCov25", &m_ptrTrkCov25);
  m_inputChain->SetBranchAddress("trkCov26", &m_ptrTrkCov26);

  m_inputChain->SetBranchAddress("trkCov31", &m_ptrTrkCov31);
  m_inputChain->SetBranchAddress("trkCov32", &m_ptrTrkCov32);
  m_inputChain->SetBranchAddress("trkCov33", &m_ptrTrkCov33);
  m_inputChain->SetBranchAddress("trkCov34", &m_ptrTrkCov34);
  m_inputChain->SetBranchAddress("trkCov35", &m_ptrTrkCov35);
  m_inputChain->SetBranchAddress("trkCov36", &m_ptrTrkCov36);

  m_inputChain->SetBranchAddress("trkCov41", &m_ptrTrkCov41);
  m_inputChain->SetBranchAddress("trkCov42", &m_ptrTrkCov42);
  m_inputChain->SetBranchAddress("trkCov43", &m_ptrTrkCov43);
  m_inputChain->SetBranchAddress("trkCov44", &m_ptrTrkCov44);
  m_inputChain->SetBranchAddress("trkCov45", &m_ptrTrkCov45);
  m_inputChain->SetBranchAddress("trkCov46", &m_ptrTrkCov46);

  m_inputChain->SetBranchAddress("trkCov51", &m_ptrTrkCov51);
  m_inputChain->SetBranchAddress("trkCov52", &m_ptrTrkCov52);
  m_inputChain->SetBranchAddress("trkCov53", &m_ptrTrkCov53);
  m_inputChain->SetBranchAddress("trkCov54", &m_ptrTrkCov54);
  m_inputChain->SetBranchAddress("trkCov55", &m_ptrTrkCov55);
  m_inputChain->SetBranchAddress("trkCov56", &m_ptrTrkCov56);

  m_inputChain->SetBranchAddress("trkCov61", &m_ptrTrkCov61);
  m_inputChain->SetBranchAddress("trkCov62", &m_ptrTrkCov62);
  m_inputChain->SetBranchAddress("trkCov63", &m_ptrTrkCov63);
  m_inputChain->SetBranchAddress("trkCov64", &m_ptrTrkCov64);
  m_inputChain->SetBranchAddress("trkCov65", &m_ptrTrkCov65);
  m_inputChain->SetBranchAddress("trkCov66", &m_ptrTrkCov66);

  // loop over the input files
  for (auto inputFile : m_cfg.fileList) {
    // add file to the input chain
    m_inputChain->Add(inputFile.c_str());
    ACTS_DEBUG("Adding File " << inputFile << " to tree '" << m_cfg.treeName
                              << "'.");
  }

  m_events = m_inputChain->GetEntries();
  ACTS_DEBUG("The full chain has " << m_events << " entries.");
}

FW::RootVertexAndTracksReader::~RootVertexAndTracksReader() {
  delete m_ptrVx;
  delete m_ptrVy;
  delete m_ptrVz;
  delete m_ptrD0;
  delete m_ptrZ0;
  delete m_ptrPhi;
  delete m_ptrTheta;
  delete m_ptrQP;
  delete m_ptrTime;
  delete m_ptrVtxID;
  // delete m_ptrTrkCov;
}

std::string FW::RootVertexAndTracksReader::name() const {
  return "RootVertexAndTracksReader";
}

std::pair<size_t, size_t> FW::RootVertexAndTracksReader::availableEvents()
    const {
  return {0u, m_events};
}

FW::ProcessCode FW::RootVertexAndTracksReader::read(
    const FW::AlgorithmContext& context) {
  ACTS_DEBUG("Trying to read vertex and tracks.");

  if (m_inputChain && context.eventNumber < m_events) {
    // Lock the mutex
    std::lock_guard<std::mutex> lock(m_read_mutex);

    // The collection to be written
    std::vector<FW::VertexAndTracks> mCollection;

    for (size_t ib = 0; ib < m_cfg.batchSize; ++ib) {
      // Read the correct entry: batch size * event_number + ib
      m_inputChain->GetEntry(m_cfg.batchSize * context.eventNumber + ib);
      ACTS_VERBOSE("Reading entry: " << m_cfg.batchSize * context.eventNumber +
                                            ib);

      // Loop over all vertices
      for (size_t idx = 0; idx < m_ptrVx->size(); ++idx) {
        FW::VertexAndTracks vtxAndTracks;
        vtxAndTracks.vertex.position4[0] = (*m_ptrVx)[idx];
        vtxAndTracks.vertex.position4[1] = (*m_ptrVy)[idx];
        vtxAndTracks.vertex.position4[2] = (*m_ptrVz)[idx];
        vtxAndTracks.vertex.position4[3] = 0;

        std::vector<Acts::BoundParameters> tracks;
        // Loop over all tracks in current event
        for (size_t trkId = 0; trkId < m_ptrD0->size(); ++trkId) {
          // Take only tracks that belong to current vertex
          if (static_cast<size_t>((*m_ptrVtxID)[trkId]) == idx) {
            // Get track parameter
            Acts::BoundVector newTrackParams;
            newTrackParams << (*m_ptrD0)[trkId], (*m_ptrZ0)[trkId],
                (*m_ptrPhi)[trkId], (*m_ptrTheta)[trkId], (*m_ptrQP)[trkId],
                (*m_ptrTime)[trkId];

            // Get track covariance vector
            // std::vector<double> trkCovVec = (*m_ptrTrkCov)[trkId];

            // Construct track covariance
            Acts::BoundSymMatrix covMat;

            covMat << m_cov11[trkId], m_cov12[trkId], m_cov13[trkId],
                m_cov14[trkId], m_cov15[trkId], m_cov16[trkId], m_cov21[trkId],
                m_cov22[trkId], m_cov23[trkId], m_cov24[trkId], m_cov25[trkId],
                m_cov26[trkId], m_cov31[trkId], m_cov32[trkId], m_cov33[trkId],
                m_cov34[trkId], m_cov35[trkId], m_cov36[trkId], m_cov41[trkId],
                m_cov42[trkId], m_cov43[trkId], m_cov44[trkId], m_cov45[trkId],
                m_cov46[trkId], m_cov51[trkId], m_cov52[trkId], m_cov53[trkId],
                m_cov54[trkId], m_cov55[trkId], m_cov56[trkId], m_cov61[trkId],
                m_cov62[trkId], m_cov63[trkId], m_cov64[trkId], m_cov65[trkId],
                m_cov66[trkId];

            // std::cout << covMat <<std::endl;

            // Create track parameters and add to track list
            std::shared_ptr<Acts::PerigeeSurface> perigeeSurface =
                Acts::Surface::makeShared<Acts::PerigeeSurface>(
                    Acts::Vector3D(0., 0., 0.));
            tracks.push_back(
                Acts::BoundParameters(context.geoContext, std::move(covMat),
                                      newTrackParams, perigeeSurface));
          }
        }  // End loop over all tracks
        // Set tracks
        vtxAndTracks.tracks = tracks;
        // Add to collection
        mCollection.push_back(std::move(vtxAndTracks));
      }
    }

    // Write to the collection to the EventStore
    context.eventStore.add(m_cfg.outputCollection, std::move(mCollection));
  }
  // Return success flag
  return FW::ProcessCode::SUCCESS;
}
