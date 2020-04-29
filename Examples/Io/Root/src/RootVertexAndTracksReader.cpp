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
  m_inputChain->SetBranchAddress("trkCov", &m_ptrTrkCov);

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
  delete m_ptrTrkCov;
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
            std::vector<double> trkCovVec = (*m_ptrTrkCov)[trkId];

            // Construct track covariance
            Acts::BoundSymMatrix covMat =
                Eigen::Map<Acts::BoundSymMatrix>(trkCovVec.data());

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
