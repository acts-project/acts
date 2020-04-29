// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Root/RootMaterialTrackReader.hpp"

#include <TChain.h>
#include <TFile.h>
#include <iostream>

#include "ACTFW/Framework/WhiteBoard.hpp"

FW::RootMaterialTrackReader::RootMaterialTrackReader(
    const FW::RootMaterialTrackReader::Config& cfg)
    : FW::IReader(), m_cfg(cfg), m_events(0), m_inputChain(nullptr) {
  m_inputChain = new TChain(m_cfg.treeName.c_str());

  // Set the branches
  m_inputChain->SetBranchAddress("v_x", &m_v_x);
  m_inputChain->SetBranchAddress("v_y", &m_v_y);
  m_inputChain->SetBranchAddress("v_z", &m_v_z);
  m_inputChain->SetBranchAddress("v_px", &m_v_px);
  m_inputChain->SetBranchAddress("v_py", &m_v_py);
  m_inputChain->SetBranchAddress("v_pz", &m_v_pz);
  m_inputChain->SetBranchAddress("v_phi", &m_v_phi);
  m_inputChain->SetBranchAddress("v_eta", &m_v_eta);
  m_inputChain->SetBranchAddress("t_X0", &m_tX0);
  m_inputChain->SetBranchAddress("t_L0", &m_tL0);
  m_inputChain->SetBranchAddress("mat_x", &m_step_x);
  m_inputChain->SetBranchAddress("mat_y", &m_step_y);
  m_inputChain->SetBranchAddress("mat_z", &m_step_z);
  m_inputChain->SetBranchAddress("mat_step_length", &m_step_length);
  m_inputChain->SetBranchAddress("mat_X0", &m_step_X0);
  m_inputChain->SetBranchAddress("mat_L0", &m_step_L0);
  m_inputChain->SetBranchAddress("mat_A", &m_step_A);
  m_inputChain->SetBranchAddress("mat_Z", &m_step_Z);
  m_inputChain->SetBranchAddress("mat_rho", &m_step_rho);

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

FW::RootMaterialTrackReader::~RootMaterialTrackReader() {
  delete m_step_x;
  delete m_step_y;
  delete m_step_z;
  delete m_step_length;
  delete m_step_X0;
  delete m_step_L0;
  delete m_step_A;
  delete m_step_Z;
  delete m_step_rho;
}

std::string FW::RootMaterialTrackReader::name() const {
  return m_cfg.name;
}

std::pair<size_t, size_t> FW::RootMaterialTrackReader::availableEvents() const {
  return {0u, m_events};
}

FW::ProcessCode FW::RootMaterialTrackReader::read(
    const FW::AlgorithmContext& context) {
  ACTS_DEBUG("Trying to read recorded material from tracks.");
  // read in the material track
  if (m_inputChain && context.eventNumber < m_events) {
    // lock the mutex
    std::lock_guard<std::mutex> lock(m_read_mutex);
    // now read

    // The collection to be written
    std::vector<Acts::RecordedMaterialTrack> mtrackCollection;

    for (size_t ib = 0; ib < m_cfg.batchSize; ++ib) {
      // Read the correct entry: batch size * event_number + ib
      m_inputChain->GetEntry(m_cfg.batchSize * context.eventNumber + ib);
      ACTS_VERBOSE("Reading entry: " << m_cfg.batchSize * context.eventNumber +
                                            ib);

      Acts::RecordedMaterialTrack rmTrack;
      // Fill the position and momentum
      rmTrack.first.first = Acts::Vector3D(m_v_x, m_v_y, m_v_z);
      rmTrack.first.second = Acts::Vector3D(m_v_px, m_v_py, m_v_pz);

      // Fill the individual steps
      size_t msteps = m_step_length->size();
      ACTS_VERBOSE("Reading " << msteps << " material steps.");
      rmTrack.second.materialInteractions.reserve(msteps);
      rmTrack.second.materialInX0 = 0.;
      rmTrack.second.materialInL0 = 0.;

      for (size_t is = 0; is < msteps; ++is) {
        double mX0 = (*m_step_X0)[is];
        double mL0 = (*m_step_L0)[is];
        double s = (*m_step_length)[is];

        rmTrack.second.materialInX0 += s / mX0;
        rmTrack.second.materialInL0 += s / mL0;

        /// Fill the position & the material
        Acts::MaterialInteraction mInteraction;
        mInteraction.position =
            Acts::Vector3D((*m_step_x)[is], (*m_step_y)[is], (*m_step_z)[is]);
        mInteraction.materialProperties = Acts::MaterialProperties(
            mX0, mL0, (*m_step_A)[is], (*m_step_Z)[is], (*m_step_rho)[is], s);
        rmTrack.second.materialInteractions.push_back(std::move(mInteraction));
      }
      mtrackCollection.push_back(std::move(rmTrack));
    }

    // Write to the collection to the EventStore
    context.eventStore.add(m_cfg.collection, std::move(mtrackCollection));
  }
  // Return success flag
  return FW::ProcessCode::SUCCESS;
}
