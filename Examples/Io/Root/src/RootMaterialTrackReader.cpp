// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMaterialTrackReader.hpp"

#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>

ActsExamples::RootMaterialTrackReader::RootMaterialTrackReader(
    const Config& config, Acts::Logging::Level level)
    : ActsExamples::IReader(),
      m_logger{Acts::getDefaultLogger(name(), level)},
      m_cfg(config) {
  m_inputChain = new TChain(m_cfg.treeName.c_str());

  // Set the branches
  m_inputChain->SetBranchAddress("event_id", &m_eventId);
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
  m_inputChain->SetBranchAddress("mat_dx", &m_step_dx);
  m_inputChain->SetBranchAddress("mat_dy", &m_step_dy);
  m_inputChain->SetBranchAddress("mat_dz", &m_step_dz);
  m_inputChain->SetBranchAddress("mat_step_length", &m_step_length);
  m_inputChain->SetBranchAddress("mat_X0", &m_step_X0);
  m_inputChain->SetBranchAddress("mat_L0", &m_step_L0);
  m_inputChain->SetBranchAddress("mat_A", &m_step_A);
  m_inputChain->SetBranchAddress("mat_Z", &m_step_Z);
  m_inputChain->SetBranchAddress("mat_rho", &m_step_rho);
  if (m_cfg.readCachedSurfaceInformation) {
    m_inputChain->SetBranchAddress("sur_id", &m_sur_id);
    m_inputChain->SetBranchAddress("sur_x", &m_sur_x);
    m_inputChain->SetBranchAddress("sur_y", &m_sur_y);
    m_inputChain->SetBranchAddress("sur_z", &m_sur_z);
    m_inputChain->SetBranchAddress("sur_pathCorrection", &m_sur_pathCorrection);
  }
  if (m_cfg.fileList.empty()) {
    throw std::invalid_argument{"No input files given"};
  }

  // loop over the input files
  for (const auto& inputFile : m_cfg.fileList) {
    // add file to the input chain
    m_inputChain->Add(inputFile.c_str());
    ACTS_DEBUG("Adding File " << inputFile << " to tree '" << m_cfg.treeName
                              << "'.");
  }

  m_events = static_cast<size_t>(m_inputChain->GetMaximum("event_id") + 1);
  size_t nentries = m_inputChain->GetEntries();
  m_batchSize = nentries / m_events;
  ACTS_DEBUG("The full chain has "
             << nentries << " entries for " << m_events
             << " events this corresponds to a batch size of: " << m_batchSize);
  std::cout << "The full chain has " << nentries << " entries for " << m_events
            << " events this corresponds to a batch size of: " << m_batchSize
            << std::endl;

  // If the events are not in order, get the entry numbers for ordered events
  if (not m_cfg.orderedEvents) {
    m_entryNumbers.resize(nentries);
    m_inputChain->Draw("event_id", "", "goff");
    // Sort to get the entry numbers of the ordered events
    TMath::Sort(m_inputChain->GetEntries(), m_inputChain->GetV1(),
                m_entryNumbers.data(), false);
  }
  m_outputMaterialTracks.initialize(m_cfg.collection);
}

ActsExamples::RootMaterialTrackReader::~RootMaterialTrackReader() {
  delete m_inputChain;

  delete m_step_x;
  delete m_step_y;
  delete m_step_z;
  delete m_step_dx;
  delete m_step_dy;
  delete m_step_dz;
  delete m_step_length;
  delete m_step_X0;
  delete m_step_L0;
  delete m_step_A;
  delete m_step_Z;
  delete m_step_rho;

  delete m_sur_id;
  delete m_sur_x;
  delete m_sur_y;
  delete m_sur_z;
  delete m_sur_pathCorrection;
}

std::string ActsExamples::RootMaterialTrackReader::name() const {
  return "RootMaterialTrackReader";
}

std::pair<size_t, size_t>
ActsExamples::RootMaterialTrackReader::availableEvents() const {
  return {0u, m_events};
}

ActsExamples::ProcessCode ActsExamples::RootMaterialTrackReader::read(
    const ActsExamples::AlgorithmContext& context) {
  ACTS_DEBUG("Trying to read recorded material from tracks.");
  // read in the material track
  if (m_inputChain != nullptr && context.eventNumber < m_events) {
    // lock the mutex
    std::lock_guard<std::mutex> lock(m_read_mutex);
    // now read

    // The collection to be written
    std::unordered_map<size_t, Acts::RecordedMaterialTrack> mtrackCollection;

    // Loop over the entries for this event
    for (size_t ib = 0; ib < m_batchSize; ++ib) {
      // Read the correct entry: startEntry + ib
      auto entry = m_batchSize * context.eventNumber + ib;
      if (not m_cfg.orderedEvents and entry < m_entryNumbers.size()) {
        entry = m_entryNumbers[entry];
      }
      ACTS_VERBOSE("Reading event: " << context.eventNumber
                                     << " with stored entry: " << entry);
      m_inputChain->GetEntry(entry);

      Acts::RecordedMaterialTrack rmTrack;
      // Fill the position and momentum
      rmTrack.first.first = Acts::Vector3(m_v_x, m_v_y, m_v_z);
      rmTrack.first.second = Acts::Vector3(m_v_px, m_v_py, m_v_pz);

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
            Acts::Vector3((*m_step_x)[is], (*m_step_y)[is], (*m_step_z)[is]);
        mInteraction.direction =
            Acts::Vector3((*m_step_dx)[is], (*m_step_dy)[is], (*m_step_dz)[is]);
        mInteraction.materialSlab = Acts::MaterialSlab(
            Acts::Material::fromMassDensity(mX0, mL0, (*m_step_A)[is],
                                            (*m_step_Z)[is], (*m_step_rho)[is]),
            s);
        if (m_cfg.readCachedSurfaceInformation) {
          // add the surface information to the interaction this allows the
          // mapping to be speed up
          mInteraction.intersectionID =
              Acts::GeometryIdentifier((*m_sur_id)[is]);
          mInteraction.intersection =
              Acts::Vector3((*m_sur_x)[is], (*m_sur_y)[is], (*m_sur_z)[is]);
          mInteraction.pathCorrection = (*m_sur_pathCorrection)[is];
        } else {
          mInteraction.intersectionID = Acts::GeometryIdentifier();
          mInteraction.intersection = Acts::Vector3(0, 0, 0);
        }
        rmTrack.second.materialInteractions.push_back(std::move(mInteraction));
      }
      mtrackCollection[ib] = (std::move(rmTrack));
    }
    // Write to the collection to the EventStore
    m_outputMaterialTracks(context, std::move(mtrackCollection));
  }
  // Return success flag
  return ActsExamples::ProcessCode::SUCCESS;
}
