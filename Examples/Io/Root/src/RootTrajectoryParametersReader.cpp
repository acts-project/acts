// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootTrajectoryParametersReader.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <iostream>

#include <TChain.h>
#include <TFile.h>

ActsExamples::RootTrajectoryParametersReader::RootTrajectoryParametersReader(
    const ActsExamples::RootTrajectoryParametersReader::Config& cfg)
    : ActsExamples::IReader(), m_cfg(cfg), m_events(0), m_inputChain(nullptr) {
  m_inputChain = new TChain(m_cfg.treeName.c_str());

  // Set the branches
  m_inputChain->SetBranchAddress("event_nr", &m_eventNr);
  m_inputChain->SetBranchAddress("multiTraj_nr", &m_multiTrajNr);
  m_inputChain->SetBranchAddress("subTraj_nr", &m_subTrajNr);
  m_inputChain->SetBranchAddress("t_barcode", &m_t_barcode);
  m_inputChain->SetBranchAddress("t_charge", &m_t_charge);
  m_inputChain->SetBranchAddress("t_time", &m_t_time);
  m_inputChain->SetBranchAddress("t_vx", &m_t_vx);
  m_inputChain->SetBranchAddress("t_vy", &m_t_vy);
  m_inputChain->SetBranchAddress("t_vz", &m_t_vz);
  m_inputChain->SetBranchAddress("t_px", &m_t_px);
  m_inputChain->SetBranchAddress("t_py", &m_t_py);
  m_inputChain->SetBranchAddress("t_pz", &m_t_pz);
  m_inputChain->SetBranchAddress("t_theta", &m_t_theta);
  m_inputChain->SetBranchAddress("t_phi", &m_t_phi);
  m_inputChain->SetBranchAddress("t_eta", &m_t_eta);
  m_inputChain->SetBranchAddress("t_pT", &m_t_pT);

  m_inputChain->SetBranchAddress("hasFittedParams", &m_hasFittedParams);
  m_inputChain->SetBranchAddress("eLOC0_fit", &m_eLOC0_fit);
  m_inputChain->SetBranchAddress("eLOC1_fit", &m_eLOC1_fit);
  m_inputChain->SetBranchAddress("ePHI_fit", &m_ePHI_fit);
  m_inputChain->SetBranchAddress("eTHETA_fit", &m_eTHETA_fit);
  m_inputChain->SetBranchAddress("eQOP_fit", &m_eQOP_fit);
  m_inputChain->SetBranchAddress("eT_fit", &m_eT_fit);
  m_inputChain->SetBranchAddress("err_eLOC0_fit", &m_err_eLOC0_fit);
  m_inputChain->SetBranchAddress("err_eLOC1_fit", &m_err_eLOC1_fit);
  m_inputChain->SetBranchAddress("err_ePHI_fit", &m_err_ePHI_fit);
  m_inputChain->SetBranchAddress("err_eTHETA_fit", &m_err_eTHETA_fit);
  m_inputChain->SetBranchAddress("err_eQOP_fit", &m_err_eQOP_fit);
  m_inputChain->SetBranchAddress("err_eT_fit", &m_err_eT_fit);

  // add file to the input chain
  m_inputChain->Add(m_cfg.inputFile.c_str());
  ACTS_DEBUG("Adding File " << m_cfg.inputFile << " to tree '" << m_cfg.treeName
                            << "'.");


  m_events = m_inputChain->GetEntries();
  ACTS_DEBUG("The full chain has " << m_events << " entries.");
}

std::string ActsExamples::RootTrajectoryParametersReader::name() const {
  return m_cfg.name;
}

std::pair<size_t, size_t>
ActsExamples::RootTrajectoryParametersReader::availableEvents() const {
  return {0u, m_events};
}

ActsExamples::RootTrajectoryParametersReader::~RootTrajectoryParametersReader() {
  delete m_multiTrajNr;
  delete m_subTrajNr;
  delete m_t_barcode;
  delete m_t_charge;
  delete m_t_time;
  delete m_t_vx;
  delete m_t_vy;
  delete m_t_vz;
  delete m_t_px;
  delete m_t_py;
  delete m_t_pz;
  delete m_t_theta;
  delete m_t_phi;
  delete m_t_pT;
  delete m_t_eta;
  delete m_hasFittedParams;
  delete m_eLOC0_fit;
  delete m_eLOC1_fit;
  delete m_ePHI_fit;
  delete m_eTHETA_fit;
  delete m_eQOP_fit;
  delete m_eT_fit;
  delete m_err_eLOC0_fit;
  delete m_err_eLOC1_fit;
  delete m_err_ePHI_fit;
  delete m_err_eTHETA_fit;
  delete m_err_eQOP_fit;
  delete m_err_eT_fit;
}

ActsExamples::ProcessCode ActsExamples::RootTrajectoryParametersReader::read(
    const ActsExamples::AlgorithmContext& context) {
  ACTS_DEBUG("Trying to read recorded tracks.");

  std::cout << "Trying to read recorded tracks now." << std::endl;
  // read in the material track
  if (m_inputChain && context.eventNumber < m_events) {
    // lock the mutex
    std::lock_guard<std::mutex> lock(m_read_mutex);
    // now read

    std::shared_ptr<Acts::PerigeeSurface> perigeeSurface =
        Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3(0., 0., 0.));

    // The collection to be written
    std::vector<Acts::BoundTrackParameters> mtrackCollection;
    std::vector<ActsFatras::Particle> mTruthParticleCollection;
  
    // Read the correct entry: batch size * event_number + ib
    m_inputChain->GetEntry(context.eventNumber);
    std::cout << "Reading track from event " << m_eventNr << std::endl;
    ACTS_VERBOSE("Reading entry: " << context.eventNumber);

    unsigned int nTracks = m_eLOC0_fit->size();
    for(unsigned int i = 0; i < nTracks; i++){

      Acts::BoundVector paramVec;
      paramVec << (*m_eLOC0_fit)[i], (*m_eLOC1_fit)[i], (*m_ePHI_fit)[i], (*m_eTHETA_fit)[i], (*m_eQOP_fit)[i], (*m_eT_fit)[i];

      //std::cout << "param vector:  " <<  paramVec << std::endl;

      // Resolutions
      double resD0 = (*m_err_eLOC0_fit)[i];
      double resZ0 = (*m_err_eLOC1_fit)[i];
      double resPh = (*m_err_ePHI_fit)[i];
      double resTh = (*m_err_eTHETA_fit)[i];
      double resQp = (*m_err_eQOP_fit)[i];
      double resT  = (*m_err_eT_fit)[i];

      // Fill vector of track objects with simple covariance matrix
      Acts::BoundSymMatrix covMat;

      covMat << resD0 * resD0, 0., 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0.,
          0., 0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., 0., resTh * resTh,
          0., 0., 0., 0., 0., 0., resQp * resQp, 0., 0., 0., 0., 0., 0., resT*resT;
   
      
      mtrackCollection.push_back(Acts::BoundTrackParameters(perigeeSurface, paramVec, std::move(covMat)));
    }

    unsigned int nTruthParticles = m_t_vx->size();
    for(unsigned int i = 0; i < nTruthParticles; i++){
      ActsFatras::Particle truthParticle;

      truthParticle.setPosition4((*m_t_vx)[i], (*m_t_vy)[i], (*m_t_vz)[i], (*m_t_time)[i]);
      truthParticle.setDirection((*m_t_px)[i], (*m_t_py)[i], (*m_t_pz)[i]);

      mTruthParticleCollection.push_back(truthParticle);
    }

    std::cout << "event done... writing " << mtrackCollection.size() << "/" << mTruthParticleCollection.size() << " tracks/particles to whiteboard..." << std::endl;
    

    // Write the collections to the EventStore
    context.eventStore.add(m_cfg.trackCollection, std::move(mtrackCollection));
    context.eventStore.add(m_cfg.particleCollection, std::move(mTruthParticleCollection));
  }
  std::cout << "exiting read method.."<< std::endl;
  // Return success flag
  return ActsExamples::ProcessCode::SUCCESS;
}
