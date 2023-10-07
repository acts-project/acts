// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootTrajectorySummaryReader.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>

ActsExamples::RootTrajectorySummaryReader::RootTrajectorySummaryReader(
    const ActsExamples::RootTrajectorySummaryReader::Config& config,
    Acts::Logging::Level level)
    : ActsExamples::IReader(),
      m_logger{Acts::getDefaultLogger(name(), level)},
      m_cfg(config) {
  m_inputChain = new TChain(m_cfg.treeName.c_str());

  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing input filename");
  }

  m_outputTrackParameters.initialize(m_cfg.outputTracks);
  m_outputParticles.initialize(m_cfg.outputParticles);

  // Set the branches
  m_inputChain->SetBranchAddress("event_nr", &m_eventNr);
  m_inputChain->SetBranchAddress("multiTraj_nr", &m_multiTrajNr);
  m_inputChain->SetBranchAddress("subTraj_nr", &m_subTrajNr);

  // These info is not really stored in the event store, but still read in
  m_inputChain->SetBranchAddress("nStates", &m_nStates);
  m_inputChain->SetBranchAddress("nMeasurements", &m_nMeasurements);
  m_inputChain->SetBranchAddress("nOutliers", &m_nOutliers);
  m_inputChain->SetBranchAddress("nHoles", &m_nHoles);
  m_inputChain->SetBranchAddress("chi2Sum", &m_chi2Sum);
  m_inputChain->SetBranchAddress("NDF", &m_NDF);
  m_inputChain->SetBranchAddress("measurementChi2", &m_measurementChi2);
  m_inputChain->SetBranchAddress("outlierChi2", &m_outlierChi2);
  m_inputChain->SetBranchAddress("measurementVolume", &m_measurementVolume);
  m_inputChain->SetBranchAddress("measurementLayer", &m_measurementLayer);
  m_inputChain->SetBranchAddress("outlierVolume", &m_outlierVolume);
  m_inputChain->SetBranchAddress("outlierLayer", &m_outlierLayer);

  m_inputChain->SetBranchAddress("majorityParticleId", &m_majorityParticleId);
  m_inputChain->SetBranchAddress("nMajorityHits", &m_nMajorityHits);
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

  auto path = m_cfg.filePath;

  // add file to the input chain
  m_inputChain->Add(path.c_str());
  ACTS_DEBUG("Adding File " << path << " to tree '" << m_cfg.treeName << "'.");

  m_events = m_inputChain->GetEntries();
  ACTS_DEBUG("The full chain has " << m_events << " entries.");

  // If the events are not in order, get the entry numbers for ordered events
  if (not m_cfg.orderedEvents) {
    m_entryNumbers.resize(m_events);
    m_inputChain->Draw("event_nr", "", "goff");
    // Sort to get the entry numbers of the ordered events
    TMath::Sort(m_inputChain->GetEntries(), m_inputChain->GetV1(),
                m_entryNumbers.data(), false);
  }
}

std::pair<size_t, size_t>
ActsExamples::RootTrajectorySummaryReader::availableEvents() const {
  return {0u, m_events};
}

ActsExamples::RootTrajectorySummaryReader::~RootTrajectorySummaryReader() {
  delete m_multiTrajNr;
  delete m_subTrajNr;
  delete m_nStates;
  delete m_nMeasurements;
  delete m_nOutliers;
  delete m_nHoles;
  delete m_chi2Sum;
  delete m_NDF;
  delete m_measurementChi2;
  delete m_outlierChi2;
  delete m_measurementVolume;
  delete m_measurementLayer;
  delete m_outlierVolume;
  delete m_outlierLayer;
  delete m_majorityParticleId;
  delete m_nMajorityHits;
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

ActsExamples::ProcessCode ActsExamples::RootTrajectorySummaryReader::read(
    const ActsExamples::AlgorithmContext& context) {
  ACTS_DEBUG("Trying to read recorded tracks.");

  // read in the fitted track parameters and particles
  if (m_inputChain != nullptr && context.eventNumber < m_events) {
    // lock the mutex
    std::lock_guard<std::mutex> lock(m_read_mutex);
    // now read

    std::shared_ptr<Acts::PerigeeSurface> perigeeSurface =
        Acts::Surface::makeShared<Acts::PerigeeSurface>(
            Acts::Vector3(0., 0., 0.));

    // The collection to be written
    std::vector<Acts::BoundTrackParameters> trackParameterCollection;
    SimParticleContainer truthParticleCollection;

    // Read the correct entry
    auto entry = context.eventNumber;
    if (not m_cfg.orderedEvents and entry < m_entryNumbers.size()) {
      entry = m_entryNumbers[entry];
    }
    m_inputChain->GetEntry(entry);
    ACTS_INFO("Reading event: " << context.eventNumber
                                << " stored as entry: " << entry);

    unsigned int nTracks = m_eLOC0_fit->size();
    for (unsigned int i = 0; i < nTracks; i++) {
      Acts::BoundVector paramVec;
      paramVec << (*m_eLOC0_fit)[i], (*m_eLOC1_fit)[i], (*m_ePHI_fit)[i],
          (*m_eTHETA_fit)[i], (*m_eQOP_fit)[i], (*m_eT_fit)[i];

      // Resolutions
      double resD0 = (*m_err_eLOC0_fit)[i];
      double resZ0 = (*m_err_eLOC1_fit)[i];
      double resPh = (*m_err_ePHI_fit)[i];
      double resTh = (*m_err_eTHETA_fit)[i];
      double resQp = (*m_err_eQOP_fit)[i];
      double resT = (*m_err_eT_fit)[i];

      // Fill vector of track objects with simple covariance matrix
      Acts::BoundSymMatrix covMat;

      covMat << resD0 * resD0, 0., 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0.,
          0., 0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., 0., resTh * resTh,
          0., 0., 0., 0., 0., 0., resQp * resQp, 0., 0., 0., 0., 0., 0.,
          resT * resT;

      trackParameterCollection.push_back(Acts::BoundTrackParameters(
          perigeeSurface, paramVec, std::move(covMat)));
    }

    unsigned int nTruthParticles = m_t_vx->size();
    for (unsigned int i = 0; i < nTruthParticles; i++) {
      ActsFatras::Particle truthParticle;

      truthParticle.setPosition4((*m_t_vx)[i], (*m_t_vy)[i], (*m_t_vz)[i],
                                 (*m_t_time)[i]);
      truthParticle.setDirection((*m_t_px)[i], (*m_t_py)[i], (*m_t_pz)[i]);
      truthParticle.setParticleId((*m_majorityParticleId)[i]);

      truthParticleCollection.insert(truthParticleCollection.end(),
                                     truthParticle);
    }
    // Write the collections to the EventStore
    m_outputTrackParameters(context, std::move(trackParameterCollection));
    m_outputParticles(context, std::move(truthParticleCollection));
  } else {
    ACTS_WARNING("Could not read in event.");
  }
  // Return success flag
  return ActsExamples::ProcessCode::SUCCESS;
}
