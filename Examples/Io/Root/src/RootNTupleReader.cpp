// This file is part of the Acts project.
//
// Copyright (C) 2017-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootNTupleReader.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>

ActsExamples::RootNTupleReader::RootNTupleReader(
    const ActsExamples::RootNTupleReader::Config& config,
    Acts::Logging::Level level)
    : ActsExamples::IReader(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.inputFilePath.empty()) {
    throw std::invalid_argument("Missing input filename");
  }
  if (m_cfg.inputTreeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_inputChain = new TChain(m_cfg.inputTreeName.c_str());

  // unused event identifier
  std::int32_t eventNumber;

  // Set the branches
  m_inputChain->SetBranchAddress("EventNumber", &eventNumber);
  m_inputChain->SetBranchAddress("track_d0", &m_branches.m_track_d0);
  m_inputChain->SetBranchAddress("track_z0", &m_branches.m_track_z0);
  m_inputChain->SetBranchAddress("track_theta", &m_branches.m_track_theta);
  m_inputChain->SetBranchAddress("track_phi", &m_branches.m_track_phi);
  m_inputChain->SetBranchAddress("track_qOverP", &m_branches.m_track_qOverP);
  m_inputChain->SetBranchAddress("track_t", &m_branches.m_track_t);
  m_inputChain->SetBranchAddress("track_z", &m_branches.m_track_z);

  // Covariance stuff
  m_inputChain->SetBranchAddress("track_var_d0", &m_branches.m_track_var_d0);
  m_inputChain->SetBranchAddress("track_var_z0", &m_branches.m_track_var_z0);
  m_inputChain->SetBranchAddress("track_var_phi", &m_branches.m_track_var_phi);
  m_inputChain->SetBranchAddress("track_var_theta",
                                 &m_branches.m_track_var_theta);
  m_inputChain->SetBranchAddress("track_var_qOverP",
                                 &m_branches.m_track_var_qOverP);
  m_inputChain->SetBranchAddress("track_cov_d0z0",
                                 &m_branches.m_track_cov_d0z0);
  m_inputChain->SetBranchAddress("track_cov_d0phi",
                                 &m_branches.m_track_cov_d0phi);
  m_inputChain->SetBranchAddress("track_cov_d0theta",
                                 &m_branches.m_track_cov_d0theta);
  m_inputChain->SetBranchAddress("track_cov_d0qOverP",
                                 &m_branches.m_track_cov_d0qOverP);
  m_inputChain->SetBranchAddress("track_cov_z0phi",
                                 &m_branches.m_track_cov_z0phi);
  m_inputChain->SetBranchAddress("track_cov_z0theta",
                                 &m_branches.m_track_cov_z0theta);
  m_inputChain->SetBranchAddress("track_cov_z0qOverP",
                                 &m_branches.m_track_cov_z0qOverP);
  m_inputChain->SetBranchAddress("track_cov_phitheta",
                                 &m_branches.m_track_cov_phitheta);
  m_inputChain->SetBranchAddress("track_cov_phiqOverP",
                                 &m_branches.m_track_cov_phiqOverP);
  m_inputChain->SetBranchAddress("track_cov_tehtaqOverP",
                                 &m_branches.m_track_cov_tehtaqOverP);

  // Truth vertex
  m_inputChain->SetBranchAddress("truthvertex_x", &m_branches.m_truthvertex_x);
  m_inputChain->SetBranchAddress("truthvertex_y", &m_branches.m_truthvertex_y);
  m_inputChain->SetBranchAddress("truthvertex_z", &m_branches.m_truthvertex_z);
  m_inputChain->SetBranchAddress("truthvertex_t", &m_branches.m_truthvertex_t);

  m_inputChain->SetBranchAddress("recovertex_x", &m_branches.m_recovertex_x);
  m_inputChain->SetBranchAddress("recovertex_y", &m_branches.m_recovertex_y);
  m_inputChain->SetBranchAddress("recovertex_z", &m_branches.m_recovertex_z);
  m_inputChain->SetBranchAddress("truthvertex_tracks_idx",
                                 &m_branches.m_truthvertex_tracks_idx);

  m_inputChain->SetBranchAddress("beamspot_x", &m_branches.m_beamspot_x);
  m_inputChain->SetBranchAddress("beamspot_y", &m_branches.m_beamspot_y);
  m_inputChain->SetBranchAddress("beamspot_z", &m_branches.m_beamspot_z);
  m_inputChain->SetBranchAddress("beamspot_sigX", &m_branches.m_beamspot_sigX);
  m_inputChain->SetBranchAddress("beamspot_sigY", &m_branches.m_beamspot_sigY);
  m_inputChain->SetBranchAddress("beamspot_sigZ", &m_branches.m_beamspot_sigZ);

  auto path = m_cfg.inputFilePath;

  // add file to the input chain
  m_inputChain->Add(path.c_str());
  ACTS_DEBUG("Adding File " << path << " to tree '" << m_cfg.inputTreeName
                            << "'.");

  m_events = m_inputChain->GetEntries();
  ACTS_DEBUG("The full chain has " << m_events << " entries.");

  {
    // The entry numbers for accessing events in increased order (there could be
    // multiple entries corresponding to one event number)
    std::vector<long long> m_entryNumbers;
    // If the events are not in order, get the entry numbers for ordered events
    m_entryNumbers.resize(m_events);
    m_inputChain->Draw("EventNumber", "", "goff");
    // Sort to get the entry numbers of the ordered events
    TMath::Sort(m_inputChain->GetEntries(), m_inputChain->GetV1(),
                m_entryNumbers.data(), false);
  }
}

std::pair<std::size_t, std::size_t>
ActsExamples::RootNTupleReader::availableEvents() const {
  return {0u, m_events};
}

ActsExamples::ProcessCode ActsExamples::RootNTupleReader::read(
    const ActsExamples::AlgorithmContext& context) {
  ACTS_DEBUG("Trying to read track parameters from ntuple.");

  Acts::Vector3 pos(0, 0, 0);
  std::shared_ptr<Acts::PerigeeSurface> surface =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(pos);

  if (context.eventNumber >= m_events) {
    ACTS_ERROR("event out of bounds");
  }

  std::lock_guard<std::mutex> lock(m_read_mutex);

  auto entry = context.eventNumber;
  m_inputChain->GetEntry(entry);
  ACTS_INFO("Reading event: " << context.eventNumber
                              << " stored as entry: " << entry);

  const unsigned int nTracks = m_branches.m_track_d0->size();
  const unsigned int nTruthVtx = m_branches.m_truthvertex_z->size();
  const unsigned int nRecoVtx = m_branches.m_recovertex_z->size();

  ACTS_DEBUG("nTracks = " << nTracks);
  ACTS_DEBUG("nTruthVtx = " << nTruthVtx);
  ACTS_DEBUG("nRecoVtx = " << nRecoVtx);

  TrackParametersContainer trackContainer;
  trackContainer.reserve(nTracks);

  for (unsigned int i = 0; i < nTracks; i++) {
    Acts::BoundVector params;

    params[Acts::BoundIndices::eBoundLoc0] = (*m_branches.m_track_d0)[i];
    params[Acts::BoundIndices::eBoundLoc1] = (*m_branches.m_track_z0)[i];
    params[Acts::BoundIndices::eBoundPhi] = (*m_branches.m_track_phi)[i];
    params[Acts::BoundIndices::eBoundTheta] = (*m_branches.m_track_theta)[i];
    params[Acts::BoundIndices::eBoundQOverP] = (*m_branches.m_track_qOverP)[i];
    params[Acts::BoundIndices::eBoundTime] = (*m_branches.m_track_t)[i];

    const double q = 1;

    // Construct and fill covariance matrix
    Acts::BoundSymMatrix cov;

    // Variances
    cov(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundLoc0) =
        (*m_branches.m_track_var_d0)[i];
    cov(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundLoc1) =
        (*m_branches.m_track_var_z0)[i];
    cov(Acts::BoundIndices::eBoundPhi, Acts::BoundIndices::eBoundPhi) =
        (*m_branches.m_track_var_phi)[i];
    cov(Acts::BoundIndices::eBoundTheta, Acts::BoundIndices::eBoundTheta) =
        (*m_branches.m_track_var_theta)[i];
    cov(Acts::BoundIndices::eBoundQOverP, Acts::BoundIndices::eBoundQOverP) =
        (*m_branches.m_track_var_qOverP)[i];

    cov(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundLoc1) =
        (*m_branches.m_track_cov_d0z0)[i];
    cov(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundPhi) =
        (*m_branches.m_track_cov_d0phi)[i];
    cov(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundTheta) =
        (*m_branches.m_track_cov_d0theta)[i];
    cov(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundQOverP) =
        (*m_branches.m_track_cov_d0qOverP)[i];
    cov(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundPhi) =
        (*m_branches.m_track_cov_z0phi)[i];
    cov(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundTheta) =
        (*m_branches.m_track_cov_z0theta)[i];
    cov(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundQOverP) =
        (*m_branches.m_track_cov_z0qOverP)[i];
    cov(Acts::BoundIndices::eBoundPhi, Acts::BoundIndices::eBoundTheta) =
        (*m_branches.m_track_cov_phitheta)[i];
    cov(Acts::BoundIndices::eBoundPhi, Acts::BoundIndices::eBoundQOverP) =
        (*m_branches.m_track_cov_phiqOverP)[i];
    cov(Acts::BoundIndices::eBoundTheta, Acts::BoundIndices::eBoundQOverP) =
        (*m_branches.m_track_cov_tehtaqOverP)[i];

    cov(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundLoc0) =
        (*m_branches.m_track_cov_d0z0)[i];
    cov(Acts::BoundIndices::eBoundPhi, Acts::BoundIndices::eBoundLoc0) =
        (*m_branches.m_track_cov_d0phi)[i];
    cov(Acts::BoundIndices::eBoundTheta, Acts::BoundIndices::eBoundLoc0) =
        (*m_branches.m_track_cov_d0theta)[i];
    cov(Acts::BoundIndices::eBoundQOverP, Acts::BoundIndices::eBoundLoc0) =
        (*m_branches.m_track_cov_d0qOverP)[i];
    cov(Acts::BoundIndices::eBoundPhi, Acts::BoundIndices::eBoundLoc1) =
        (*m_branches.m_track_cov_z0phi)[i];
    cov(Acts::BoundIndices::eBoundTheta, Acts::BoundIndices::eBoundLoc1) =
        (*m_branches.m_track_cov_z0theta)[i];
    cov(Acts::BoundIndices::eBoundQOverP, Acts::BoundIndices::eBoundLoc1) =
        (*m_branches.m_track_cov_z0qOverP)[i];
    cov(Acts::BoundIndices::eBoundTheta, Acts::BoundIndices::eBoundPhi) =
        (*m_branches.m_track_cov_phitheta)[i];
    cov(Acts::BoundIndices::eBoundQOverP, Acts::BoundIndices::eBoundPhi) =
        (*m_branches.m_track_cov_phiqOverP)[i];
    cov(Acts::BoundIndices::eBoundQOverP, Acts::BoundIndices::eBoundTheta) =
        (*m_branches.m_track_cov_tehtaqOverP)[i];

    Acts::BoundTrackParameters tc(surface, params, q, cov);
    trackContainer.push_back(tc);
  }

  std::vector<Acts::Vector4> truthVertexContainer;
  for (unsigned int i = 0; i < nTruthVtx; i++) {
    Acts::Vector4 vtx(
        (*m_branches.m_truthvertex_x)[i], (*m_branches.m_truthvertex_y)[i],
        (*m_branches.m_truthvertex_z)[i], (*m_branches.m_truthvertex_t)[i]);
    truthVertexContainer.push_back(vtx);
  }
  std::vector<Acts::Vector4> recoVertexContainer;
  for (unsigned int i = 0; i < nRecoVtx; i++) {
    Acts::Vector4 vtx((*m_branches.m_recovertex_x)[i],
                      (*m_branches.m_recovertex_y)[i],
                      (*m_branches.m_recovertex_z)[i], 0);
    recoVertexContainer.push_back(vtx);
  }

  Acts::Vertex<Acts::BoundTrackParameters> beamspotConstraint;
  Acts::Vector3 beamspotPos;
  Acts::SymMatrix3 beamspotCov;

  beamspotPos << m_branches.m_beamspot_x, m_branches.m_beamspot_y,
      m_branches.m_beamspot_z;
  beamspotCov << m_branches.m_beamspot_sigX * m_branches.m_beamspot_sigX, 0, 0,
      0, m_branches.m_beamspot_sigY * m_branches.m_beamspot_sigY, 0, 0, 0,
      m_branches.m_beamspot_sigZ * m_branches.m_beamspot_sigZ;

  beamspotConstraint.setPosition(beamspotPos);
  beamspotConstraint.setCovariance(beamspotCov);

  context.eventStore.add(m_cfg.outputTrackParameters,
                         std::move(trackContainer));
  context.eventStore.add(m_cfg.outputTruthVtxParameters,
                         std::move(truthVertexContainer));
  context.eventStore.add(m_cfg.outputRecoVtxParameters,
                         std::move(recoVertexContainer));
  context.eventStore.add(m_cfg.outputBeamspotConstraint,
                         std::move(beamspotConstraint));

  // Return success flag
  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::RootNTupleReader::~RootNTupleReader() = default;
