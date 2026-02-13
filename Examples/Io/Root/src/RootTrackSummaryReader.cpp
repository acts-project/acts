// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootTrackSummaryReader.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/Root/RootUtility.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <iostream>
#include <stdexcept>

#include <TChain.h>

namespace ActsExamples {

RootTrackSummaryReader::RootTrackSummaryReader(
    const RootTrackSummaryReader::Config& config, Acts::Logging::Level level)
    : IReader(),
      m_logger{Acts::getDefaultLogger(name(), level)},
      m_cfg(config) {
  m_inputChain = std::make_unique<TChain>(m_cfg.treeName.c_str());

  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing input filename");
  }

  m_outputTrackParameters.initialize(m_cfg.outputTracks);
  m_outputParticles.initialize(m_cfg.outputParticles);

  // Set the branches
  m_inputChain->SetBranchAddress("event_nr", &m_eventNr);
  m_inputChain->SetBranchAddress("multiTraj_nr", &m_multiTrajNr.get());
  m_inputChain->SetBranchAddress("subTraj_nr", &m_subTrajNr.get());

  // These info is not really stored in the event store, but still read in
  m_inputChain->SetBranchAddress("nStates", &m_nStates.get());
  m_inputChain->SetBranchAddress("nMeasurements", &m_nMeasurements.get());
  m_inputChain->SetBranchAddress("nOutliers", &m_nOutliers.get());
  m_inputChain->SetBranchAddress("nHoles", &m_nHoles.get());
  m_inputChain->SetBranchAddress("chi2Sum", &m_chi2Sum.get());
  m_inputChain->SetBranchAddress("NDF", &m_NDF.get());
  m_inputChain->SetBranchAddress("measurementChi2", &m_measurementChi2.get());
  m_inputChain->SetBranchAddress("outlierChi2", &m_outlierChi2.get());
  m_inputChain->SetBranchAddress("measurementVolume",
                                 &m_measurementVolume.get());
  m_inputChain->SetBranchAddress("measurementLayer", &m_measurementLayer.get());
  m_inputChain->SetBranchAddress("outlierVolume", &m_outlierVolume.get());
  m_inputChain->SetBranchAddress("outlierLayer", &m_outlierLayer.get());

  if (m_inputChain->GetBranch("majorityParticleId") != nullptr) {
    m_hasCombinedMajorityParticleId = true;
    m_majorityParticleId.allocate();
    m_inputChain->SetBranchAddress("majorityParticleId",
                                   &m_majorityParticleId.get());
  } else {
    m_hasCombinedMajorityParticleId = false;
    m_majorityParticleVertexPrimary.allocate();
    m_majorityParticleVertexSecondary.allocate();
    m_majorityParticleParticle.allocate();
    m_majorityParticleGeneration.allocate();
    m_majorityParticleSubParticle.allocate();
    m_inputChain->SetBranchAddress("majorityParticleId_vertex_primary",
                                   &m_majorityParticleVertexPrimary.get());
    m_inputChain->SetBranchAddress("majorityParticleId_vertex_secondary",
                                   &m_majorityParticleVertexSecondary.get());
    m_inputChain->SetBranchAddress("majorityParticleId_particle",
                                   &m_majorityParticleParticle.get());
    m_inputChain->SetBranchAddress("majorityParticleId_generation",
                                   &m_majorityParticleGeneration.get());
    m_inputChain->SetBranchAddress("majorityParticleId_sub_particle",
                                   &m_majorityParticleSubParticle.get());
  }
  m_inputChain->SetBranchAddress("nMajorityHits", &m_nMajorityHits.get());
  m_inputChain->SetBranchAddress("t_charge", &m_t_charge.get());
  m_inputChain->SetBranchAddress("t_time", &m_t_time.get());
  m_inputChain->SetBranchAddress("t_vx", &m_t_vx.get());
  m_inputChain->SetBranchAddress("t_vy", &m_t_vy.get());
  m_inputChain->SetBranchAddress("t_vz", &m_t_vz.get());
  m_inputChain->SetBranchAddress("t_px", &m_t_px.get());
  m_inputChain->SetBranchAddress("t_py", &m_t_py.get());
  m_inputChain->SetBranchAddress("t_pz", &m_t_pz.get());
  m_inputChain->SetBranchAddress("t_theta", &m_t_theta.get());
  m_inputChain->SetBranchAddress("t_phi", &m_t_phi.get());
  m_inputChain->SetBranchAddress("t_eta", &m_t_eta.get());
  m_inputChain->SetBranchAddress("t_pT", &m_t_pT.get());

  m_inputChain->SetBranchAddress("hasFittedParams", &m_hasFittedParams.get());
  m_inputChain->SetBranchAddress("eLOC0_fit", &m_eLOC0_fit.get());
  m_inputChain->SetBranchAddress("eLOC1_fit", &m_eLOC1_fit.get());
  m_inputChain->SetBranchAddress("ePHI_fit", &m_ePHI_fit.get());
  m_inputChain->SetBranchAddress("eTHETA_fit", &m_eTHETA_fit.get());
  m_inputChain->SetBranchAddress("eQOP_fit", &m_eQOP_fit.get());
  m_inputChain->SetBranchAddress("eT_fit", &m_eT_fit.get());
  m_inputChain->SetBranchAddress("err_eLOC0_fit", &m_err_eLOC0_fit.get());
  m_inputChain->SetBranchAddress("err_eLOC1_fit", &m_err_eLOC1_fit.get());
  m_inputChain->SetBranchAddress("err_ePHI_fit", &m_err_ePHI_fit.get());
  m_inputChain->SetBranchAddress("err_eTHETA_fit", &m_err_eTHETA_fit.get());
  m_inputChain->SetBranchAddress("err_eQOP_fit", &m_err_eQOP_fit.get());
  m_inputChain->SetBranchAddress("err_eT_fit", &m_err_eT_fit.get());

  auto path = m_cfg.filePath;

  // add file to the input chain
  m_inputChain->Add(path.c_str());
  ACTS_DEBUG("Adding File " << path << " to tree '" << m_cfg.treeName << "'.");

  m_events = m_inputChain->GetEntries();
  ACTS_DEBUG("The full chain has " << m_events << " entries.");

  // Sort the entry numbers of the events
  {
    // necessary to guarantee that m_inputChain->GetV1() is valid for the
    // entire range
    m_inputChain->SetEstimate(m_events + 1);

    m_entryNumbers.resize(m_events);
    m_inputChain->Draw("event_nr", "", "goff");
    RootUtility::stableSort(m_inputChain->GetEntries(), m_inputChain->GetV1(),
                            m_entryNumbers.data(), false);
  }
}

std::pair<std::size_t, std::size_t> RootTrackSummaryReader::availableEvents()
    const {
  return {0u, m_events};
}

ProcessCode RootTrackSummaryReader::read(const AlgorithmContext& context) {
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
    TrackParametersContainer trackParameterCollection;
    SimParticleContainer truthParticleCollection;

    // Read the correct entry
    auto entry = m_entryNumbers.at(context.eventNumber);
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
      Acts::BoundMatrix covMat;

      covMat << resD0 * resD0, 0., 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0.,
          0., 0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., 0., resTh * resTh,
          0., 0., 0., 0., 0., 0., resQp * resQp, 0., 0., 0., 0., 0., 0.,
          resT * resT;

      // TODO we do not have a hypothesis at hand here. defaulting to pion
      trackParameterCollection.push_back(Acts::BoundTrackParameters(
          perigeeSurface, paramVec, std::move(covMat),
          Acts::ParticleHypothesis::pion()));
    }

    unsigned int nTruthParticles = m_t_vx->size();
    for (unsigned int i = 0; i < nTruthParticles; i++) {
      SimParticleState truthParticle;

      truthParticle.setPosition4((*m_t_vx)[i], (*m_t_vy)[i], (*m_t_vz)[i],
                                 (*m_t_time)[i]);
      truthParticle.setDirection((*m_t_px)[i], (*m_t_py)[i], (*m_t_pz)[i]);
      auto makeBarcode = [](std::uint32_t vp, std::uint32_t vs,
                            std::uint32_t particle, std::uint32_t generation,
                            std::uint32_t subParticle) {
        SimBarcode barcode = SimBarcode::Invalid();
        return barcode
            .withVertexPrimary(static_cast<SimBarcode::PrimaryVertexId>(vp))
            .withVertexSecondary(static_cast<SimBarcode::SecondaryVertexId>(vs))
            .withParticle(static_cast<SimBarcode::ParticleId>(particle))
            .withGeneration(static_cast<SimBarcode::GenerationId>(generation))
            .withSubParticle(
                static_cast<SimBarcode::SubParticleId>(subParticle));
      };

      SimBarcode barcode = SimBarcode::Invalid();
      if (m_hasCombinedMajorityParticleId && m_majorityParticleId.hasValue()) {
        const auto& components = (*m_majorityParticleId)[i];
        auto comp = [&](std::size_t idx) -> std::uint32_t {
          return (components.size() > idx) ? components[idx] : 0u;
        };
        barcode = makeBarcode(comp(0), comp(1), comp(2), comp(3), comp(4));
      } else {
        auto safeAt = [](const auto& branch, std::size_t idx) -> std::uint32_t {
          const auto* vec = branch.get();
          return (vec != nullptr && vec->size() > idx) ? vec->at(idx) : 0u;
        };
        barcode = makeBarcode(safeAt(m_majorityParticleVertexPrimary, i),
                              safeAt(m_majorityParticleVertexSecondary, i),
                              safeAt(m_majorityParticleParticle, i),
                              safeAt(m_majorityParticleGeneration, i),
                              safeAt(m_majorityParticleSubParticle, i));
      }
      truthParticle.setParticleId(barcode);

      truthParticleCollection.insert(truthParticleCollection.end(),
                                     SimParticle(truthParticle, truthParticle));
    }
    // Write the collections to the EventStore
    m_outputTrackParameters(context, std::move(trackParameterCollection));
    m_outputParticles(context, std::move(truthParticleCollection));
  } else {
    ACTS_WARNING("Could not read in event.");
  }
  // Return success flag
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
