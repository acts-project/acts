// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootTrackParameterWriter.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/MultiIndex.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <cmath>
#include <cstddef>
#include <ios>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

ActsExamples::RootTrackParameterWriter::RootTrackParameterWriter(
    const ActsExamples::RootTrackParameterWriter::Config& config,
    Acts::Logging::Level level)
    : TrackParameterWriter(config.inputTrackParameters,
                           "RootTrackParameterWriter", level),
      m_cfg(config) {
  if (m_cfg.inputProtoTracks.empty()) {
    throw std::invalid_argument("Missing proto tracks input collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-simulated-hits map input collection");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    auto path = m_cfg.filePath;
    m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + path + "'");
    }
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  } else {
    // The estimated track parameters
    m_outputTree->Branch("event_nr", &m_eventNr);
    m_outputTree->Branch("loc0", &m_loc0);
    m_outputTree->Branch("loc1", &m_loc1);
    m_outputTree->Branch("phi", &m_phi);
    m_outputTree->Branch("theta", &m_theta);
    m_outputTree->Branch("qop", &m_qop);
    m_outputTree->Branch("time", &m_time);
    m_outputTree->Branch("p", &m_p);
    m_outputTree->Branch("pt", &m_pt);
    m_outputTree->Branch("eta", &m_eta);
    // The truth track parameters
    m_outputTree->Branch("eventNr", &m_eventNr);
    m_outputTree->Branch("t_loc0", &m_t_loc0);
    m_outputTree->Branch("t_loc1", &m_t_loc1);
    m_outputTree->Branch("t_phi", &m_t_phi);
    m_outputTree->Branch("t_theta", &m_t_theta);
    m_outputTree->Branch("t_qop", &m_t_qop);
    m_outputTree->Branch("t_time", &m_t_time);
    m_outputTree->Branch("truthMatched", &m_truthMatched);
  }
}

ActsExamples::RootTrackParameterWriter::~RootTrackParameterWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootTrackParameterWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  ACTS_INFO("Wrote estimated parameters from seed to tree '"
            << m_cfg.treeName << "' in '" << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootTrackParameterWriter::writeT(
    const ActsExamples::AlgorithmContext& ctx,
    const TrackParametersContainer& trackParams) {
  // Read additional input collections
  const auto& protoTracks = m_inputProtoTracks(ctx);
  const auto& particles = m_inputParticles(ctx);
  const auto& simHits = m_inputSimHits(ctx);
  const auto& hitParticlesMap = m_inputMeasurementParticlesMap(ctx);
  const auto& hitSimHitsMap = m_inputMeasurementSimHitsMap(ctx);

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = ctx.eventNumber;

  ACTS_VERBOSE("Writing " << trackParams.size() << " track parameters");

  // Loop over the estimated track parameters
  for (std::size_t iparams = 0; iparams < trackParams.size(); ++iparams) {
    // The reference surface of the parameters, i.e. also the reference surface
    // of the first space point
    const auto& surface = trackParams[iparams].referenceSurface();
    // The estimated bound parameters vector
    const auto params = trackParams[iparams].parameters();
    m_loc0 = params[Acts::eBoundLoc0];
    m_loc1 = params[Acts::eBoundLoc1];
    m_phi = params[Acts::eBoundPhi];
    m_theta = params[Acts::eBoundTheta];
    m_qop = params[Acts::eBoundQOverP];
    m_time = params[Acts::eBoundTime];
    m_p = std::abs(1.0 / m_qop);
    m_pt = m_p * std::sin(m_theta);
    m_eta = std::atanh(std::cos(m_theta));

    // Get the proto track from which the track parameters are estimated
    const auto& ptrack = protoTracks[iparams];
    std::vector<ParticleHitCount> particleHitCounts;
    identifyContributingParticles(hitParticlesMap, ptrack, particleHitCounts);
    m_truthMatched = false;
    if (particleHitCounts.size() == 1) {
      m_truthMatched = true;
    }
    // Get the index of the first space point
    const auto& hitIdx = ptrack.front();
    // Get the sim hits via the measurement to sim hits map
    auto indices = makeRange(hitSimHitsMap.equal_range(hitIdx));
    auto [truthLocal, truthPos4, truthUnitDir] =
        averageSimHits(ctx.geoContext, surface, simHits, indices, logger());
    // Get the truth track parameter at the first space point
    m_t_loc0 = truthLocal[Acts::ePos0];
    m_t_loc1 = truthLocal[Acts::ePos1];
    m_t_phi = phi(truthUnitDir);
    m_t_theta = theta(truthUnitDir);
    m_t_time = truthPos4[Acts::eTime];
    // momentum averaging makes even less sense than averaging position and
    // direction. use the first momentum or set q/p to zero
    if (!indices.empty()) {
      // we assume that the indices are within valid ranges so we do not
      // need to check their validity again.
      const auto simHitIdx0 = indices.begin()->second;
      const auto& simHit0 = *simHits.nth(simHitIdx0);
      const auto p =
          simHit0.momentum4Before().template segment<3>(Acts::eMom0).norm();
      const auto& particleId = simHit0.particleId();
      // The truth charge has to be retrieved from the sim particle
      auto ip = particles.find(particleId);
      if (ip != particles.end()) {
        const auto& particle = *ip;
        m_t_charge = static_cast<int>(particle.charge());
        m_t_qop = m_t_charge / p;
      } else {
        ACTS_DEBUG("Truth particle with barcode "
                   << particleId << "=" << particleId.value() << " not found!");
      }
    }

    m_outputTree->Fill();
  }

  return ProcessCode::SUCCESS;
}
