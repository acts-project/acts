// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootEstimatedParametersWriter.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <ios>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

#include "detail/AverageSimHits.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

ActsExamples::RootEstimatedParametersWriter::RootEstimatedParametersWriter(
    const ActsExamples::RootEstimatedParametersWriter::Config& cfg,
    Acts::Logging::Level lvl)
    : WriterT(cfg.inputTrackParameters, "RootEstimatedParametersWriter", lvl),
      m_cfg(cfg),
      m_outputFile(cfg.rootFile) {
  if (m_cfg.inputTrackParamsSeedMap.empty()) {
    throw std::invalid_argument(
        "Missing parameters-to-seed map input collection");
  }
  if (m_cfg.inputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds input collection");
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
  if (m_cfg.outputFilename.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.outputTreename.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    auto path = joinPaths(m_cfg.outputDir, m_cfg.outputFilename);
    m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + path);
    }
  }
  m_outputFile->cd();
  m_outputTree =
      new TTree(m_cfg.outputTreename.c_str(), m_cfg.outputTreename.c_str());
  if (m_outputTree == nullptr)
    throw std::bad_alloc();
  else {
    // I/O parameters
    m_outputTree->Branch("eventNr", &m_eventNr);
    m_outputTree->Branch("t_loc0", &m_t_loc0);
    m_outputTree->Branch("t_loc1", &m_t_loc1);
    m_outputTree->Branch("t_phi", &m_t_phi);
    m_outputTree->Branch("t_theta", &m_t_theta);
    m_outputTree->Branch("t_qop", &m_t_qop);
    m_outputTree->Branch("t_time", &m_t_time);
    m_outputTree->Branch("truthMatched", &m_truthMatched);

    m_outputTree->Branch("loc0_est", &m_loc0_est);
    m_outputTree->Branch("loc1_est", &m_loc1_est);
    m_outputTree->Branch("phi_est", &m_phi_est);
    m_outputTree->Branch("theta_est", &m_theta_est);
    m_outputTree->Branch("qop_est", &m_qop_est);
    m_outputTree->Branch("time_est", &m_time_est);
    m_outputTree->Branch("pt_est", &m_pt_est);
    m_outputTree->Branch("p_est", &m_p_est);
    m_outputTree->Branch("eta_est", &m_eta_est);
  }
}

ActsExamples::RootEstimatedParametersWriter::~RootEstimatedParametersWriter() {
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode
ActsExamples::RootEstimatedParametersWriter::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_outputTree->Write();
    ACTS_INFO("Write estimated parameters from seed to tree '"
              << m_cfg.outputTreename << "' in '"
              << joinPaths(m_cfg.outputDir, m_cfg.outputFilename) << "'");
  }
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootEstimatedParametersWriter::writeT(
    const AlgorithmContext& ctx, const TrackParametersContainer& parameters) {
  using TrackParamsSeedMap = std::map<Index, std::pair<Index, Index>>;
  using SeedContainer = std::vector<std::vector<Acts::Seed<SimSpacePoint>>>;
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;
  using HitSimHitsMap = IndexMultimap<Index>;

  if (m_outputFile == nullptr)
    return ProcessCode::SUCCESS;

  // Read additional input collections
  const auto& trackParamsSeedMap =
      ctx.eventStore.get<TrackParamsSeedMap>(m_cfg.inputTrackParamsSeedMap);
  const auto& seeds = ctx.eventStore.get<SeedContainer>(m_cfg.inputSeeds);
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimHits);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);
  const auto& hitSimHitsMap =
      ctx.eventStore.get<HitSimHitsMap>(m_cfg.inputMeasurementSimHitsMap);

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = ctx.eventNumber;

  // Loop over the estimated track parameters
  for (size_t iparams = 0; iparams < parameters.size(); ++iparams) {
    // The reference surface of the parameters, i.e. also the reference surface
    // of the first space point
    const auto& surface = parameters[iparams].referenceSurface();
    // The estimated bound parameters vector
    const auto params = parameters[iparams].parameters();
    m_loc0_est = params[Acts::eBoundLoc0];
    m_loc1_est = params[Acts::eBoundLoc1];
    m_phi_est = params[Acts::eBoundPhi];
    m_theta_est = params[Acts::eBoundTheta];
    m_qop_est = params[Acts::eBoundQOverP];
    m_time_est = params[Acts::eBoundTime];
    m_p_est = std::abs(1.0 / m_qop_est);
    m_pt_est = m_p_est * std::sin(m_theta_est);
    m_eta_est = std::atanh(std::cos(m_theta_est));

    // Get the seed via the estimated parameters to seed map
    const auto& [iregion, iseed] = trackParamsSeedMap.at(iparams);
    const auto& seed = seeds[iregion][iseed];
    // check the seed quality
    ProtoTrack ptrack{seed.sp()[0]->measurementIndex(),
                      seed.sp()[1]->measurementIndex(),
                      seed.sp()[2]->measurementIndex()};
    std::vector<ParticleHitCount> particleHitCounts;
    identifyContributingParticles(hitParticlesMap, ptrack, particleHitCounts);
    m_truthMatched = false;
    if (particleHitCounts.size() == 1) {
      m_truthMatched = true;
    }

    // Get the index of the first space point
    const auto& firstSP = seed.sp().front();
    const auto& hitIdx = firstSP->measurementIndex();
    // Get the sim hits via the measurement to sim hits map
    auto indices = makeRange(hitSimHitsMap.equal_range(hitIdx));
    auto [truthLocal, truthPos4, truthUnitDir] =
        detail::averageSimHits(ctx.geoContext, surface, simHits, indices);
    // Get the truth track parameter at the first space point
    m_t_loc0 = truthLocal[Acts::ePos0];
    m_t_loc1 = truthLocal[Acts::ePos1];
    m_t_phi = phi(truthUnitDir);
    m_t_theta = theta(truthUnitDir);
    m_t_time = truthPos4[Acts::eTime];
    // momemtum averaging makes even less sense than averaging position and
    // direction. use the first momentum or set q/p to zero
    if (not indices.empty()) {
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
        m_t_charge = particle.charge();
        m_t_qop = m_t_charge / p;
      } else {
        ACTS_WARNING("Truth particle with barcode = " << particleId
                                                      << " not found!");
      }
    }

    m_outputTree->Fill();
  }

  return ProcessCode::SUCCESS;
}
