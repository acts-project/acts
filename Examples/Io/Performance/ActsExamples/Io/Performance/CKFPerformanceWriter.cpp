// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"

#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <numeric>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>
#include <TVectorF.h>

ActsExamples::CKFPerformanceWriter::CKFPerformanceWriter(
    ActsExamples::CKFPerformanceWriter::Config cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.inputTrajectories, "CKFPerformanceWriter", lvl),
      m_cfg(std::move(cfg)),
      m_effPlotTool(m_cfg.effPlotToolConfig, lvl),
      m_fakeRatePlotTool(m_cfg.fakeRatePlotToolConfig, lvl),
      m_duplicationPlotTool(m_cfg.duplicationPlotToolConfig, lvl),
      m_trackSummaryPlotTool(m_cfg.trackSummaryPlotToolConfig, lvl) {
  // trajectories collection name is already checked by base ctor
  if (m_cfg.inputTrackParametersTips.empty()) {
    throw std::invalid_argument(
        "Missing track parameters tips input collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::invalid_argument("Could not open '" + m_cfg.filePath + "'");
  }

  // initialize the plot tools
  m_effPlotTool.book(m_effPlotCache);
  m_fakeRatePlotTool.book(m_fakeRatePlotCache);
  m_duplicationPlotTool.book(m_duplicationPlotCache);
  m_trackSummaryPlotTool.book(m_trackSummaryPlotCache);
}

ActsExamples::CKFPerformanceWriter::~CKFPerformanceWriter() {
  m_effPlotTool.clear(m_effPlotCache);
  m_fakeRatePlotTool.clear(m_fakeRatePlotCache);
  m_duplicationPlotTool.clear(m_duplicationPlotCache);
  m_trackSummaryPlotTool.clear(m_trackSummaryPlotCache);
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::CKFPerformanceWriter::endRun() {
  float eff_tracks = float(m_nTotalMatchedTracks) / m_nTotalTracks;
  float fakeRate_tracks = float(m_nTotalFakeTracks) / m_nTotalTracks;
  float duplicationRate_tracks =
      float(m_nTotalDuplicateTracks) / m_nTotalTracks;

  float eff_particle = float(m_nTotalMatchedParticles) / m_nTotalParticles;
  float fakeRate_particle = float(m_nTotalFakeParticles) / m_nTotalParticles;
  float duplicationRate_particle =
      float(m_nTotalDuplicateParticles) / m_nTotalParticles;

  ACTS_DEBUG("nTotalTracks                = " << m_nTotalTracks);
  ACTS_DEBUG("nTotalMatchedTracks         = " << m_nTotalMatchedTracks);
  ACTS_DEBUG("nTotalDuplicateTracks       = " << m_nTotalDuplicateTracks);
  ACTS_DEBUG("nTotalFakeTracks            = " << m_nTotalFakeTracks);

  ACTS_INFO(
      "Efficiency with tracks (nMatchedTracks/ nAllTracks) = " << eff_tracks);
  ACTS_INFO(
      "Fake rate with tracks (nFakeTracks/nAllTracks) = " << fakeRate_tracks);
  ACTS_INFO("Duplicate rate with tracks (nDuplicateTracks/nAllTracks) = "
            << duplicationRate_tracks);
  ACTS_INFO("Efficiency with particles (nMatchedParticles/nTrueParticles) = "
            << eff_particle);
  ACTS_INFO("Fake rate with particles (nFakeParticles/nTrueParticles) = "
            << fakeRate_particle);
  ACTS_INFO(
      "Duplicate rate with particles (nDuplicateParticles/nTrueParticles) = "
      << duplicationRate_particle);

  auto write_float = [&](float f, const char* name) {
    TVectorF v(1);
    v[0] = f;
    m_outputFile->WriteObject(&v, name);
  };

  if (m_outputFile != nullptr) {
    m_outputFile->cd();
    m_effPlotTool.write(m_effPlotCache);
    m_fakeRatePlotTool.write(m_fakeRatePlotCache);
    m_duplicationPlotTool.write(m_duplicationPlotCache);
    m_trackSummaryPlotTool.write(m_trackSummaryPlotCache);
    write_float(eff_tracks, "eff_tracks");
    write_float(fakeRate_tracks, "fakerate_tracks");
    write_float(duplicationRate_tracks, "duplicaterate_tracks");
    write_float(eff_particle, "eff_particles");
    write_float(fakeRate_particle, "fakerate_particles");
    write_float(duplicationRate_particle, "duplicaterate_particles");
    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::CKFPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const TrajectoriesContainer& trajectories) {
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;
  // The number of majority particle hits and fitted track parameters
  using RecoTrackInfo = std::pair<size_t, Acts::BoundTrackParameters>;
  using Acts::VectorHelpers::perp;

  const auto& trackTips =
      ctx.eventStore.get<std::vector<std::pair<size_t, size_t>>>(
          m_cfg.inputTrackParametersTips);
  // Read truth input collections
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);

  // Counter of truth-matched reco tracks
  std::map<ActsFatras::Barcode, std::vector<RecoTrackInfo>> matched;
  // Counter of truth-unmatched reco tracks
  std::map<ActsFatras::Barcode, size_t> unmatched;
  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCounts;

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Vector of input features for neural network classification
  std::vector<float> inputFeatures(3);

  // Loop over all trajectories
  for (auto [itraj, trackTip] : trackTips) {
    const auto& traj = trajectories[itraj];
    const auto& mj = traj.multiTrajectory();
    auto trajState =
        Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);

    // Reco track selection
    //@TODO: add interface for applying others cuts on reco tracks:
    // -> pT, d0, z0, detector-specific hits/holes number cut
    if (trajState.nMeasurements < m_cfg.nMeasurementsMin) {
      continue;
    }
    // Check if the reco track has fitted track parameters
    if (not traj.hasTrackParameters(trackTip)) {
      ACTS_WARNING(
          "No fitted track parameters for trajectory with entry index = "
          << trackTip);
      continue;
    }
    const auto& fittedParameters = traj.trackParameters(trackTip);
    // Requirement on the pT of the track
    const auto& momentum = fittedParameters.momentum();
    const auto pT = perp(momentum);
    if (pT < m_cfg.ptMin) {
      continue;
    }
    // Fill the trajectory summary info
    m_trackSummaryPlotTool.fill(m_trackSummaryPlotCache, fittedParameters,
                                trajState.nStates, trajState.nMeasurements,
                                trajState.nOutliers, trajState.nHoles,
                                trajState.nSharedHits);

    // Get the majority truth particle to this track
    identifyContributingParticles(hitParticlesMap, traj, trackTip,
                                  particleHitCounts);
    if (particleHitCounts.empty()) {
      ACTS_WARNING(
          "No truth particle associated with this trajectory with entry "
          "index = "
          << trackTip);
      continue;
    }
    // Get the majority particleId and majority particle counts
    // Note that the majority particle might be not in the truth seeds
    // collection
    ActsFatras::Barcode majorityParticleId =
        particleHitCounts.front().particleId;
    size_t nMajorityHits = particleHitCounts.front().hitCount;

    // Check if the trajectory is matched with truth.
    // If not, it will be classified as 'fake'
    bool isFake = false;
    if (nMajorityHits * 1. / trajState.nMeasurements >=
        m_cfg.truthMatchProbMin) {
      matched[majorityParticleId].push_back({nMajorityHits, fittedParameters});
    } else {
      isFake = true;
      unmatched[majorityParticleId]++;
    }
    // Fill fake rate plots
    m_fakeRatePlotTool.fill(m_fakeRatePlotCache, fittedParameters, isFake);

    // Use neural network classification for duplication rate plots
    // Currently, the network used for this example can only handle
    // good/duplicate classification, so need to manually exclude fake tracks
    if (m_cfg.duplicatedPredictor && !isFake) {
      inputFeatures[0] = trajState.nMeasurements;
      inputFeatures[1] = trajState.nOutliers;
      inputFeatures[2] = trajState.chi2Sum * 1.0 / trajState.NDF;
      // predict if current trajectory is 'duplicate'
      bool isDuplicated = m_cfg.duplicatedPredictor(inputFeatures);
      // Add to number of duplicated particles
      if (isDuplicated) {
        m_nTotalDuplicateTracks++;
      }
      // Fill the duplication rate
      m_duplicationPlotTool.fill(m_duplicationPlotCache, fittedParameters,
                                 isDuplicated);
    }
    // Counting number of total trajectories
    m_nTotalTracks++;
  }

  // Use truth-based classification for duplication rate plots
  if (!m_cfg.duplicatedPredictor) {
    // Loop over all truth-matched reco tracks for duplication rate plots
    for (auto& [particleId, matchedTracks] : matched) {
      // Sort the reco tracks matched to this particle by the number of majority
      // hits
      std::sort(matchedTracks.begin(), matchedTracks.end(),
                [](const RecoTrackInfo& lhs, const RecoTrackInfo& rhs) {
                  return lhs.first > rhs.first;
                });
      for (size_t itrack = 0; itrack < matchedTracks.size(); itrack++) {
        const auto& [nMajorityHits, fittedParameters] =
            matchedTracks.at(itrack);
        // The tracks with maximum number of majority hits is taken as the
        // 'real' track; others are as 'duplicated'
        bool isDuplicated = (itrack != 0);
        // the track is associated to the same particle
        if (isDuplicated) {
          m_nTotalDuplicateTracks++;
        }
        // Fill the duplication rate
        m_duplicationPlotTool.fill(m_duplicationPlotCache, fittedParameters,
                                   isDuplicated);
      }
    }
  }

  // Loop over all truth particle seeds for efficiency plots and reco details.
  // These are filled w.r.t. truth particle seed info
  for (const auto& particle : particles) {
    if (particle.transverseMomentum() < m_cfg.ptMin) {
      continue;
    }
    auto particleId = particle.particleId();
    // Investigate the truth-matched tracks
    size_t nMatchedTracks = 0;
    bool isReconstructed = false;
    auto imatched = matched.find(particleId);
    if (imatched != matched.end()) {
      nMatchedTracks = imatched->second.size();
      // Add number for total matched tracks here
      m_nTotalMatchedTracks += nMatchedTracks;
      m_nTotalMatchedParticles += 1;
      // Check if the particle has more than one matched track for the duplicate
      // rate
      if (nMatchedTracks > 1) {
        m_nTotalDuplicateParticles += 1;
      }
      isReconstructed = true;
    }
    // Fill efficiency plots
    m_effPlotTool.fill(m_effPlotCache, particle, isReconstructed);
    // Fill number of duplicated tracks for this particle
    m_duplicationPlotTool.fill(m_duplicationPlotCache, particle,
                               nMatchedTracks - 1);

    // Investigate the fake (i.e. truth-unmatched) tracks
    size_t nFakeTracks = 0;
    auto ifake = unmatched.find(particleId);
    if (ifake != unmatched.end()) {
      nFakeTracks = ifake->second;
      m_nTotalFakeTracks += nFakeTracks;
      // unmatched is a map of majority particle id to # of tracks associated
      // with that particle
      m_nTotalFakeParticles += 1;
    }
    // Fill number of reconstructed/truth-matched/fake tracks for this particle
    m_fakeRatePlotTool.fill(m_fakeRatePlotCache, particle, nMatchedTracks,
                            nFakeTracks);
    m_nTotalParticles += 1;
  }  // end all truth particles

  return ProcessCode::SUCCESS;
}
