// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <numeric>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Io/Performance/CKFPerformanceWriter.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"

FW::CKFPerformanceWriter::CKFPerformanceWriter(
    FW::CKFPerformanceWriter::Config cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.inputTrajectories, "CKFPerformanceWriter", lvl),
      m_cfg(std::move(cfg)),
      m_effPlotTool(m_cfg.effPlotToolConfig, lvl),
      m_fakeRatePlotTool(m_cfg.fakeRatePlotToolConfig, lvl),
      m_duplicationPlotTool(m_cfg.duplicationPlotToolConfig, lvl),
      m_trackSummaryPlotTool(m_cfg.trackSummaryPlotToolConfig, lvl) {
  // Input track and truth collection name
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing input trajectories collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles collection");
  }
  if (m_cfg.outputFilename.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  auto path = joinPaths(m_cfg.outputDir, m_cfg.outputFilename);
  m_outputFile = TFile::Open(path.c_str(), "RECREATE");
  if (not m_outputFile) {
    throw std::invalid_argument("Could not open '" + path + "'");
  }

  // initialize the plot tools
  m_effPlotTool.book(m_effPlotCache);
  m_fakeRatePlotTool.book(m_fakeRatePlotCache);
  m_duplicationPlotTool.book(m_duplicationPlotCache);
  m_trackSummaryPlotTool.book(m_trackSummaryPlotCache);
}

FW::CKFPerformanceWriter::~CKFPerformanceWriter() {
  m_effPlotTool.clear(m_effPlotCache);
  m_fakeRatePlotTool.clear(m_fakeRatePlotCache);
  m_duplicationPlotTool.clear(m_duplicationPlotCache);
  m_trackSummaryPlotTool.clear(m_trackSummaryPlotCache);
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

FW::ProcessCode FW::CKFPerformanceWriter::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_effPlotTool.write(m_effPlotCache);
    m_fakeRatePlotTool.write(m_fakeRatePlotCache);
    m_duplicationPlotTool.write(m_duplicationPlotCache);
    m_trackSummaryPlotTool.write(m_trackSummaryPlotCache);
    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode FW::CKFPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const TrajectoryContainer& trajectories) {
  // The number of majority particle hits and fitted track parameters
  using RecoTrackInfo = std::pair<size_t, Acts::BoundParameters>;

  // Read truth particles from input collection
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Counter of truth-matched reco tracks
  std::map<ActsFatras::Barcode, std::vector<RecoTrackInfo>> matched;
  // Counter of truth-unmatched reco tracks
  std::map<ActsFatras::Barcode, size_t> unmatched;

  // Loop over all trajectories
  for (const auto& traj : trajectories) {
    // The trajectory entry indices and the multiTrajectory
    const auto& [trackTips, mj] = traj.trajectory();
    if (trackTips.empty()) {
      ACTS_WARNING("Empty multiTrajectory.");
      continue;
    }

    // Loop over all trajectories in a multiTrajectory
    for (const size_t& trackTip : trackTips) {
      // Collect the trajectory summary info
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
      // Fill the trajectory summary info
      m_trackSummaryPlotTool.fill(m_trackSummaryPlotCache, fittedParameters,
                                  trajState.nStates, trajState.nMeasurements,
                                  trajState.nOutliers, trajState.nHoles);

      // Get the majority truth particle to this track
      std::vector<ParticleHitCount> particleHitCount =
          traj.identifyMajorityParticle(trackTip);
      if (particleHitCount.empty()) {
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
          particleHitCount.front().particleId;
      size_t nMajorityHits = particleHitCount.front().hitCount;

      // Check if the trajectory is matched with truth.
      // If not, it will be classified as 'fake'
      bool isFake = false;
      if (nMajorityHits * 1. / trajState.nMeasurements >=
          m_cfg.truthMatchProbMin) {
        matched[majorityParticleId].push_back(
            {nMajorityHits, fittedParameters});
      } else {
        isFake = true;
        unmatched[majorityParticleId]++;
      }
      // Fill fake rate plots
      m_fakeRatePlotTool.fill(m_fakeRatePlotCache, fittedParameters, isFake);
    }  // end all trajectories in a multiTrajectory
  }    // end all multiTrajectories

  // Loop over all truth-matched reco tracks for duplication rate plots
  for (auto& [particleId, matchedTracks] : matched) {
    // Sort the reco tracks matched to this particle by the number of majority
    // hits
    std::sort(matchedTracks.begin(), matchedTracks.end(),
              [](const RecoTrackInfo& lhs, const RecoTrackInfo& rhs) {
                return lhs.first > rhs.first;
              });
    for (size_t itrack = 0; itrack < matchedTracks.size(); itrack++) {
      const auto& [nMajorityHits, fittedParameters] = matchedTracks.at(itrack);
      // The tracks with maximum number of majority hits is taken as the 'real'
      // track; others are as 'duplicated'
      bool isDuplicated = (itrack != 0);
      // Fill the duplication rate
      m_duplicationPlotTool.fill(m_duplicationPlotCache, fittedParameters,
                                 isDuplicated);
    }
  }

  // Loop over all truth particle seeds for efficiency plots and reco details.
  // These are filled w.r.t. truth particle seed info
  for (const auto& particle : particles) {
    auto particleId = particle.particleId();
    // Investigate the truth-matched tracks
    size_t nMatchedTracks = 0;
    bool isReconstructed = false;
    auto imatched = matched.find(particleId);
    if (imatched != matched.end()) {
      nMatchedTracks = imatched->second.size();
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
    }
    // Fill number of reconstructed/truth-matched/fake tracks for this particle
    m_fakeRatePlotTool.fill(m_fakeRatePlotCache, particle, nMatchedTracks,
                            nFakeTracks);
  }  // end all truth particles

  return ProcessCode::SUCCESS;
}
