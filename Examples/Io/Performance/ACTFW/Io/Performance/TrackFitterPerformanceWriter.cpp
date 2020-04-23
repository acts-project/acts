// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Performance/TrackFitterPerformanceWriter.hpp"

#include <Acts/Utilities/Helpers.hpp>
#include <TFile.h>
#include <TTree.h>
#include <stdexcept>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Utilities/Paths.hpp"

using Acts::VectorHelpers::eta;

FW::TrackFitterPerformanceWriter::TrackFitterPerformanceWriter(
    FW::TrackFitterPerformanceWriter::Config cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.inputTrajectories, "TrackFitterPerformanceWriter", lvl),
      m_cfg(std::move(cfg)),
      m_resPlotTool(m_cfg.resPlotToolConfig, lvl),
      m_effPlotTool(m_cfg.effPlotToolConfig, lvl),
      m_trackSummaryPlotTool(m_cfg.trackSummaryPlotToolConfig, lvl)

{
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

  // initialize the residual and efficiency plots tool
  m_resPlotTool.book(m_resPlotCache);
  m_effPlotTool.book(m_effPlotCache);
  m_trackSummaryPlotTool.book(m_trackSummaryPlotCache);
}

FW::TrackFitterPerformanceWriter::~TrackFitterPerformanceWriter() {
  m_resPlotTool.clear(m_resPlotCache);
  m_effPlotTool.clear(m_effPlotCache);
  m_trackSummaryPlotTool.clear(m_trackSummaryPlotCache);

  if (m_outputFile) {
    m_outputFile->Close();
  }
}

FW::ProcessCode FW::TrackFitterPerformanceWriter::endRun() {
  // fill residual and pull details into additional hists
  m_resPlotTool.refinement(m_resPlotCache);

  if (m_outputFile) {
    m_outputFile->cd();
    m_resPlotTool.write(m_resPlotCache);
    m_effPlotTool.write(m_effPlotCache);
    m_trackSummaryPlotTool.write(m_trackSummaryPlotCache);

    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode FW::TrackFitterPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const TrajectoryContainer& trajectories) {
  // Read truth particles from input collection
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // All reconstructed trajectories with truth info
  std::map<ActsFatras::Barcode, TruthFitTrack> reconTrajectories;

  // Loop over all trajectories
  for (const auto& traj : trajectories) {
    if (not traj.hasTrajectory()) {
      continue;
    }
    const auto& [trackTip, track] = traj.trajectory();

    // get the majority truth particle to this track
    const auto particleHitCount = traj.identifyMajorityParticle();
    if (particleHitCount.empty()) {
      continue;
    }

    // find the truth particle for the majority barcode
    const auto ip = particles.find(particleHitCount.front().particleId);
    if (ip == particles.end()) {
      continue;
    }

    // record this trajectory with its truth info
    reconTrajectories.emplace(ip->particleId(), traj);

    // count the total number of hits and hits from the majority truth
    // particle
    size_t nTotalStates = 0, nHits = 0, nOutliers = 0, nHoles = 0;
    track.visitBackwards(trackTip, [&](const auto& state) {
      nTotalStates++;
      auto typeFlags = state.typeFlags();
      if (typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
        nHits++;
      } else if (typeFlags.test(Acts::TrackStateFlag::OutlierFlag)) {
        nOutliers++;
      } else if (typeFlags.test(Acts::TrackStateFlag::HoleFlag)) {
        nHoles++;
      }
    });

    // fill the track detailed info
    m_trackSummaryPlotTool.fill(m_trackSummaryPlotCache, *ip, nTotalStates,
                                nHits, nOutliers, nHoles);

    // fill the residual plots it the track has fitted parameter
    if (traj.hasTrackParameters()) {
      m_resPlotTool.fill(m_resPlotCache, ctx.geoContext, *ip,
                         traj.trackParameters());
    }
  }

  // Fill the efficiency, defined as the ratio between number of tracks with
  // fitted parameter and total truth tracks (assumes one truth partilce means
  // one truth track)
  // @Todo: add fake rate plots
  for (const auto& particle : particles) {
    const auto it = reconTrajectories.find(particle.particleId());
    if (it != reconTrajectories.end()) {
      // when the trajectory is reconstructed
      m_effPlotTool.fill(m_effPlotCache, particle,
                         it->second.hasTrackParameters());
    } else {
      // when the trajectory is NOT reconstructed
      m_effPlotTool.fill(m_effPlotCache, particle, false);
    }
  }

  return ProcessCode::SUCCESS;
}
