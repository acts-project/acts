// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootTrackFinderPerformanceWriter.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsPlugins/Root/HistogramConverter.hpp"

#include <cstddef>
#include <map>
#include <ostream>
#include <stdexcept>
#include <utility>

#include <TEfficiency.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TTree.h>
#include <TVectorFfwd.h>
#include <TVectorT.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::phi;
using ActsPlugins::toRoot;

namespace ActsExamples {

namespace {

void writeTrackSummaryPlots(const TrackSummaryPlotTool& tool) {
  for (const auto& [name, prof] : tool.profiles()) {
    toRoot(prof)->Write();
  }
}

}  // namespace

RootTrackFinderPerformanceWriter::RootTrackFinderPerformanceWriter(
    RootTrackFinderPerformanceWriter::Config cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.inputTracks, "RootTrackFinderPerformanceWriter", lvl),
      m_cfg(std::move(cfg)),
      m_effPlotTool(m_cfg.effPlotToolConfig, lvl),
      m_fakePlotTool(m_cfg.fakePlotToolConfig, lvl),
      m_duplicationPlotTool(m_cfg.duplicationPlotToolConfig, lvl),
      m_trackSummaryPlotTool(m_cfg.trackSummaryPlotToolConfig, lvl),
      m_trackQualityPlotTool(m_cfg.trackQualityPlotToolConfig, lvl) {
  // tracks collection name is already checked by base ctor
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.inputTrackParticleMatching.empty()) {
    throw std::invalid_argument("Missing input track particles matching");
  }
  if (m_cfg.inputParticleTrackMatching.empty()) {
    throw std::invalid_argument("Missing input particle track matching");
  }
  if (m_cfg.inputParticleMeasurementsMap.empty()) {
    throw std::invalid_argument("Missing input measurement particles map");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);
  m_inputParticleTrackMatching.initialize(m_cfg.inputParticleTrackMatching);
  m_inputParticleMeasurementsMap.initialize(m_cfg.inputParticleMeasurementsMap);

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::invalid_argument("Could not open '" + m_cfg.filePath + "'");
  }

  if (m_cfg.writeMatchingDetails) {
    m_matchingTree = new TTree("matchingdetails", "matchingdetails");

    m_matchingTree->Branch("event_nr", &m_treeEventNr);
    m_matchingTree->Branch("particle_id_vertex_primary",
                           &m_treeParticleVertexPrimary);
    m_matchingTree->Branch("particle_id_vertex_secondary",
                           &m_treeParticleVertexSecondary);
    m_matchingTree->Branch("particle_id_particle", &m_treeParticleParticle);
    m_matchingTree->Branch("particle_id_generation", &m_treeParticleGeneration);
    m_matchingTree->Branch("particle_id_sub_particle",
                           &m_treeParticleSubParticle);
    m_matchingTree->Branch("matched", &m_treeIsMatched);
  }

  // Create subdetector track summary tools with prefixes
  for (const auto& [key, _] : m_cfg.subDetectorTrackSummaryVolumes) {
    TrackSummaryPlotTool::Config subConfig = m_cfg.trackSummaryPlotToolConfig;
    subConfig.prefix = key;
    m_subDetectorSummaryTools.emplace(std::piecewise_construct,
                                      std::forward_as_tuple(key),
                                      std::forward_as_tuple(subConfig, lvl));
  }
}

RootTrackFinderPerformanceWriter::~RootTrackFinderPerformanceWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode RootTrackFinderPerformanceWriter::finalize() {
  float eff_tracks = static_cast<float>(m_nTotalMatchedTracks) / m_nTotalTracks;
  float fakeRatio_tracks =
      static_cast<float>(m_nTotalFakeTracks) / m_nTotalTracks;
  float duplicationRatio_tracks =
      static_cast<float>(m_nTotalDuplicateTracks) / m_nTotalTracks;

  float eff_particle =
      static_cast<float>(m_nTotalMatchedParticles) / m_nTotalParticles;
  float fakeRatio_particle =
      static_cast<float>(m_nTotalFakeParticles) / m_nTotalParticles;
  float duplicationRatio_particle =
      static_cast<float>(m_nTotalDuplicateParticles) / m_nTotalParticles;

  ACTS_DEBUG("nTotalTracks                = " << m_nTotalTracks);
  ACTS_DEBUG("nTotalMatchedTracks         = " << m_nTotalMatchedTracks);
  ACTS_DEBUG("nTotalDuplicateTracks       = " << m_nTotalDuplicateTracks);
  ACTS_DEBUG("nTotalFakeTracks            = " << m_nTotalFakeTracks);

  ACTS_INFO(
      "Efficiency with tracks (nMatchedTracks/ nAllTracks) = " << eff_tracks);
  ACTS_INFO(
      "Fake ratio with tracks (nFakeTracks/nAllTracks) = " << fakeRatio_tracks);
  ACTS_INFO("Duplicate ratio with tracks (nDuplicateTracks/nAllTracks) = "
            << duplicationRatio_tracks);
  ACTS_INFO("Efficiency with particles (nMatchedParticles/nTrueParticles) = "
            << eff_particle);
  ACTS_INFO("Fake ratio with particles (nFakeParticles/nTrueParticles) = "
            << fakeRatio_particle);
  ACTS_INFO(
      "Duplicate ratio with particles (nDuplicateParticles/nTrueParticles) = "
      << duplicationRatio_particle);

  auto writeFloat = [&](float f, const char* name) {
    TVectorF v(1);
    v[0] = f;
    m_outputFile->WriteObject(&v, name);
  };

  if (m_outputFile != nullptr) {
    m_outputFile->cd();

    // Write efficiency histograms
    for (const auto& [name, eff] : m_effPlotTool.efficiencies1D()) {
      toRoot(eff)->Write();
    }
    for (const auto& [name, eff] : m_effPlotTool.efficiencies2D()) {
      toRoot(eff)->Write();
    }

    for (const auto& eff : m_effPlotTool.trackEffVsEtaInPtRanges()) {
      toRoot(eff)->Write();
    }
    for (const auto& eff : m_effPlotTool.trackEffVsPtInAbsEtaRanges()) {
      toRoot(eff)->Write();
    }

    // Write fake ratio histograms
    for (const auto& [name, hist] : m_fakePlotTool.histograms()) {
      toRoot(hist)->Write();
    }
    for (const auto& [name, eff] : m_fakePlotTool.efficiencies()) {
      toRoot(eff)->Write();
    }

    // Write duplication ratio histograms
    for (const auto& [name, prof] : m_duplicationPlotTool.profiles()) {
      toRoot(prof)->Write();
    }
    for (const auto& [name, eff] : m_duplicationPlotTool.efficiencies()) {
      toRoot(eff)->Write();
    }

    // Write track summary histograms
    writeTrackSummaryPlots(m_trackSummaryPlotTool);
    for (const auto& [key, tool] : m_subDetectorSummaryTools) {
      writeTrackSummaryPlots(tool);
    }

    // Write track quality histograms
    for (const auto& [name, prof] : m_trackQualityPlotTool.profiles()) {
      toRoot(prof)->Write();
    }

    writeFloat(eff_tracks, "eff_tracks");
    writeFloat(fakeRatio_tracks, "fakeratio_tracks");
    writeFloat(duplicationRatio_tracks, "duplicateratio_tracks");
    writeFloat(eff_particle, "eff_particles");
    writeFloat(fakeRatio_particle, "fakeratio_particles");
    writeFloat(duplicationRatio_particle, "duplicateratio_particles");

    if (m_matchingTree != nullptr) {
      m_matchingTree->Write();
    }

    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootTrackFinderPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const ConstTrackContainer& tracks) {
  // The number of majority particle hits and fitted track parameters
  using Acts::VectorHelpers::perp;

  // Read truth input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);
  const auto& particleTrackMatching = m_inputParticleTrackMatching(ctx);
  const auto& particleMeasurementsMap = m_inputParticleMeasurementsMap(ctx);

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Vector of input features for neural network classification
  std::vector<float> inputFeatures(3);

  ACTS_DEBUG("Collect information from " << tracks.size() << " tracks");
  std::size_t unmatched = 0, missingRefSurface = 0;
  for (const auto& track : tracks) {
    // Counting number of total trajectories
    m_nTotalTracks++;

    // Check if the reco track has fitted track parameters
    if (!track.hasReferenceSurface()) {
      ACTS_VERBOSE("No fitted track parameters for track, index = "
                   << track.index() << " tip index = " << track.tipIndex());
      missingRefSurface++;
      continue;
    }

    Acts::BoundTrackParameters fittedParameters =
        track.createParametersAtReference();

    // Fill the trajectory summary info
    m_trackSummaryPlotTool.fill(fittedParameters, track.nTrackStates(),
                                track.nMeasurements(), track.nOutliers(),
                                track.nHoles(), track.nSharedHits());

    // Potentially fill other track summary caches for the given volumes
    for (const auto& [key, volumes] : m_cfg.subDetectorTrackSummaryVolumes) {
      ACTS_VERBOSE("Fill track summary stats for subset " << key);
      std::size_t nTrackStates{}, nMeasurements{}, nOutliers{}, nHoles{},
          nSharedHits{};
      for (auto state : track.trackStatesReversed()) {
        if (!state.hasReferenceSurface() ||
            !volumes.contains(state.referenceSurface().geometryId().volume())) {
          continue;
        }

        nTrackStates++;
        nMeasurements +=
            static_cast<std::size_t>(state.typeFlags().isMeasurement());
        nOutliers += static_cast<std::size_t>(state.typeFlags().isOutlier());
        nHoles += static_cast<std::size_t>(state.typeFlags().isHole());
        nSharedHits +=
            static_cast<std::size_t>(state.typeFlags().isSharedHit());
      }
      m_subDetectorSummaryTools.at(key).fill(fittedParameters, nTrackStates,
                                             nMeasurements, nOutliers, nHoles,
                                             nSharedHits);
    }

    // Get the truth matching information
    auto imatched = trackParticleMatching.find(track.index());
    if (imatched == trackParticleMatching.end()) {
      ACTS_DEBUG("No truth matching information for this track, index = "
                 << track.index() << " tip index = " << track.tipIndex());
      unmatched++;
      continue;
    }

    const auto& particleMatch = imatched->second;

    if (particleMatch.classification == TrackMatchClassification::Fake) {
      m_nTotalFakeTracks++;
    }

    if (particleMatch.classification == TrackMatchClassification::Duplicate) {
      m_nTotalDuplicateTracks++;
    }

    // Fill fake ratio plots
    m_fakePlotTool.fill(fittedParameters, particleMatch.classification ==
                                              TrackMatchClassification::Fake);

    // Fill the duplication ratio
    m_duplicationPlotTool.fill(
        fittedParameters,
        particleMatch.classification == TrackMatchClassification::Duplicate);

    if (particleMatch.particle.has_value() &&
        particleMeasurementsMap.contains(particleMatch.particle.value())) {
      const auto measurements =
          particleMeasurementsMap.equal_range(particleMatch.particle.value());

      std::size_t nTrackMeasurements =
          track.nMeasurements() + track.nOutliers();
      std::size_t nMatchedHits =
          particleMatch.contributingParticles.front().hitCount;
      std::size_t nParticleHits =
          std::distance(measurements.first, measurements.second);

      double completeness = static_cast<double>(nMatchedHits) / nParticleHits;
      double purity = static_cast<double>(nMatchedHits) / nTrackMeasurements;

      // Fill the track quality plots
      m_trackQualityPlotTool.fill(fittedParameters, completeness, purity);
    }
  }

  if (unmatched > 0) {
    ACTS_DEBUG("No matching information found for " << unmatched << " tracks");
  }
  if (missingRefSurface > 0) {
    ACTS_DEBUG("Reference surface was missing for " << missingRefSurface
                                                    << " tracks");
  }

  // Loop over all truth particles for efficiency plots and reco details.
  for (const auto& particle : particles) {
    auto particleId = particle.particleId();

    // Investigate the truth-matched tracks
    std::size_t nMatchedTracks = 0;
    std::size_t nFakeTracks = 0;
    bool isReconstructed = false;
    if (auto imatched = particleTrackMatching.find(particleId);
        imatched != particleTrackMatching.end()) {
      isReconstructed = imatched->second.track.has_value();
      nMatchedTracks = (isReconstructed ? 1 : 0) + imatched->second.duplicates;

      // Add number for total matched tracks here
      m_nTotalMatchedTracks += nMatchedTracks;
      m_nTotalMatchedParticles += isReconstructed ? 1 : 0;

      // Check if the particle has more than one matched track for the duplicate
      // rate/ratio
      if (nMatchedTracks > 1) {
        m_nTotalDuplicateParticles += 1;
      }

      nFakeTracks = imatched->second.fakes;
      if (nFakeTracks > 0) {
        m_nTotalFakeParticles += 1;
      }
    }

    // Loop over all the other truth particle and find the distance to the
    // closest one
    double minDeltaR = -1;
    for (const auto& closeParticle : particles) {
      if (closeParticle.particleId() == particleId) {
        continue;
      }
      double distance = Acts::VectorHelpers::deltaR(particle.direction(),
                                                    closeParticle.direction());
      if (minDeltaR == -1 || distance < minDeltaR) {
        minDeltaR = distance;
      }
    }

    // Fill efficiency plots
    m_effPlotTool.fill(ctx.geoContext, particle.initialState(), minDeltaR,
                       isReconstructed);
    // Fill number of duplicated tracks for this particle
    // Guard against underflow when nMatchedTracks == 0
    m_duplicationPlotTool.fill(particle.initialState(), nMatchedTracks);

    // Fill number of reconstructed/truth-matched/fake tracks for this particle
    m_fakePlotTool.fill(particle.initialState(), nMatchedTracks, nFakeTracks);

    m_nTotalParticles += 1;
  }

  // Write additional stuff to TTree
  if (m_cfg.writeMatchingDetails && m_matchingTree != nullptr) {
    for (const auto& particle : particles) {
      auto particleId = particle.particleId();

      m_treeEventNr = ctx.eventNumber;
      m_treeParticleVertexPrimary = particleId.vertexPrimary();
      m_treeParticleVertexSecondary = particleId.vertexSecondary();
      m_treeParticleParticle = particleId.particle();
      m_treeParticleGeneration = particleId.generation();
      m_treeParticleSubParticle = particleId.subParticle();

      m_treeIsMatched = false;
      if (auto imatched = particleTrackMatching.find(particleId);
          imatched != particleTrackMatching.end()) {
        m_treeIsMatched = imatched->second.track.has_value();
      }

      m_matchingTree->Fill();
    }
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
