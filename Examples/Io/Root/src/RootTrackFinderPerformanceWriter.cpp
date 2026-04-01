// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootTrackFinderPerformanceWriter.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsPlugins/Root/HistogramConverter.hpp"

#include <stdexcept>

#include <TEfficiency.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TTree.h>
#include <TVectorFfwd.h>
#include <TVectorT.h>

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
      m_collector(
          TrackFinderPerformanceCollector::Config{
              m_cfg.effPlotToolConfig, m_cfg.fakePlotToolConfig,
              m_cfg.duplicationPlotToolConfig, m_cfg.trackSummaryPlotToolConfig,
              m_cfg.trackQualityPlotToolConfig,
              m_cfg.subDetectorTrackSummaryVolumes},
          logger().clone()) {
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
}

RootTrackFinderPerformanceWriter::~RootTrackFinderPerformanceWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode RootTrackFinderPerformanceWriter::finalize() {
  m_collector.logSummary();

  auto writeFloat = [&](float f, const char* name) {
    TVectorF v(1);
    v[0] = f;
    m_outputFile->WriteObject(&v, name);
  };

  if (m_outputFile != nullptr) {
    m_outputFile->cd();

    // Write efficiency histograms
    for (const auto& [name, eff] : m_collector.effPlotTool().efficiencies1D()) {
      toRoot(eff)->Write();
    }
    for (const auto& [name, eff] : m_collector.effPlotTool().efficiencies2D()) {
      toRoot(eff)->Write();
    }
    for (const auto& eff :
         m_collector.effPlotTool().trackEffVsEtaInPtRanges()) {
      toRoot(eff)->Write();
    }
    for (const auto& eff :
         m_collector.effPlotTool().trackEffVsPtInAbsEtaRanges()) {
      toRoot(eff)->Write();
    }

    // Write fake ratio histograms
    for (const auto& [name, hist] : m_collector.fakePlotTool().histograms()) {
      toRoot(hist)->Write();
    }
    for (const auto& [name, eff] : m_collector.fakePlotTool().efficiencies()) {
      toRoot(eff)->Write();
    }

    // Write duplication ratio histograms
    for (const auto& [name, prof] :
         m_collector.duplicationPlotTool().profiles()) {
      toRoot(prof)->Write();
    }
    for (const auto& [name, eff] :
         m_collector.duplicationPlotTool().efficiencies()) {
      toRoot(eff)->Write();
    }

    // Write track summary histograms
    writeTrackSummaryPlots(m_collector.trackSummaryPlotTool());
    for (const auto& [key, tool] : m_collector.subDetectorSummaryTools()) {
      writeTrackSummaryPlots(tool);
    }

    // Write track quality histograms
    for (const auto& [name, prof] :
         m_collector.trackQualityPlotTool().profiles()) {
      toRoot(prof)->Write();
    }

    // Write summary scalars derived from the collector's accumulated counts.
    const auto s = m_collector.stats();
    float eff_tracks =
        static_cast<float>(s.nTotalMatchedTracks) / s.nTotalTracks;
    float fakeRatio_tracks =
        static_cast<float>(s.nTotalFakeTracks) / s.nTotalTracks;
    float duplicationRatio_tracks =
        static_cast<float>(s.nTotalDuplicateTracks) / s.nTotalTracks;
    float eff_particle =
        static_cast<float>(s.nTotalMatchedParticles) / s.nTotalParticles;
    float fakeRatio_particle =
        static_cast<float>(s.nTotalFakeParticles) / s.nTotalParticles;
    float duplicationRatio_particle =
        static_cast<float>(s.nTotalDuplicateParticles) / s.nTotalParticles;

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
  // Read truth input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);
  const auto& particleTrackMatching = m_inputParticleTrackMatching(ctx);
  const auto& particleMeasurementsMap = m_inputParticleMeasurementsMap(ctx);

  // Exclusive access to the histograms while filling
  std::lock_guard<std::mutex> lock(m_writeMutex);

  m_collector.fill(ctx, tracks, particles, trackParticleMatching,
                   particleTrackMatching, particleMeasurementsMap);

  // Write additional matching details to TTree
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
