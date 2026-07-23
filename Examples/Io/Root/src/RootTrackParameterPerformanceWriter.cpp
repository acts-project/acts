// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootTrackParameterPerformanceWriter.hpp"

#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <stdexcept>
#include <utility>

#include <TFile.h>

namespace ActsExamples {

RootTrackParameterPerformanceWriter::RootTrackParameterPerformanceWriter(
    RootTrackParameterPerformanceWriter::Config config,
    Acts::Logging::Level level)
    : WriterT(config.inputTracks, "RootTrackParameterPerformanceWriter", level),
      m_cfg(std::move(config)),
      m_collector(
          TrackParameterPerformanceCollector::Config{m_cfg.resPlotToolConfig,
                                                     m_cfg.parameterType,
                                                     m_cfg.geometrySelection},
          logger().clone()) {
  // tracks collection name is already checked by base ctor
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.inputTrackParticleMatching.empty()) {
    throw std::invalid_argument("Missing input track particles matching");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument("Missing input measurement simulated hits map");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), "RECREATE");
  if (m_outputFile == nullptr) {
    throw std::invalid_argument("Could not open '" + path + "'");
  }
}

RootTrackParameterPerformanceWriter::~RootTrackParameterPerformanceWriter() {
  delete m_outputFile;
}

ProcessCode RootTrackParameterPerformanceWriter::finalize() {
  if (m_outputFile == nullptr) {
    return ProcessCode::SUCCESS;
  }

  m_collector.logSummary();

  m_outputFile->cd();

  writeResPlots(m_collector.resPlotTool(), m_cfg.resPlotRefinement, logger());

  ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");

  m_outputFile->Close();
  return ProcessCode::SUCCESS;
}

ProcessCode RootTrackParameterPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const ConstTrackContainer& tracks) {
  // Read truth input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);
  const auto& simHits = m_inputSimHits(ctx);
  const auto& measurementSimHitsMap = m_inputMeasurementSimHitsMap(ctx);

  // Exclusive access to the histograms while filling
  std::lock_guard<std::mutex> lock(m_writeMutex);

  m_collector.fill(ctx.geoContext, tracks, particles, trackParticleMatching,
                   simHits, measurementSimHitsMap);

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
