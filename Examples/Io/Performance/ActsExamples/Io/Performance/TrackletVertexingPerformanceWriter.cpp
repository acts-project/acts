// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TrackletVertexingPerformanceWriter.hpp"

#include "Acts/Utilities/MultiIndex.hpp"
#include "ActsExamples/Utilities/EventDataTransforms.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>
#include <string>
#include <TFile.h>

namespace ActsExamples {
struct AlgorithmContext;
}  // namespace ActsExamples

ActsExamples::TrackletVertexingPerformanceWriter::TrackletVertexingPerformanceWriter(
    ActsExamples::TrackletVertexingPerformanceWriter::Config config,
    Acts::Logging::Level level)
    : WriterT(config.inputSeeds, "TrackletVertexingPerformanceWriter", level),
      m_cfg(std::move(config)),
      m_trkVtxPlotTool(m_cfg.trkVtxPlotToolConfig, level) {
  if (m_cfg.inputRecPrimaryVertex.empty()) {
    throw std::invalid_argument("Missing input z rec");
  }
  if (m_cfg.inputGenPrimaryVertex.empty()) {
    throw std::invalid_argument("Missing input z gen");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  m_inputRecPrimaryVertex.initialize(m_cfg.inputRecPrimaryVertex);
  m_inputGenPrimaryVertex.initialize(m_cfg.inputGenPrimaryVertex);

  m_inputFitFunction.initialize(m_cfg.inputFitFunction);
  m_inputZTracklets.initialize(m_cfg.inputZTracklets);
  m_inputZTrackletsPeak.initialize(m_cfg.inputZTrackletsPeak);
  ev_counter = 0;
  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::invalid_argument("Could not open '" + path + "'");
  }

  // Create a subdirectory in the file
  subdirHist = m_outputFile->mkdir("all");
  subdirPeak = m_outputFile->mkdir("peak");

  // initialize the plot tools
  m_trkVtxPlotTool.book(m_trkVtxPlotCache);
}

ActsExamples::TrackletVertexingPerformanceWriter::~TrackletVertexingPerformanceWriter() {
  m_trkVtxPlotTool.clear(m_trkVtxPlotCache);
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::TrackletVertexingPerformanceWriter::finalize() {
  if (m_outputFile != nullptr) {
    m_outputFile->cd();
    m_trkVtxPlotTool.write(m_trkVtxPlotCache);
    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::TrackletVertexingPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const SimSeedContainer& seeds) {

  const auto& recPrimaryVertex = m_inputRecPrimaryVertex(ctx);
  const auto& genPrimaryVertex = m_inputGenPrimaryVertex(ctx);
  const auto& binZTracklets = m_inputZTracklets(ctx);
  const auto& binZTrackletsPeak = m_inputZTrackletsPeak(ctx);
  ev_counter += 1;
  m_trkVtxPlotTool.fill(m_trkVtxPlotCache, recPrimaryVertex, genPrimaryVertex); 
  if (m_outputFile != nullptr) {
    std::string histName = "histName" + std::to_string(ev_counter);
    std::string histPeakName = "histPeakName" + std::to_string(ev_counter);
    TH1D hist(histName.c_str(),";;",60,-65,10);
    TH1D histPeak(histPeakName.c_str(),";;",60,-65,10);
    for(int i=1;i<=60;i++){
      hist.SetBinContent(i, binZTracklets[i-1]);
      histPeak.SetBinContent(i, binZTrackletsPeak[i-1]);
    }
    subdirHist->cd();
    hist.Write();
    subdirPeak->cd();
    histPeak.Write();
  }
  return ProcessCode::SUCCESS;
}
