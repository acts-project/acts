// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootTrackFitterPerformanceWriter.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsPlugins/Root/HistogramConverter.hpp"

#include <stdexcept>
#include <utility>

#include <TEfficiency.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>

using ActsPlugins::toRoot;

namespace ActsExamples {

RootTrackFitterPerformanceWriter::RootTrackFitterPerformanceWriter(
    RootTrackFitterPerformanceWriter::Config config, Acts::Logging::Level level)
    : WriterT(config.inputTracks, "RootTrackFitterPerformanceWriter", level),
      m_cfg(std::move(config)),
      m_collector(
          TrackFitterPerformanceCollector::Config{
              m_cfg.resPlotToolConfig, m_cfg.effPlotToolConfig,
              m_cfg.trackSummaryPlotToolConfig, m_cfg.fitMinEntries,
              m_cfg.fitSigmaRange, m_cfg.fitIterations},
          logger().clone()) {
  // trajectories collection name is already checked by base ctor
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.inputTrackParticleMatching.empty()) {
    throw std::invalid_argument("Missing input track particles matching");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), "RECREATE");
  if (m_outputFile == nullptr) {
    throw std::invalid_argument("Could not open '" + path + "'");
  }
}

RootTrackFitterPerformanceWriter::~RootTrackFitterPerformanceWriter() {
  delete m_outputFile;
}

ProcessCode RootTrackFitterPerformanceWriter::finalize() {
  if (m_outputFile == nullptr) {
    return ProcessCode::SUCCESS;
  }

  m_collector.logSummary();

  m_outputFile->cd();

  const auto& resPlotTool = m_collector.resPlotTool();
  const auto& effPlotTool = m_collector.effPlotTool();
  const auto& trackSummaryPlotTool = m_collector.trackSummaryPlotTool();

  // Helper lambda to write 2D histogram and extract mean/width profiles
  const auto writeWithRefinement = [this](auto& hist,
                                          const std::string& meanPrefix,
                                          const std::string& widthPrefix) {
    hist.Write();

    // Get the histogram name and extract the suffix (e.g., "_d0_vs_eta")
    const std::string baseName = hist.GetName();
    const std::string suffix = baseName.substr(baseName.find('_'));

    auto [meanHist, widthHist, fitFailureFraction] =
        ActsPlugins::extractMeanWidthProfiles(
            hist, meanPrefix + suffix, widthPrefix + suffix,
            m_cfg.fitMinEntries, m_cfg.fitSigmaRange, m_cfg.fitIterations,
            logger());
    if (fitFailureFraction >= m_cfg.warningThresholdFitFailureFraction) {
      ACTS_WARNING("Fit failures for " << baseName << ": "
                                       << fitFailureFraction * 100 << "%");
    }

    meanHist->Write();
    widthHist->Write();
  };

  // Write residual histograms
  for (const auto& [name, hist] : resPlotTool.res()) {
    toRoot(hist)->Write();
  }
  for (const auto& [name, hist] : resPlotTool.resVsEta()) {
    writeWithRefinement(*toRoot(hist), "resmean", "reswidth");
  }
  for (const auto& [name, hist] : resPlotTool.resVsPt()) {
    writeWithRefinement(*toRoot(hist), "resmean", "reswidth");
  }
  for (const auto& [name, hist] : resPlotTool.resVsEtaPhi()) {
    writeWithRefinement(*toRoot(hist), "resmean", "reswidth");
  }
  for (const auto& [name, hist] : resPlotTool.resVsEtaPt()) {
    writeWithRefinement(*toRoot(hist), "resmean", "reswidth");
  }

  // Write pull histograms
  for (const auto& [name, hist] : resPlotTool.pull()) {
    toRoot(hist)->Write();
  }
  for (const auto& [name, hist] : resPlotTool.pullVsEta()) {
    writeWithRefinement(*toRoot(hist), "pullmean", "pullwidth");
  }
  for (const auto& [name, hist] : resPlotTool.pullVsPt()) {
    writeWithRefinement(*toRoot(hist), "pullmean", "pullwidth");
  }
  for (const auto& [name, hist] : resPlotTool.pullVsEtaPhi()) {
    writeWithRefinement(*toRoot(hist), "pullmean", "pullwidth");
  }
  for (const auto& [name, hist] : resPlotTool.pullVsEtaPt()) {
    writeWithRefinement(*toRoot(hist), "pullmean", "pullwidth");
  }

  // Write efficiency histograms
  for (const auto& [name, eff] : effPlotTool.efficiencies1D()) {
    toRoot(eff)->Write();
  }
  for (const auto& [name, eff] : effPlotTool.efficiencies2D()) {
    toRoot(eff)->Write();
  }

  for (const auto& eff : effPlotTool.trackEffVsEtaInPtRanges()) {
    toRoot(eff)->Write();
  }
  for (const auto& eff : effPlotTool.trackEffVsPtInAbsEtaRanges()) {
    toRoot(eff)->Write();
  }

  // Write track summary histograms
  for (const auto& [name, prof] : trackSummaryPlotTool.profiles()) {
    toRoot(prof)->Write();
  }

  ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");

  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode RootTrackFitterPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const ConstTrackContainer& tracks) {
  // Read truth input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);

  // Exclusive access to the histograms while filling
  std::lock_guard<std::mutex> lock(m_writeMutex);

  m_collector.fill(ctx.geoContext, tracks, particles, trackParticleMatching);

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
