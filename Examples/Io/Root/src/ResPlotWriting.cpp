// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/ResPlotWriting.hpp"

#include "ActsPlugins/Root/HistogramConverter.hpp"

#include <string>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

namespace ActsExamples {

void writeResPlots(const ResPlotTool& resPlotTool,
                   const ResPlotRefinementConfig& config,
                   const Acts::Logger& logger) {
  using ActsPlugins::toRoot;

  // Helper lambda to write 2D histogram and extract mean/width profiles
  const auto writeWithRefinement = [&](auto& hist,
                                       const std::string& meanPrefix,
                                       const std::string& widthPrefix) {
    hist.Write();

    // Get the histogram name and extract the suffix (e.g., "_d0_vs_eta")
    const std::string baseName = hist.GetName();
    const std::string suffix = baseName.substr(baseName.find('_'));

    auto [meanHist, widthHist, fitFailureFraction] =
        ActsPlugins::extractMeanWidthProfiles(
            hist, meanPrefix + suffix, widthPrefix + suffix,
            config.fitMinEntries, config.fitSigmaRange, config.fitIterations,
            logger);
    if (fitFailureFraction >= config.warningThresholdFitFailureFraction) {
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
}

}  // namespace ActsExamples
