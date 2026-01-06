// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"

#include <TProfile.h>

namespace ActsExamples {

TrackSummaryPlotTool::TrackSummaryPlotTool(
    const TrackSummaryPlotTool::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("TrackSummaryPlotTool", lvl)) {}

void TrackSummaryPlotTool::book(Cache& cache, const std::string& prefix) const {
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  PlotHelpers::Binning bNum = m_cfg.varBinning.at("Num");
  ACTS_DEBUG("Initialize the histograms for track info plots, use prefix '"
             << prefix << "'");
  auto addPrefix = [&](const std::string& name) {
    return prefix.empty() ? name : prefix + "_" + name;
  };
  // number of track states versus eta
  cache.nStates_vs_eta =
      PlotHelpers::bookProf(addPrefix("nStates_vs_eta").c_str(),
                            "Number of total states vs. #eta", bEta, bNum);
  // number of measurements versus eta
  cache.nMeasurements_vs_eta =
      PlotHelpers::bookProf(addPrefix("nMeasurements_vs_eta").c_str(),
                            "Number of measurements vs. #eta", bEta, bNum);
  // number of holes versus eta
  cache.nHoles_vs_eta =
      PlotHelpers::bookProf(addPrefix("nHoles_vs_eta").c_str(),
                            "Number of holes vs. #eta", bEta, bNum);
  // number of outliers versus eta
  cache.nOutliers_vs_eta =
      PlotHelpers::bookProf(addPrefix("nOutliers_vs_eta").c_str(),
                            "Number of outliers vs. #eta", bEta, bNum);
  // number of Shared Hits versus eta
  cache.nSharedHits_vs_eta =
      PlotHelpers::bookProf(addPrefix("nSharedHits_vs_eta").c_str(),
                            "Number of Shared Hits vs. #eta", bEta, bNum);
  // number of track states versus pt
  cache.nStates_vs_pt =
      PlotHelpers::bookProf(addPrefix("nStates_vs_pT").c_str(),
                            "Number of total states vs. pT", bPt, bNum);
  // number of measurements versus pt
  cache.nMeasurements_vs_pt =
      PlotHelpers::bookProf(addPrefix("nMeasurements_vs_pT").c_str(),
                            "Number of measurements vs. pT", bPt, bNum);
  // number of holes versus pt
  cache.nHoles_vs_pt = PlotHelpers::bookProf(
      addPrefix("nHoles_vs_pT").c_str(), "Number of holes vs. pT", bPt, bNum);
  // number of outliers versus pt
  cache.nOutliers_vs_pt =
      PlotHelpers::bookProf(addPrefix("nOutliers_vs_pT").c_str(),
                            "Number of outliers vs. pT", bPt, bNum);
  // number of Shared Hits versus pt
  cache.nSharedHits_vs_pt =
      PlotHelpers::bookProf(addPrefix("nSharedHits_vs_pT").c_str(),
                            "Number of Shared Hits vs. pT", bPt, bNum);
}

void TrackSummaryPlotTool::clear(Cache& cache) const {
  delete cache.nStates_vs_eta;
  delete cache.nMeasurements_vs_eta;
  delete cache.nOutliers_vs_eta;
  delete cache.nHoles_vs_eta;
  delete cache.nSharedHits_vs_eta;
  delete cache.nStates_vs_pt;
  delete cache.nMeasurements_vs_pt;
  delete cache.nOutliers_vs_pt;
  delete cache.nHoles_vs_pt;
  delete cache.nSharedHits_vs_pt;
}

void TrackSummaryPlotTool::write(const Cache& cache) const {
  ACTS_DEBUG("Write the plots to output file.");
  cache.nStates_vs_eta->Write();
  cache.nMeasurements_vs_eta->Write();
  cache.nOutliers_vs_eta->Write();
  cache.nHoles_vs_eta->Write();
  cache.nSharedHits_vs_eta->Write();
  cache.nStates_vs_pt->Write();
  cache.nMeasurements_vs_pt->Write();
  cache.nOutliers_vs_pt->Write();
  cache.nHoles_vs_pt->Write();
  cache.nSharedHits_vs_pt->Write();
}

void TrackSummaryPlotTool::fill(
    Cache& cache, const Acts::BoundTrackParameters& fittedParameters,
    std::size_t nStates, std::size_t nMeasurements, std::size_t nOutliers,
    std::size_t nHoles, std::size_t nSharedHits) const {
  using Acts::VectorHelpers::eta;
  using Acts::VectorHelpers::perp;
  const auto momentum = fittedParameters.momentum();
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  PlotHelpers::fillProf(cache.nStates_vs_eta, fit_eta, nStates);
  PlotHelpers::fillProf(cache.nMeasurements_vs_eta, fit_eta, nMeasurements);
  PlotHelpers::fillProf(cache.nOutliers_vs_eta, fit_eta, nOutliers);
  PlotHelpers::fillProf(cache.nHoles_vs_eta, fit_eta, nHoles);
  PlotHelpers::fillProf(cache.nSharedHits_vs_eta, fit_eta, nSharedHits);

  PlotHelpers::fillProf(cache.nStates_vs_pt, fit_pT, nStates);
  PlotHelpers::fillProf(cache.nMeasurements_vs_pt, fit_pT, nMeasurements);
  PlotHelpers::fillProf(cache.nOutliers_vs_pt, fit_pT, nOutliers);
  PlotHelpers::fillProf(cache.nHoles_vs_pt, fit_pT, nHoles);
  PlotHelpers::fillProf(cache.nSharedHits_vs_pt, fit_pT, nSharedHits);
}

}  // namespace ActsExamples
