// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"

#include <format>

namespace ActsExamples {

TrackSummaryPlotTool::TrackSummaryPlotTool(
    const TrackSummaryPlotTool::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("TrackSummaryPlotTool", lvl)) {}

void TrackSummaryPlotTool::book(Cache& cache, const std::string& prefix) const {
  const auto& etaAxis = m_cfg.varBinning.at("Eta");
  const auto& ptAxis = m_cfg.varBinning.at("Pt");
  ACTS_DEBUG("Initialize the histograms for track info plots, use prefix '"
             << prefix << "'");
  auto addPrefix = [&](const std::string& name) {
    return prefix.empty() ? name : std::format("{}_{}", prefix, name);
  };

  // number of track states versus eta
  cache.nStates_vs_eta = Acts::Experimental::ProfileHistogram1(
      addPrefix("nStates_vs_eta"), "Number of total states vs. #eta",
      std::array{etaAxis}, "N");

  // number of measurements versus eta
  cache.nMeasurements_vs_eta = Acts::Experimental::ProfileHistogram1(
      addPrefix("nMeasurements_vs_eta"), "Number of measurements vs. #eta",
      std::array{etaAxis}, "N");

  // number of holes versus eta
  cache.nHoles_vs_eta = Acts::Experimental::ProfileHistogram1(
      addPrefix("nHoles_vs_eta"), "Number of holes vs. #eta",
      std::array{etaAxis}, "N");

  // number of outliers versus eta
  cache.nOutliers_vs_eta = Acts::Experimental::ProfileHistogram1(
      addPrefix("nOutliers_vs_eta"), "Number of outliers vs. #eta",
      std::array{etaAxis}, "N");

  // number of Shared Hits versus eta
  cache.nSharedHits_vs_eta = Acts::Experimental::ProfileHistogram1(
      addPrefix("nSharedHits_vs_eta"), "Number of Shared Hits vs. #eta",
      std::array{etaAxis}, "N");

  // number of track states versus pt
  cache.nStates_vs_pt = Acts::Experimental::ProfileHistogram1(
      addPrefix("nStates_vs_pT"), "Number of total states vs. pT",
      std::array{ptAxis}, "N");

  // number of measurements versus pt
  cache.nMeasurements_vs_pt = Acts::Experimental::ProfileHistogram1(
      addPrefix("nMeasurements_vs_pT"), "Number of measurements vs. pT",
      std::array{ptAxis}, "N");

  // number of holes versus pt
  cache.nHoles_vs_pt = Acts::Experimental::ProfileHistogram1(
      addPrefix("nHoles_vs_pT"), "Number of holes vs. pT", std::array{ptAxis},
      "N");

  // number of outliers versus pt
  cache.nOutliers_vs_pt = Acts::Experimental::ProfileHistogram1(
      addPrefix("nOutliers_vs_pT"), "Number of outliers vs. pT",
      std::array{ptAxis}, "N");

  // number of Shared Hits versus pt
  cache.nSharedHits_vs_pt = Acts::Experimental::ProfileHistogram1(
      addPrefix("nSharedHits_vs_pT"), "Number of Shared Hits vs. pT",
      std::array{ptAxis}, "N");
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

  cache.nStates_vs_eta.fill({fit_eta}, static_cast<double>(nStates));
  cache.nMeasurements_vs_eta.fill({fit_eta},
                                  static_cast<double>(nMeasurements));
  cache.nOutliers_vs_eta.fill({fit_eta}, static_cast<double>(nOutliers));
  cache.nHoles_vs_eta.fill({fit_eta}, static_cast<double>(nHoles));
  cache.nSharedHits_vs_eta.fill({fit_eta}, static_cast<double>(nSharedHits));

  cache.nStates_vs_pt.fill({fit_pT}, static_cast<double>(nStates));
  cache.nMeasurements_vs_pt.fill({fit_pT}, static_cast<double>(nMeasurements));
  cache.nOutliers_vs_pt.fill({fit_pT}, static_cast<double>(nOutliers));
  cache.nHoles_vs_pt.fill({fit_pT}, static_cast<double>(nHoles));
  cache.nSharedHits_vs_pt.fill({fit_pT}, static_cast<double>(nSharedHits));
}

}  // namespace ActsExamples
