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

using namespace Acts::Experimental;
using namespace ActsExamples;

namespace { 

ProfileHistogram1 makeProfile(
    const TrackSummaryPlotTool::TrackSummaryPlotTool::Config cfg,
    std::string name, const std::string& title,
    const Acts::Experimental::AxisVariant& ax) {
  if( !cfg.prefix.empty()) {
    name = std::format("{}_{}", cfg.prefix, name);
  }
  const auto &yAxis = cfg.varBinning.at("Num");
  Acts::Range1D<double> yRange{yAxis.bin(0).lower(), yAxis.bin(yAxis.size() - 1).upper()};
  return ProfileHistogram1(name, title, {ax}, yAxis.metadata(), yRange);
}

}

namespace ActsExamples {

TrackSummaryPlotTool::TrackSummaryPlotTool(
    const TrackSummaryPlotTool::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("TrackSummaryPlotTool", lvl)),
      m_nStatesVsEta(
          makeProfile(m_cfg, "nStates_vs_eta", "Number of total states vs. #eta",
                      m_cfg.varBinning.at("Eta"))),
      m_nMeasurementsVsEta(makeProfile(
          m_cfg, "nMeasurements_vs_eta", "Number of measurements vs. #eta",
          m_cfg.varBinning.at("Eta"))),
      m_nHolesVsEta(makeProfile(
          m_cfg, "nHoles_vs_eta", "Number of holes vs. #eta",
          m_cfg.varBinning.at("Eta"))),
      m_nOutliersVsEta(makeProfile(
          m_cfg, "nOutliers_vs_eta", "Number of outliers vs. #eta",
          m_cfg.varBinning.at("Eta"))),
      m_nSharedHitsVsEta(makeProfile(
          m_cfg, "nSharedHits_vs_eta", "Number of Shared Hits vs. #eta",
          m_cfg.varBinning.at("Eta"))),
      m_nStatesVsPt(makeProfile(
          m_cfg, "nStates_vs_pT", "Number of total states vs. pT",
          m_cfg.varBinning.at("Pt"))),
      m_nMeasurementsVsPt(makeProfile(
          m_cfg, "nMeasurements_vs_pT", "Number of measurements vs. pT",
          m_cfg.varBinning.at("Pt"))),
      m_nHolesVsPt(makeProfile(
          m_cfg, "nHoles_vs_pT", "Number of holes vs. pT",
          m_cfg.varBinning.at("Pt"))),
      m_nOutliersVsPt(makeProfile(
          m_cfg, "nOutliers_vs_pT", "Number of outliers vs. pT",
          m_cfg.varBinning.at("Pt"))),
      m_nSharedHitsVsPt(makeProfile(
          m_cfg, "nSharedHits_vs_pT", "Number of Shared Hits vs. pT",
          m_cfg.varBinning.at("Pt"))) {
    ACTS_DEBUG(
      "Initialize the histograms for track info plots"
      << (m_cfg.prefix.empty() ? "" : ", use prefix '" + m_cfg.prefix + "'"));
}

void TrackSummaryPlotTool::fill(
    const Acts::BoundTrackParameters& fittedParameters, std::size_t nStates,
    std::size_t nMeasurements, std::size_t nOutliers, std::size_t nHoles,
    std::size_t nSharedHits) {
  using Acts::VectorHelpers::eta;
  using Acts::VectorHelpers::perp;
  const auto momentum = fittedParameters.momentum();
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  m_nStatesVsEta.fill({fit_eta}, static_cast<double>(nStates));
  m_nMeasurementsVsEta.fill({fit_eta}, static_cast<double>(nMeasurements));
  m_nOutliersVsEta.fill({fit_eta}, static_cast<double>(nOutliers));
  m_nHolesVsEta.fill({fit_eta}, static_cast<double>(nHoles));
  m_nSharedHitsVsEta.fill({fit_eta}, static_cast<double>(nSharedHits));

  m_nStatesVsPt.fill({fit_pT}, static_cast<double>(nStates));
  m_nMeasurementsVsPt.fill({fit_pT}, static_cast<double>(nMeasurements));
  m_nOutliersVsPt.fill({fit_pT}, static_cast<double>(nOutliers));
  m_nHolesVsPt.fill({fit_pT}, static_cast<double>(nHoles));
  m_nSharedHitsVsPt.fill({fit_pT}, static_cast<double>(nSharedHits));
}

}  // namespace ActsExamples
