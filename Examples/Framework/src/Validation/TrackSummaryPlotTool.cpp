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
    const TrackSummaryPlotTool::TrackSummaryPlotTool::Config& cfg,
    std::string name, const std::string& title,
    const Acts::Experimental::AxisVariant& ax) {
  if (!cfg.prefix.empty()) {
    name = std::format("{}_{}", cfg.prefix, name);
  }
  const auto& yAxis = cfg.varBinning.at("Num");
  Acts::Range1D<double> yRange{yAxis.bin(0).lower(),
                               yAxis.bin(yAxis.size() - 1).upper()};
  return ProfileHistogram1(name, title, {ax}, yAxis.metadata(), yRange);
}

}  // namespace

namespace ActsExamples {

TrackSummaryPlotTool::TrackSummaryPlotTool(
    const TrackSummaryPlotTool::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("TrackSummaryPlotTool", lvl)) {
  ACTS_DEBUG(
      "Initialize the histograms for track info plots"
      << (m_cfg.prefix.empty() ? "" : ", use prefix '" + m_cfg.prefix + "'"));

  m_profiles.insert(
      {"nStates_vs_eta",
       makeProfile(m_cfg, "nStates_vs_eta", "Number of total states vs. #eta",
                   m_cfg.varBinning.at("Eta"))});
  m_profiles.insert(
      {"nMeasurements_vs_eta", makeProfile(m_cfg, "nMeasurements_vs_eta",
                                           "Number of measurements vs. #eta",
                                           m_cfg.varBinning.at("Eta"))});
  m_profiles.insert({"nHoles_vs_eta", makeProfile(m_cfg, "nHoles_vs_eta",
                                                  "Number of holes vs. #eta",
                                                  m_cfg.varBinning.at("Eta"))});
  m_profiles.insert(
      {"nOutliers_vs_eta",
       makeProfile(m_cfg, "nOutliers_vs_eta", "Number of outliers vs. #eta",
                   m_cfg.varBinning.at("Eta"))});
  m_profiles.insert(
      {"nSharedHits_vs_eta", makeProfile(m_cfg, "nSharedHits_vs_eta",
                                         "Number of Shared Hits vs. #eta",
                                         m_cfg.varBinning.at("Eta"))});
  m_profiles.insert(
      {"nStates_vs_pT",
       makeProfile(m_cfg, "nStates_vs_pT", "Number of total states vs. pT",
                   m_cfg.varBinning.at("Pt"))});
  m_profiles.insert(
      {"nMeasurements_vs_pT", makeProfile(m_cfg, "nMeasurements_vs_pT",
                                          "Number of measurements vs. pT",
                                          m_cfg.varBinning.at("Pt"))});
  m_profiles.insert({"nHoles_vs_pT", makeProfile(m_cfg, "nHoles_vs_pT",
                                                 "Number of holes vs. pT",
                                                 m_cfg.varBinning.at("Pt"))});
  m_profiles.insert(
      {"nOutliers_vs_pT",
       makeProfile(m_cfg, "nOutliers_vs_pT", "Number of outliers vs. pT",
                   m_cfg.varBinning.at("Pt"))});
  m_profiles.insert(
      {"nSharedHits_vs_pT",
       makeProfile(m_cfg, "nSharedHits_vs_pT", "Number of Shared Hits vs. pT",
                   m_cfg.varBinning.at("Pt"))});
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

  m_profiles.at("nStates_vs_eta").fill({fit_eta}, static_cast<double>(nStates));
  m_profiles.at("nMeasurements_vs_eta")
      .fill({fit_eta}, static_cast<double>(nMeasurements));
  m_profiles.at("nOutliers_vs_eta")
      .fill({fit_eta}, static_cast<double>(nOutliers));
  m_profiles.at("nHoles_vs_eta").fill({fit_eta}, static_cast<double>(nHoles));
  m_profiles.at("nSharedHits_vs_eta")
      .fill({fit_eta}, static_cast<double>(nSharedHits));

  m_profiles.at("nStates_vs_pT").fill({fit_pT}, static_cast<double>(nStates));
  m_profiles.at("nMeasurements_vs_pT")
      .fill({fit_pT}, static_cast<double>(nMeasurements));
  m_profiles.at("nOutliers_vs_pT")
      .fill({fit_pT}, static_cast<double>(nOutliers));
  m_profiles.at("nHoles_vs_pT").fill({fit_pT}, static_cast<double>(nHoles));
  m_profiles.at("nSharedHits_vs_pT")
      .fill({fit_pT}, static_cast<double>(nSharedHits));
}

}  // namespace ActsExamples
