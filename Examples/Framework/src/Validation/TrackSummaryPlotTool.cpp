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
      m_logger(Acts::getDefaultLogger("TrackSummaryPlotTool", lvl)),
      m_nStatesVsEta(m_cfg.prefix.empty()
                         ? "nStates_vs_eta"
                         : std::format("{}_nStates_vs_eta", m_cfg.prefix),
                     "Number of total states vs. #eta",
                     std::array{m_cfg.varBinning.at("Eta")}, "N"),
      m_nMeasurementsVsEta(
          m_cfg.prefix.empty()
              ? "nMeasurements_vs_eta"
              : std::format("{}_nMeasurements_vs_eta", m_cfg.prefix),
          "Number of measurements vs. #eta",
          std::array{m_cfg.varBinning.at("Eta")}, "N"),
      m_nHolesVsEta(m_cfg.prefix.empty()
                        ? "nHoles_vs_eta"
                        : std::format("{}_nHoles_vs_eta", m_cfg.prefix),
                    "Number of holes vs. #eta",
                    std::array{m_cfg.varBinning.at("Eta")}, "N"),
      m_nOutliersVsEta(m_cfg.prefix.empty()
                           ? "nOutliers_vs_eta"
                           : std::format("{}_nOutliers_vs_eta", m_cfg.prefix),
                       "Number of outliers vs. #eta",
                       std::array{m_cfg.varBinning.at("Eta")}, "N"),
      m_nSharedHitsVsEta(
          m_cfg.prefix.empty()
              ? "nSharedHits_vs_eta"
              : std::format("{}_nSharedHits_vs_eta", m_cfg.prefix),
          "Number of Shared Hits vs. #eta",
          std::array{m_cfg.varBinning.at("Eta")}, "N"),
      m_nStatesVsPt(m_cfg.prefix.empty()
                        ? "nStates_vs_pT"
                        : std::format("{}_nStates_vs_pT", m_cfg.prefix),
                    "Number of total states vs. pT",
                    std::array{m_cfg.varBinning.at("Pt")}, "N"),
      m_nMeasurementsVsPt(
          m_cfg.prefix.empty()
              ? "nMeasurements_vs_pT"
              : std::format("{}_nMeasurements_vs_pT", m_cfg.prefix),
          "Number of measurements vs. pT",
          std::array{m_cfg.varBinning.at("Pt")}, "N"),
      m_nHolesVsPt(
          m_cfg.prefix.empty() ? "nHoles_vs_pT"
                               : std::format("{}_nHoles_vs_pT", m_cfg.prefix),
          "Number of holes vs. pT", std::array{m_cfg.varBinning.at("Pt")}, "N"),
      m_nOutliersVsPt(m_cfg.prefix.empty()
                          ? "nOutliers_vs_pT"
                          : std::format("{}_nOutliers_vs_pT", m_cfg.prefix),
                      "Number of outliers vs. pT",
                      std::array{m_cfg.varBinning.at("Pt")}, "N"),
      m_nSharedHitsVsPt(m_cfg.prefix.empty()
                            ? "nSharedHits_vs_pT"
                            : std::format("{}_nSharedHits_vs_pT", m_cfg.prefix),
                        "Number of Shared Hits vs. pT",
                        std::array{m_cfg.varBinning.at("Pt")}, "N") {
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
