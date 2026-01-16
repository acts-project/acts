// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/TrackQualityPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace ActsExamples {

TrackQualityPlotTool::TrackQualityPlotTool(const Config& cfg,
                                           Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("TrackQualityPlotTool", lvl)),
      m_completenessVsPt("completeness_vs_pT",
                         "Completeness;pT [GeV/c];Completeness",
                         std::array{m_cfg.varBinning.at("Pt")}, "Completeness"),
      m_completenessVsEta(
          "completeness_vs_eta", "Completeness;#eta;Completeness",
          std::array{m_cfg.varBinning.at("Eta")}, "Completeness"),
      m_completenessVsPhi(
          "completeness_vs_phi", "Completeness;#phi;Completeness",
          std::array{m_cfg.varBinning.at("Phi")}, "Completeness"),
      m_purityVsPt("purity_vs_pT", "Purity;pT [GeV/c];Purity",
                   std::array{m_cfg.varBinning.at("Pt")}, "Purity"),
      m_purityVsEta("purity_vs_eta", "Purity;#eta;Purity",
                    std::array{m_cfg.varBinning.at("Eta")}, "Purity"),
      m_purityVsPhi("purity_vs_phi", "Purity;#phi;Purity",
                    std::array{m_cfg.varBinning.at("Phi")}, "Purity") {
  ACTS_DEBUG("Initialize the histograms for track quality plots");
}

void TrackQualityPlotTool::fill(
    const Acts::BoundTrackParameters& fittedParameters, double completeness,
    double purity) {
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  m_completenessVsPt.fill({fit_pT}, completeness);
  m_completenessVsEta.fill({fit_eta}, completeness);
  m_completenessVsPhi.fill({fit_phi}, completeness);

  m_purityVsPt.fill({fit_pT}, purity);
  m_purityVsEta.fill({fit_eta}, purity);
  m_purityVsPhi.fill({fit_phi}, purity);
}

}  // namespace ActsExamples
