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
      m_logger(Acts::getDefaultLogger("TrackCompletenessPlotTool", lvl)) {}

void TrackQualityPlotTool::book(Cache& cache) const {
  const auto& ptAxis = m_cfg.varBinning.at("Pt");
  const auto& etaAxis = m_cfg.varBinning.at("Eta");
  const auto& phiAxis = m_cfg.varBinning.at("Phi");
  ACTS_DEBUG("Initialize the histograms for completeness plots");

  // completeness vs pT
  cache.completeness_vs_pT.emplace("completeness_vs_pT",
                                   "Completeness;pT [GeV/c];Completeness",
                                   std::array{ptAxis}, "Completeness");

  // completeness vs eta
  cache.completeness_vs_eta.emplace("completeness_vs_eta",
                                    "Completeness;#eta;Completeness",
                                    std::array{etaAxis}, "Completeness");

  // completeness vs phi
  cache.completeness_vs_phi.emplace("completeness_vs_phi",
                                    "Completeness;#phi;Completeness",
                                    std::array{phiAxis}, "Completeness");

  // purity vs pT
  cache.purity_vs_pT.emplace("purity_vs_pT", "Purity;pT [GeV/c];Purity",
                             std::array{ptAxis}, "Purity");

  // purity vs eta
  cache.purity_vs_eta.emplace("purity_vs_eta", "Purity;#eta;Purity",
                              std::array{etaAxis}, "Purity");

  // purity vs phi
  cache.purity_vs_phi.emplace("purity_vs_phi", "Purity;#phi;Purity",
                              std::array{phiAxis}, "Purity");
}

void TrackQualityPlotTool::fill(
    Cache& cache, const Acts::BoundTrackParameters& fittedParameters,
    double completeness, double purity) const {
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  cache.completeness_vs_pT->fill({fit_pT}, completeness);
  cache.completeness_vs_eta->fill({fit_eta}, completeness);
  cache.completeness_vs_phi->fill({fit_phi}, completeness);

  cache.purity_vs_pT->fill({fit_pT}, purity);
  cache.purity_vs_eta->fill({fit_eta}, purity);
  cache.purity_vs_phi->fill({fit_phi}, purity);
}

}  // namespace ActsExamples
