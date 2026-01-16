// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/DuplicationPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace ActsExamples {

DuplicationPlotTool::DuplicationPlotTool(const DuplicationPlotTool::Config& cfg,
                                         Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("DuplicationPlotTool", lvl)),
      m_nDuplicatedVsPt("nDuplicated_vs_pT",
                        "Number of duplicated track candidates",
                        std::array{m_cfg.varBinning.at("Pt")}, "N"),
      m_nDuplicatedVsEta("nDuplicated_vs_eta",
                         "Number of duplicated track candidates",
                         std::array{m_cfg.varBinning.at("Eta")}, "N"),
      m_nDuplicatedVsPhi("nDuplicated_vs_phi",
                         "Number of duplicated track candidates",
                         std::array{m_cfg.varBinning.at("Phi")}, "N"),
      m_duplicationRatioVsPt(
          "duplicationratio_vs_pT",
          "Duplication ratio;pT [GeV/c];Duplication ratio",
          std::array{m_cfg.varBinning.at("Pt")}),
      m_duplicationRatioVsEta("duplicationratio_vs_eta",
                              "Duplication ratio;#eta;Duplication ratio",
                              std::array{m_cfg.varBinning.at("Eta")}),
      m_duplicationRatioVsPhi("duplicationratio_vs_phi",
                              "Duplication ratio;#phi;Duplication ratio",
                              std::array{m_cfg.varBinning.at("Phi")}) {
  ACTS_DEBUG("Initialize the histograms for duplication ratio plots");
}

void DuplicationPlotTool::fill(
    const Acts::BoundTrackParameters& fittedParameters, bool status) {
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  m_duplicationRatioVsPt.fill({fit_pT}, status);
  m_duplicationRatioVsEta.fill({fit_eta}, status);
  m_duplicationRatioVsPhi.fill({fit_phi}, status);
}

void DuplicationPlotTool::fill(const SimParticleState& truthParticle,
                               std::size_t nDuplicatedTracks) {
  const auto t_phi = phi(truthParticle.direction());
  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();

  m_nDuplicatedVsPt.fill({t_pT}, static_cast<double>(nDuplicatedTracks));
  m_nDuplicatedVsEta.fill({t_eta}, static_cast<double>(nDuplicatedTracks));
  m_nDuplicatedVsPhi.fill({t_phi}, static_cast<double>(nDuplicatedTracks));
}

}  // namespace ActsExamples
