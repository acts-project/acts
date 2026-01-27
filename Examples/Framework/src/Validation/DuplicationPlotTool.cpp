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

using namespace Acts::Experimental;

namespace {

ProfileHistogram1 makeProfile(
    const ActsExamples::DuplicationPlotTool::Config& cfg, std::string name,
    const std::string& title, const AxisVariant& ax) {
  const auto& yAxis = cfg.varBinning.at("Num");
  Acts::Range1D<double> yRange{yAxis.bin(0).lower(),
                               yAxis.bin(yAxis.size() - 1).upper()};
  return ProfileHistogram1(name, title, {ax}, yAxis.metadata(), yRange);
}

}  // namespace

namespace ActsExamples {

DuplicationPlotTool::DuplicationPlotTool(const DuplicationPlotTool::Config& cfg,
                                         Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("DuplicationPlotTool", lvl)),
      m_nDuplicatedVsPt(makeProfile(m_cfg, "nDuplicated_vs_pT",
                                    "Number of duplicated track candidates",
                                    m_cfg.varBinning.at("Pt"))),
      m_nDuplicatedVsEta(makeProfile(m_cfg, "nDuplicated_vs_eta",
                                     "Number of duplicated track candidates",
                                     m_cfg.varBinning.at("Eta"))),
      m_nDuplicatedVsPhi(makeProfile(m_cfg, "nDuplicated_vs_phi",
                                     "Number of duplicated track candidates",
                                     m_cfg.varBinning.at("Phi"))),
      m_duplicationRatioVsPt("duplicationRatio_vs_pT",
                             "Duplication ratio;pT [GeV/c];Duplication ratio",
                             std::array{m_cfg.varBinning.at("Pt")}),
      m_duplicationRatioVsEta("duplicationRatio_vs_eta",
                              "Duplication ratio;#eta;Duplication ratio",
                              std::array{m_cfg.varBinning.at("Eta")}),
      m_duplicationRatioVsPhi("duplicationRatio_vs_phi",
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
                               std::size_t nMatchedTracks) {
  const auto t_phi = phi(truthParticle.direction());
  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();

  const auto nDuplicatedTracks = nMatchedTracks == 0 ? 0 : nMatchedTracks - 1;

  m_nDuplicatedVsPt.fill({t_pT}, static_cast<double>(nDuplicatedTracks));
  m_nDuplicatedVsEta.fill({t_eta}, static_cast<double>(nDuplicatedTracks));
  m_nDuplicatedVsPhi.fill({t_phi}, static_cast<double>(nDuplicatedTracks));
}

}  // namespace ActsExamples
