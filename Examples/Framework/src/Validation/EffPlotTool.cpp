// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Validation/EffPlotTool.hpp"

#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

FW::EffPlotTool::EffPlotTool(const FW::EffPlotTool::Config& cfg,
                             Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("EffPlotTool", lvl)) {}

void FW::EffPlotTool::book(EffPlotTool::EffPlotCache& effPlotCache) const {
  PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  ACTS_DEBUG("Initialize the histograms for efficiency plots");
  // efficiency vs pT
  effPlotCache.trackeff_vs_pT = PlotHelpers::bookEff(
      "trackeff_vs_pT", "Tracking efficiency;pT [GeV/c];Efficiency", bPt);
  // efficiency vs eta
  effPlotCache.trackeff_vs_eta = PlotHelpers::bookEff(
      "trackeff_vs_eta", "Tracking efficiency;#eta;Efficiency", bEta);
  // efficiency vs phi
  effPlotCache.trackeff_vs_phi = PlotHelpers::bookEff(
      "trackeff_vs_phi", "Tracking efficiency;#phi;Efficiency", bPhi);
}

void FW::EffPlotTool::clear(EffPlotCache& effPlotCache) const {
  delete effPlotCache.trackeff_vs_pT;
  delete effPlotCache.trackeff_vs_eta;
  delete effPlotCache.trackeff_vs_phi;
}

void FW::EffPlotTool::write(
    const EffPlotTool::EffPlotCache& effPlotCache) const {
  ACTS_DEBUG("Write the plots to output file.");
  effPlotCache.trackeff_vs_pT->Write();
  effPlotCache.trackeff_vs_eta->Write();
  effPlotCache.trackeff_vs_phi->Write();
}

void FW::EffPlotTool::fill(EffPlotTool::EffPlotCache& effPlotCache,
                           const ActsFatras::Particle& truthParticle,
                           bool status) const {
  const auto t_phi = phi(truthParticle.unitDirection());
  const auto t_eta = eta(truthParticle.unitDirection());
  const auto t_pT = truthParticle.transverseMomentum();

  PlotHelpers::fillEff(effPlotCache.trackeff_vs_pT, t_pT, status);
  PlotHelpers::fillEff(effPlotCache.trackeff_vs_eta, t_eta, status);
  PlotHelpers::fillEff(effPlotCache.trackeff_vs_phi, t_phi, status);
}
