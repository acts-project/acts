// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/EffPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <TEfficiency.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace ActsExamples {

EffPlotTool::EffPlotTool(const EffPlotTool::Config& cfg,
                         Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("EffPlotTool", lvl)) {}

void EffPlotTool::book(Cache& cache) const {
  PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  PlotHelpers::Binning bDeltaR = m_cfg.varBinning.at("DeltaR");
  PlotHelpers::Binning bZ0 = m_cfg.varBinning.at("Z0");
  PlotHelpers::Binning bProdR = m_cfg.varBinning.at("prodR");
  ACTS_DEBUG("Initialize the histograms for efficiency plots");
  // efficiency vs pT
  cache.trackEff_vs_pT = PlotHelpers::bookEff(
      "trackeff_vs_pT", "Tracking efficiency;Truth pT [GeV/c];Efficiency", bPt);
  // efficiency vs eta
  cache.trackEff_vs_eta = PlotHelpers::bookEff(
      "trackeff_vs_eta", "Tracking efficiency;Truth #eta;Efficiency", bEta);
  // efficiency vs phi
  cache.trackEff_vs_phi = PlotHelpers::bookEff(
      "trackeff_vs_phi", "Tracking efficiency;Truth #phi;Efficiency", bPhi);
  // efficiency vs z0
  cache.trackEff_vs_z0 = PlotHelpers::bookEff(
      "trackeff_vs_z0", "Tracking efficiency;Truth z_0 [mm];Efficiency", bZ0);
  // efficiancy vs distance to the closest truth particle
  cache.trackEff_vs_DeltaR = PlotHelpers::bookEff(
      "trackeff_vs_DeltaR",
      "Tracking efficiency;Closest track #Delta R;Efficiency", bDeltaR);
  cache.trackEff_vs_prodR = PlotHelpers::bookEff(
      "trackeff_vs_prodR",
      "Tracking efficiency;Production radius [mm];Efficiency", bProdR);
}

void EffPlotTool::clear(Cache& cache) const {
  delete cache.trackEff_vs_pT;
  delete cache.trackEff_vs_eta;
  delete cache.trackEff_vs_phi;
  delete cache.trackEff_vs_z0;
  delete cache.trackEff_vs_DeltaR;
  delete cache.trackEff_vs_prodR;
}

void EffPlotTool::write(const Cache& cache) const {
  ACTS_DEBUG("Write the plots to output file.");
  cache.trackEff_vs_pT->Write();
  cache.trackEff_vs_eta->Write();
  cache.trackEff_vs_phi->Write();
  cache.trackEff_vs_z0->Write();
  cache.trackEff_vs_DeltaR->Write();
  cache.trackEff_vs_prodR->Write();
}

void EffPlotTool::fill(Cache& cache, const SimParticleState& truthParticle,
                       double deltaR, bool status) const {
  const auto t_phi = phi(truthParticle.direction());
  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();
  const auto t_z0 = truthParticle.position().z();
  const auto t_deltaR = deltaR;
  const auto t_prodR =
      std::sqrt(truthParticle.position().x() * truthParticle.position().x() +
                truthParticle.position().y() * truthParticle.position().y());

  PlotHelpers::fillEff(cache.trackEff_vs_pT, t_pT, status);
  PlotHelpers::fillEff(cache.trackEff_vs_eta, t_eta, status);
  PlotHelpers::fillEff(cache.trackEff_vs_phi, t_phi, status);
  PlotHelpers::fillEff(cache.trackEff_vs_z0, t_z0, status);
  PlotHelpers::fillEff(cache.trackEff_vs_DeltaR, t_deltaR, status);
  PlotHelpers::fillEff(cache.trackEff_vs_prodR, t_prodR, status);
}

}  // namespace ActsExamples
