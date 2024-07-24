// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/TrkVtxPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <TEfficiency.h>


ActsExamples::TrkVtxPlotTool::TrkVtxPlotTool(
    const ActsExamples::TrkVtxPlotTool::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("TrkVtxPlotTool", lvl)) {}

void ActsExamples::TrkVtxPlotTool::book(
    TrkVtxPlotTool::TrkVtxPlotCache& trkVtxPlotCache) const {
  PlotHelpers::Binning bZgen = m_cfg.varBinning.at("Zgen");
  PlotHelpers::Binning bZrec = m_cfg.varBinning.at("Zrec");
  PlotHelpers::Binning bRes = m_cfg.varBinning.at("Res");
  PlotHelpers::Binning bTar = m_cfg.varBinning.at("Tar");
  ACTS_DEBUG("Initialize the histograms for efficiency plots");
  // efficiency vs pT
  trkVtxPlotCache.resVtxz = PlotHelpers::bookHisto(
      "resVtxz", "Residuals; z_{rec}-z_{gen};Counts", bRes);
  trkVtxPlotCache.resVtxz_vs_zgen = PlotHelpers::bookHisto(
      "resVtxz_vs_zgen", "Residuals; z_{gen};z_{rec}-z_{gen};Counts", bZgen, bRes);
  trkVtxPlotCache.zrec_vs_zgen = PlotHelpers::bookHisto(
      "zrec_vs_zgen", "; z_{gen};z_{rec};Counts", bZgen, bZrec);
  trkVtxPlotCache.eff_vs_zgen = PlotHelpers::bookEff(
      "eff_vs_zgen", "Tracklet vertexing efficiency; z_{gen};Efficiency", bTar);
}

void ActsExamples::TrkVtxPlotTool::clear(TrkVtxPlotCache& trkVtxPlotCache) const {
  delete trkVtxPlotCache.resVtxz;
  delete trkVtxPlotCache.resVtxz_vs_zgen;
  delete trkVtxPlotCache.zrec_vs_zgen;
  delete trkVtxPlotCache.eff_vs_zgen;
}

void ActsExamples::TrkVtxPlotTool::write(
    const TrkVtxPlotTool::TrkVtxPlotCache& trkVtxPlotCache) const {
  ACTS_DEBUG("Write the plots to output file.");
  trkVtxPlotCache.resVtxz->Write();
  trkVtxPlotCache.resVtxz_vs_zgen->Write();
  trkVtxPlotCache.zrec_vs_zgen->Write();
  trkVtxPlotCache.eff_vs_zgen->Write();
}

void ActsExamples::TrkVtxPlotTool::fill(TrkVtxPlotTool::TrkVtxPlotCache& trkVtxPlotCache, double zrec, double zgen) const {
  PlotHelpers::fillHisto(trkVtxPlotCache.resVtxz, zrec-zgen, 1);
  PlotHelpers::fillHisto(trkVtxPlotCache.resVtxz_vs_zgen, zgen, zrec-zgen, 1);
  PlotHelpers::fillHisto(trkVtxPlotCache.zrec_vs_zgen, zgen, zrec, 1);
  PlotHelpers::fillEff(trkVtxPlotCache.eff_vs_zgen, zgen, abs(zrec-zgen)<6);
}
