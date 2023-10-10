// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/DuplicationPlotTool.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <TEfficiency.h>
#include <TProfile.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

ActsExamples::DuplicationPlotTool::DuplicationPlotTool(
    const ActsExamples::DuplicationPlotTool::Config& cfg,
    Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("DuplicationPlotTool", lvl)) {}

void ActsExamples::DuplicationPlotTool::book(
    DuplicationPlotTool::DuplicationPlotCache& duplicationPlotCache) const {
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
  PlotHelpers::Binning bNum = m_cfg.varBinning.at("Num");
  ACTS_DEBUG("Initialize the histograms for duplication rate plots");

  // duplication rate vs pT
  duplicationPlotCache.duplicationRate_vs_pT =
      PlotHelpers::bookEff("duplicationRate_vs_pT",
                           "Duplication rate;pT [GeV/c];Duplication rate", bPt);
  // duplication rate vs eta
  duplicationPlotCache.duplicationRate_vs_eta = PlotHelpers::bookEff(
      "duplicationRate_vs_eta", "Duplication rate;#eta;Duplication rate", bEta);
  // duplication rate vs phi
  duplicationPlotCache.duplicationRate_vs_phi = PlotHelpers::bookEff(
      "duplicationRate_vs_phi", "Duplication rate;#phi;Duplication rate", bPhi);

  // duplication number vs pT
  duplicationPlotCache.nDuplicated_vs_pT = PlotHelpers::bookProf(
      "nDuplicated_vs_pT", "Number of duplicated track candidates", bPt, bNum);
  // duplication number vs eta
  duplicationPlotCache.nDuplicated_vs_eta = PlotHelpers::bookProf(
      "nDuplicated_vs_eta", "Number of duplicated track candidates", bEta,
      bNum);
  // duplication number vs phi
  duplicationPlotCache.nDuplicated_vs_phi = PlotHelpers::bookProf(
      "nDuplicated_vs_phi", "Number of duplicated track candidates", bPhi,
      bNum);
}

void ActsExamples::DuplicationPlotTool::clear(
    DuplicationPlotCache& duplicationPlotCache) const {
  delete duplicationPlotCache.duplicationRate_vs_pT;
  delete duplicationPlotCache.duplicationRate_vs_eta;
  delete duplicationPlotCache.duplicationRate_vs_phi;
  delete duplicationPlotCache.nDuplicated_vs_pT;
  delete duplicationPlotCache.nDuplicated_vs_eta;
  delete duplicationPlotCache.nDuplicated_vs_phi;
}

void ActsExamples::DuplicationPlotTool::write(
    const DuplicationPlotTool::DuplicationPlotCache& duplicationPlotCache)
    const {
  ACTS_DEBUG("Write the plots to output file.");
  duplicationPlotCache.duplicationRate_vs_pT->Write();
  duplicationPlotCache.duplicationRate_vs_eta->Write();
  duplicationPlotCache.duplicationRate_vs_phi->Write();
  duplicationPlotCache.nDuplicated_vs_pT->Write();
  duplicationPlotCache.nDuplicated_vs_eta->Write();
  duplicationPlotCache.nDuplicated_vs_phi->Write();
}

void ActsExamples::DuplicationPlotTool::fill(
    DuplicationPlotTool::DuplicationPlotCache& duplicationPlotCache,
    const Acts::BoundTrackParameters& fittedParameters, bool status) const {
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  PlotHelpers::fillEff(duplicationPlotCache.duplicationRate_vs_pT, fit_pT,
                       status);
  PlotHelpers::fillEff(duplicationPlotCache.duplicationRate_vs_eta, fit_eta,
                       status);
  PlotHelpers::fillEff(duplicationPlotCache.duplicationRate_vs_phi, fit_phi,
                       status);
}

void ActsExamples::DuplicationPlotTool::fill(
    DuplicationPlotTool::DuplicationPlotCache& duplicationPlotCache,
    const ActsFatras::Particle& truthParticle, size_t nDuplicatedTracks) const {
  const auto t_phi = phi(truthParticle.direction());
  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();

  PlotHelpers::fillProf(duplicationPlotCache.nDuplicated_vs_pT, t_pT,
                        nDuplicatedTracks);
  PlotHelpers::fillProf(duplicationPlotCache.nDuplicated_vs_eta, t_eta,
                        nDuplicatedTracks);
  PlotHelpers::fillProf(duplicationPlotCache.nDuplicated_vs_phi, t_phi,
                        nDuplicatedTracks);
}
