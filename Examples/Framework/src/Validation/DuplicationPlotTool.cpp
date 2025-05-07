// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/DuplicationPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <TEfficiency.h>
#include <TProfile.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

namespace ActsExamples {

DuplicationPlotTool::DuplicationPlotTool(const DuplicationPlotTool::Config& cfg,
                                         Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("DuplicationPlotTool", lvl)) {}

void DuplicationPlotTool::book(Cache& cache) const {
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
  PlotHelpers::Binning bNum = m_cfg.varBinning.at("Num");
  ACTS_DEBUG("Initialize the histograms for duplication rate plots");

  // duplication rate vs pT
  cache.duplicationRate_vs_pT =
      PlotHelpers::bookEff("duplicationRate_vs_pT",
                           "Duplication rate;pT [GeV/c];Duplication rate", bPt);
  // duplication rate vs eta
  cache.duplicationRate_vs_eta = PlotHelpers::bookEff(
      "duplicationRate_vs_eta", "Duplication rate;#eta;Duplication rate", bEta);
  // duplication rate vs phi
  cache.duplicationRate_vs_phi = PlotHelpers::bookEff(
      "duplicationRate_vs_phi", "Duplication rate;#phi;Duplication rate", bPhi);

  // duplication number vs pT
  cache.nDuplicated_vs_pT = PlotHelpers::bookProf(
      "nDuplicated_vs_pT", "Number of duplicated track candidates", bPt, bNum);
  // duplication number vs eta
  cache.nDuplicated_vs_eta = PlotHelpers::bookProf(
      "nDuplicated_vs_eta", "Number of duplicated track candidates", bEta,
      bNum);
  // duplication number vs phi
  cache.nDuplicated_vs_phi = PlotHelpers::bookProf(
      "nDuplicated_vs_phi", "Number of duplicated track candidates", bPhi,
      bNum);
}

void DuplicationPlotTool::clear(Cache& cache) const {
  delete cache.duplicationRate_vs_pT;
  delete cache.duplicationRate_vs_eta;
  delete cache.duplicationRate_vs_phi;
  delete cache.nDuplicated_vs_pT;
  delete cache.nDuplicated_vs_eta;
  delete cache.nDuplicated_vs_phi;
}

void DuplicationPlotTool::write(const Cache& cache) const {
  ACTS_DEBUG("Write the plots to output file.");
  cache.duplicationRate_vs_pT->Write();
  cache.duplicationRate_vs_eta->Write();
  cache.duplicationRate_vs_phi->Write();
  cache.nDuplicated_vs_pT->Write();
  cache.nDuplicated_vs_eta->Write();
  cache.nDuplicated_vs_phi->Write();
}

void DuplicationPlotTool::fill(
    Cache& cache, const Acts::BoundTrackParameters& fittedParameters,
    bool status) const {
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  PlotHelpers::fillEff(cache.duplicationRate_vs_pT, fit_pT, status);
  PlotHelpers::fillEff(cache.duplicationRate_vs_eta, fit_eta, status);
  PlotHelpers::fillEff(cache.duplicationRate_vs_phi, fit_phi, status);
}

void DuplicationPlotTool::fill(Cache& cache,
                               const SimParticleState& truthParticle,
                               std::size_t nDuplicatedTracks) const {
  const auto t_phi = phi(truthParticle.direction());
  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();

  PlotHelpers::fillProf(cache.nDuplicated_vs_pT, t_pT, nDuplicatedTracks);
  PlotHelpers::fillProf(cache.nDuplicated_vs_eta, t_eta, nDuplicatedTracks);
  PlotHelpers::fillProf(cache.nDuplicated_vs_phi, t_phi, nDuplicatedTracks);
}

}  // namespace ActsExamples
