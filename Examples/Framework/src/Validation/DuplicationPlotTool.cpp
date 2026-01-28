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

struct DuplicationPlotTool::Impl {
  /// Number of duplicated tracks vs pT
  std::unique_ptr<TProfile> nDuplicated_vs_pT;
  /// Number of duplicated tracks vs eta
  std::unique_ptr<TProfile> nDuplicated_vs_eta;
  /// Number of duplicated tracks vs phi
  std::unique_ptr<TProfile> nDuplicated_vs_phi;
  /// Tracking duplication ratio vs pT
  std::unique_ptr<TEfficiency> duplicationRatio_vs_pT;
  /// Tracking duplication ratio vs eta
  std::unique_ptr<TEfficiency> duplicationRatio_vs_eta;
  /// Tracking duplication ratio vs phi
  std::unique_ptr<TEfficiency> duplicationRatio_vs_phi;
};

DuplicationPlotTool::DuplicationPlotTool(const DuplicationPlotTool::Config& cfg,
                                         Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("DuplicationPlotTool", lvl)),
      m_impl(std::make_unique<Impl>()) {}

DuplicationPlotTool::~DuplicationPlotTool() = default;

void DuplicationPlotTool::book() {
  Impl& cache = *m_impl;
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
  PlotHelpers::Binning bNum = m_cfg.varBinning.at("Num");
  ACTS_DEBUG("Initialize the histograms for duplication ratio plots");

  // duplication ratio vs pT
  cache.duplicationRatio_vs_pT = PlotHelpers::bookEff(
      "duplicationRatio_vs_pT",
      "Duplication ratio;pT [GeV/c];Duplication ratio", bPt);
  // duplication ratio vs eta
  cache.duplicationRatio_vs_eta =
      PlotHelpers::bookEff("duplicationRatio_vs_eta",
                           "Duplication ratio;#eta;Duplication ratio", bEta);
  // duplication ratio vs phi
  cache.duplicationRatio_vs_phi =
      PlotHelpers::bookEff("duplicationRatio_vs_phi",
                           "Duplication ratio;#phi;Duplication ratio", bPhi);

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

void DuplicationPlotTool::write() {
  Impl& cache = *m_impl;
  ACTS_DEBUG("Write the plots to output file.");
  cache.duplicationRatio_vs_pT->Write();
  cache.duplicationRatio_vs_eta->Write();
  cache.duplicationRatio_vs_phi->Write();
  cache.nDuplicated_vs_pT->Write();
  cache.nDuplicated_vs_eta->Write();
  cache.nDuplicated_vs_phi->Write();
}

void DuplicationPlotTool::fill(
    const Acts::BoundTrackParameters& fittedParameters, bool status) {
  Impl& cache = *m_impl;
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  PlotHelpers::fillEff(*cache.duplicationRatio_vs_pT, fit_pT, status);
  PlotHelpers::fillEff(*cache.duplicationRatio_vs_eta, fit_eta, status);
  PlotHelpers::fillEff(*cache.duplicationRatio_vs_phi, fit_phi, status);
}

void DuplicationPlotTool::fill(const SimParticleState& truthParticle,
                               std::size_t nMatchedTracks) {
  Impl& cache = *m_impl;
  const auto t_phi = phi(truthParticle.direction());
  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();

  const auto nDuplicatedTracks = nMatchedTracks == 0 ? 0 : nMatchedTracks - 1;

  PlotHelpers::fillProf(*cache.nDuplicated_vs_pT, t_pT, nDuplicatedTracks);
  PlotHelpers::fillProf(*cache.nDuplicated_vs_eta, t_eta, nDuplicatedTracks);
  PlotHelpers::fillProf(*cache.nDuplicated_vs_phi, t_phi, nDuplicatedTracks);
}

}  // namespace ActsExamples
