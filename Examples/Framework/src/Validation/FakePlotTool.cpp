// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/FakePlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"

#include <TEfficiency.h>
#include <TH2.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

namespace ActsExamples {

FakePlotTool::FakePlotTool(const FakePlotTool::Config& cfg,
                           Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("FakePlotTool", lvl)) {}

void FakePlotTool::book(Cache& cache) const {
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
  PlotHelpers::Binning bNum = m_cfg.varBinning.at("Num");
  ACTS_DEBUG("Initialize the histograms for fake ratio plots");

  // number of reco tracks vs pT scatter plots
  cache.nReco_vs_pT = PlotHelpers::bookHisto(
      "nRecoTracks_vs_pT", "Number of reconstructed track candidates", bPt,
      bNum);
  // number of truth-matched tracks vs pT scatter plots
  cache.nTruthMatched_vs_pT = PlotHelpers::bookHisto(
      "nTruthMatchedTracks_vs_pT", "Number of truth-matched track candidates",
      bPt, bNum);
  // number of fake tracks vs pT scatter plots
  cache.nFake_vs_pT = PlotHelpers::bookHisto(
      "nFakeTracks_vs_pT", "Number of fake track candidates", bPt, bNum);

  // number of reco tracks vs eta scatter plots
  cache.nReco_vs_eta = PlotHelpers::bookHisto(
      "nRecoTracks_vs_eta", "Number of reconstructed track candidates", bEta,
      bNum);
  // number of truth-matched tracks vs eta scatter plots
  cache.nTruthMatched_vs_eta = PlotHelpers::bookHisto(
      "nTruthMatchedTracks_vs_eta", "Number of truth-matched track candidates",
      bEta, bNum);
  // number of fake tracks vs eta scatter plots
  cache.nFake_vs_eta = PlotHelpers::bookHisto(
      "nFakeTracks_vs_eta", "Number of fake track candidates", bEta, bNum);

  // fake ratio vs pT
  cache.fakeRatio_vs_pT = PlotHelpers::bookEff(
      "fakeRatio_vs_pT", "Tracking fake ratio;pT [GeV/c];Fake ratio", bPt);
  // fake ratio vs eta
  cache.fakeRatio_vs_eta = PlotHelpers::bookEff(
      "fakeRatio_vs_eta", "Tracking fake ratio;#eta;Fake ratio", bEta);
  // fake ratio vs phi
  cache.fakeRatio_vs_phi = PlotHelpers::bookEff(
      "fakeRatio_vs_phi", "Tracking fake ratio;#phi;Fake ratio", bPhi);
}

void FakePlotTool::clear(Cache& cache) const {
  delete cache.nReco_vs_pT;
  delete cache.nTruthMatched_vs_pT;
  delete cache.nFake_vs_pT;
  delete cache.nReco_vs_eta;
  delete cache.nTruthMatched_vs_eta;
  delete cache.nFake_vs_eta;
  delete cache.fakeRatio_vs_pT;
  delete cache.fakeRatio_vs_eta;
  delete cache.fakeRatio_vs_phi;
}

void FakePlotTool::write(const Cache& cache) const {
  ACTS_DEBUG("Write the plots to output file.");
  cache.nReco_vs_pT->Write();
  cache.nTruthMatched_vs_pT->Write();
  cache.nFake_vs_pT->Write();
  cache.nReco_vs_eta->Write();
  cache.nTruthMatched_vs_eta->Write();
  cache.nFake_vs_eta->Write();
  cache.fakeRatio_vs_pT->Write();
  cache.fakeRatio_vs_eta->Write();
  cache.fakeRatio_vs_phi->Write();
}

void FakePlotTool::fill(Cache& cache,
                        const Acts::BoundTrackParameters& fittedParameters,
                        bool status) const {
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  PlotHelpers::fillEff(cache.fakeRatio_vs_pT, fit_pT, status);
  PlotHelpers::fillEff(cache.fakeRatio_vs_eta, fit_eta, status);
  PlotHelpers::fillEff(cache.fakeRatio_vs_phi, fit_phi, status);
}

void FakePlotTool::fill(Cache& cache, const SimParticleState& truthParticle,
                        std::size_t nTruthMatchedTracks,
                        std::size_t nFakeTracks) const {
  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();

  PlotHelpers::fillHisto(cache.nReco_vs_pT, t_pT,
                         nTruthMatchedTracks + nFakeTracks);
  PlotHelpers::fillHisto(cache.nTruthMatched_vs_pT, t_pT, nTruthMatchedTracks);
  PlotHelpers::fillHisto(cache.nFake_vs_pT, t_pT, nFakeTracks);

  PlotHelpers::fillHisto(cache.nReco_vs_eta, t_eta,
                         nTruthMatchedTracks + nFakeTracks);
  PlotHelpers::fillHisto(cache.nTruthMatched_vs_eta, t_eta,
                         nTruthMatchedTracks);
  PlotHelpers::fillHisto(cache.nFake_vs_eta, t_eta, nFakeTracks);
}

}  // namespace ActsExamples
