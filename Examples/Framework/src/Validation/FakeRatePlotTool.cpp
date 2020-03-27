// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Validation/FakeRatePlotTool.hpp"

#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

FW::FakeRatePlotTool::FakeRatePlotTool(const FW::FakeRatePlotTool::Config& cfg,
                                       Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("FakeRatePlotTool", lvl)) {}

void FW::FakeRatePlotTool::book(
    FakeRatePlotTool::FakeRatePlotCache& fakeRatePlotCache) const {
  PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  PlotHelpers::Binning bNum = m_cfg.varBinning.at("Num");
  ACTS_DEBUG("Initialize the histograms for fake rate plots");

  // number fo reco tracks
  fakeRatePlotCache.nRecoTracks = PlotHelpers::bookHisto(
      "nRecoTracks", "Number of reconstructed track candidates", bNum);
  // number fo truth-matched tracks
  fakeRatePlotCache.nTruthMatchedTracks = PlotHelpers::bookHisto(
      "nTruthMatchedTracks", "Number of truth-matched track candidates", bNum);
  // number fo fake tracks
  fakeRatePlotCache.nFakeTracks = PlotHelpers::bookHisto(
      "nFakeTracks", "Number of fake track candidates", bNum);

  // fake rate vs pT
  fakeRatePlotCache.fakerate_vs_pT = PlotHelpers::bookEff(
      "fakerate_vs_pT", "Tracking fake rate;pT [GeV/c];Fake rate", bPt);
  // fake rate vs eta
  fakeRatePlotCache.fakerate_vs_eta = PlotHelpers::bookEff(
      "fakerate_vs_eta", "Tracking fake rate;#eta;Fake rate", bEta);
  // fake rate vs phi
  fakeRatePlotCache.fakerate_vs_phi = PlotHelpers::bookEff(
      "fakerate_vs_phi", "Tracking fake rate;#phi;Fake rate", bPhi);

  // duplication number vs pT
  fakeRatePlotCache.duplicationNum_vs_pT = PlotHelpers::bookProf(
      "duplicationNum_vs_pT", "Duplication number vs. pT", bPt, bNum);
  // duplication number vs eta
  fakeRatePlotCache.duplicationNum_vs_eta = PlotHelpers::bookProf(
      "duplicationNum_vs_eta", "Duplication number vs. #eta", bEta, bNum);
  // duplication number vs phi
  fakeRatePlotCache.duplicationNum_vs_phi = PlotHelpers::bookProf(
      "duplicationNum_vs_phi", "Duplication number vs. #phi", bPhi, bNum);
}

void FW::FakeRatePlotTool::clear(FakeRatePlotCache& fakeRatePlotCache) const {
  delete fakeRatePlotCache.nRecoTracks;
  delete fakeRatePlotCache.nTruthMatchedTracks;
  delete fakeRatePlotCache.nFakeTracks;
  delete fakeRatePlotCache.fakerate_vs_pT;
  delete fakeRatePlotCache.fakerate_vs_eta;
  delete fakeRatePlotCache.fakerate_vs_phi;
  delete fakeRatePlotCache.duplicationNum_vs_pT;
  delete fakeRatePlotCache.duplicationNum_vs_eta;
  delete fakeRatePlotCache.duplicationNum_vs_phi;
}

void FW::FakeRatePlotTool::write(
    const FakeRatePlotTool::FakeRatePlotCache& fakeRatePlotCache) const {
  ACTS_DEBUG("Write the plots to output file.");
  fakeRatePlotCache.fakerate_vs_pT->Write();
  fakeRatePlotCache.fakerate_vs_eta->Write();
  fakeRatePlotCache.fakerate_vs_phi->Write();
  fakeRatePlotCache.nRecoTracks->Write();
  fakeRatePlotCache.nTruthMatchedTracks->Write();
  fakeRatePlotCache.nFakeTracks->Write();
  fakeRatePlotCache.duplicationNum_vs_pT->Write();
  fakeRatePlotCache.duplicationNum_vs_eta->Write();
  fakeRatePlotCache.duplicationNum_vs_phi->Write();
}

void FW::FakeRatePlotTool::fill(
    FakeRatePlotTool::FakeRatePlotCache& fakeRatePlotCache,
    const ActsFatras::Particle& truthParticle, bool status) const {
  const auto t_phi = phi(truthParticle.unitDirection());
  const auto t_eta = eta(truthParticle.unitDirection());
  const auto t_pT = truthParticle.transverseMomentum();

  PlotHelpers::fillEff(fakeRatePlotCache.fakerate_vs_pT, t_pT, status);
  PlotHelpers::fillEff(fakeRatePlotCache.fakerate_vs_eta, t_eta, status);
  PlotHelpers::fillEff(fakeRatePlotCache.fakerate_vs_phi, t_phi, status);
}

void FW::FakeRatePlotTool::fill(
    FakeRatePlotTool::FakeRatePlotCache& fakeRatePlotCache,
    const ActsFatras::Particle& truthParticle, size_t nTruthMatchedTracks,
    size_t nFakeTracks) const {
  const auto t_phi = phi(truthParticle.unitDirection());
  const auto t_eta = eta(truthParticle.unitDirection());
  const auto t_pT = truthParticle.transverseMomentum();

  PlotHelpers::fillHisto(fakeRatePlotCache.nRecoTracks,
                         nTruthMatchedTracks + nFakeTracks);
  PlotHelpers::fillHisto(fakeRatePlotCache.nTruthMatchedTracks,
                         nTruthMatchedTracks);
  PlotHelpers::fillHisto(fakeRatePlotCache.nFakeTracks, nFakeTracks);

  PlotHelpers::fillProf(fakeRatePlotCache.duplicationNum_vs_pT, t_pT,
                        nTruthMatchedTracks);
  PlotHelpers::fillProf(fakeRatePlotCache.duplicationNum_vs_eta, t_eta,
                        nTruthMatchedTracks);
  PlotHelpers::fillProf(fakeRatePlotCache.duplicationNum_vs_phi, t_phi,
                        nTruthMatchedTracks);
}
