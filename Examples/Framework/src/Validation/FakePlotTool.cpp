// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/FakePlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace ActsExamples {

FakePlotTool::FakePlotTool(const FakePlotTool::Config& cfg,
                           Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("FakePlotTool", lvl)),
      m_nRecoVsPt("nRecoTracksVsPt", "Number of reconstructed track candidates",
                  std::array{m_cfg.varBinning.at("Pt"),
                             m_cfg.varBinning.at("Num")}),
      m_nTruthMatchedVsPt("nTruthMatchedTracksVsPt",
                          "Number of truth-matched track candidates",
                          std::array{m_cfg.varBinning.at("Pt"),
                                     m_cfg.varBinning.at("Num")}),
      m_nFakeVsPt("nFakeTracksVsPt", "Number of fake track candidates",
                  std::array{m_cfg.varBinning.at("Pt"),
                             m_cfg.varBinning.at("Num")}),
      m_nRecoVsEta("nRecoTracksVsEta",
                   "Number of reconstructed track candidates",
                   std::array{m_cfg.varBinning.at("Eta"),
                              m_cfg.varBinning.at("Num")}),
      m_nTruthMatchedVsEta("nTruthMatchedTracksVsEta",
                           "Number of truth-matched track candidates",
                           std::array{m_cfg.varBinning.at("Eta"),
                                      m_cfg.varBinning.at("Num")}),
      m_nFakeVsEta("nFakeTracksVsEta", "Number of fake track candidates",
                   std::array{m_cfg.varBinning.at("Eta"),
                              m_cfg.varBinning.at("Num")}),
      m_fakeRatioVsPt("fakeRatioVsPt",
                      "Tracking fake ratio;pT [GeV/c];Fake ratio",
                      std::array{m_cfg.varBinning.at("Pt")}),
      m_fakeRatioVsEta("fakeRatioVsEta", "Tracking fake ratio;#eta;Fake ratio",
                       std::array{m_cfg.varBinning.at("Eta")}),
      m_fakeRatioVsPhi("fakeRatioVsPhi", "Tracking fake ratio;#phi;Fake ratio",
                       std::array{m_cfg.varBinning.at("Phi")}) {
  ACTS_DEBUG("Initialize the histograms for fake ratio plots");
}

void FakePlotTool::fill(const Acts::BoundTrackParameters& fittedParameters,
                        bool status) {
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  m_fakeRatioVsPt.fill({fit_pT}, status);
  m_fakeRatioVsEta.fill({fit_eta}, status);
  m_fakeRatioVsPhi.fill({fit_phi}, status);
}

void FakePlotTool::fill(const SimParticleState& truthParticle,
                        std::size_t nTruthMatchedTracks,
                        std::size_t nFakeTracks) {
  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();

  m_nRecoVsPt.fill(
      {t_pT, static_cast<double>(nTruthMatchedTracks + nFakeTracks)});
  m_nTruthMatchedVsPt.fill({t_pT, static_cast<double>(nTruthMatchedTracks)});
  m_nFakeVsPt.fill({t_pT, static_cast<double>(nFakeTracks)});

  m_nRecoVsEta.fill(
      {t_eta, static_cast<double>(nTruthMatchedTracks + nFakeTracks)});
  m_nTruthMatchedVsEta.fill({t_eta, static_cast<double>(nTruthMatchedTracks)});
  m_nFakeVsEta.fill({t_eta, static_cast<double>(nFakeTracks)});
}

}  // namespace ActsExamples
