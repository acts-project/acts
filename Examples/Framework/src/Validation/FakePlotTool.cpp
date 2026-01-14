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
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("FakePlotTool", lvl)) {}

void FakePlotTool::book(Cache& cache) const {
  const auto& ptAxis = m_cfg.varBinning.at("Pt");
  const auto& etaAxis = m_cfg.varBinning.at("Eta");
  const auto& phiAxis = m_cfg.varBinning.at("Phi");
  const auto& numAxis = m_cfg.varBinning.at("Num");
  ACTS_DEBUG("Initialize the histograms for fake ratio plots");

  // number of reco tracks vs pT scatter plots
  cache.nReco_vs_pT.emplace("nRecoTracks_vs_pT",
                            "Number of reconstructed track candidates",
                            std::array{ptAxis, numAxis});

  // number of truth-matched tracks vs pT scatter plots
  cache.nTruthMatched_vs_pT.emplace(
      "nTruthMatchedTracks_vs_pT", "Number of truth-matched track candidates",
      std::array{ptAxis, numAxis});

  // number of fake tracks vs pT scatter plots
  cache.nFake_vs_pT.emplace("nFakeTracks_vs_pT",
                            "Number of fake track candidates",
                            std::array{ptAxis, numAxis});

  // number of reco tracks vs eta scatter plots
  cache.nReco_vs_eta.emplace("nRecoTracks_vs_eta",
                             "Number of reconstructed track candidates",
                             std::array{etaAxis, numAxis});

  // number of truth-matched tracks vs eta scatter plots
  cache.nTruthMatched_vs_eta.emplace(
      "nTruthMatchedTracks_vs_eta", "Number of truth-matched track candidates",
      std::array{etaAxis, numAxis});

  // number of fake tracks vs eta scatter plots
  cache.nFake_vs_eta.emplace("nFakeTracks_vs_eta",
                             "Number of fake track candidates",
                             std::array{etaAxis, numAxis});

  // fake ratio vs pT
  cache.fakeRatio_vs_pT.emplace("fakeRatio_vs_pT",
                                "Tracking fake ratio;pT [GeV/c];Fake ratio",
                                std::array{ptAxis});

  // fake ratio vs eta
  cache.fakeRatio_vs_eta.emplace("fakeRatio_vs_eta",
                                 "Tracking fake ratio;#eta;Fake ratio",
                                 std::array{etaAxis});

  // fake ratio vs phi
  cache.fakeRatio_vs_phi.emplace("fakeRatio_vs_phi",
                                 "Tracking fake ratio;#phi;Fake ratio",
                                 std::array{phiAxis});
}

void FakePlotTool::fill(Cache& cache,
                        const Acts::BoundTrackParameters& fittedParameters,
                        bool status) const {
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  cache.fakeRatio_vs_pT->fill({fit_pT}, status);
  cache.fakeRatio_vs_eta->fill({fit_eta}, status);
  cache.fakeRatio_vs_phi->fill({fit_phi}, status);
}

void FakePlotTool::fill(Cache& cache, const SimParticleState& truthParticle,
                        std::size_t nTruthMatchedTracks,
                        std::size_t nFakeTracks) const {
  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();

  cache.nReco_vs_pT->fill(
      {t_pT, static_cast<double>(nTruthMatchedTracks + nFakeTracks)});
  cache.nTruthMatched_vs_pT->fill(
      {t_pT, static_cast<double>(nTruthMatchedTracks)});
  cache.nFake_vs_pT->fill({t_pT, static_cast<double>(nFakeTracks)});

  cache.nReco_vs_eta->fill(
      {t_eta, static_cast<double>(nTruthMatchedTracks + nFakeTracks)});
  cache.nTruthMatched_vs_eta->fill(
      {t_eta, static_cast<double>(nTruthMatchedTracks)});
  cache.nFake_vs_eta->fill({t_eta, static_cast<double>(nFakeTracks)});
}

}  // namespace ActsExamples
