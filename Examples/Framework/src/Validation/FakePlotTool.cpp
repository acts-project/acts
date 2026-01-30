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
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("FakePlotTool", lvl)) {
  ACTS_DEBUG("Initialize the histograms for fake ratio plots");

  using Histogram2 = Acts::Experimental::Histogram2;
  using Efficiency1 = Acts::Experimental::Efficiency1;

  m_histograms.insert(
      {"nRecoTracks_vs_pT",
       Histogram2(
           "nRecoTracks_vs_pT", "Number of reconstructed track candidates",
           std::array{m_cfg.varBinning.at("Pt"), m_cfg.varBinning.at("Num")})});
  m_histograms.insert(
      {"nTruthMatchedTracks_vs_pT",
       Histogram2(
           "nTruthMatchedTracks_vs_pT",
           "Number of truth-matched track candidates",
           std::array{m_cfg.varBinning.at("Pt"), m_cfg.varBinning.at("Num")})});
  m_histograms.insert(
      {"nFakeTracks_vs_pT",
       Histogram2(
           "nFakeTracks_vs_pT", "Number of fake track candidates",
           std::array{m_cfg.varBinning.at("Pt"), m_cfg.varBinning.at("Num")})});
  m_histograms.insert(
      {"nRecoTracks_vs_eta",
       Histogram2("nRecoTracks_vs_eta",
                  "Number of reconstructed track candidates",
                  std::array{m_cfg.varBinning.at("Eta"),
                             m_cfg.varBinning.at("Num")})});
  m_histograms.insert(
      {"nTruthMatchedTracks_vs_eta",
       Histogram2(
           "nTruthMatchedTracks_vs_eta",
           "Number of truth-matched track candidates",
           std::array{m_cfg.varBinning.at("Eta"), m_cfg.varBinning.at("Num")})});
  m_histograms.insert(
      {"nFakeTracks_vs_eta",
       Histogram2("nFakeTracks_vs_eta", "Number of fake track candidates",
                  std::array{m_cfg.varBinning.at("Eta"),
                             m_cfg.varBinning.at("Num")})});

  m_efficiencies.insert(
      {"fakeRatio_vs_pT",
       Efficiency1("fakeRatio_vs_pT", "Tracking fake ratio",
                   std::array{m_cfg.varBinning.at("Pt")})});
  m_efficiencies.insert(
      {"fakeRatio_vs_eta",
       Efficiency1("fakeRatio_vs_eta", "Tracking fake ratio",
                   std::array{m_cfg.varBinning.at("Eta")})});
  m_efficiencies.insert(
      {"fakeRatio_vs_phi",
       Efficiency1("fakeRatio_vs_phi", "Tracking fake ratio",
                   std::array{m_cfg.varBinning.at("Phi")})});
}

void FakePlotTool::fill(const Acts::BoundTrackParameters& fittedParameters,
                        bool status) {
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  m_efficiencies.at("fakeRatio_vs_pT").fill({fit_pT}, status);
  m_efficiencies.at("fakeRatio_vs_eta").fill({fit_eta}, status);
  m_efficiencies.at("fakeRatio_vs_phi").fill({fit_phi}, status);
}

void FakePlotTool::fill(const SimParticleState& truthParticle,
                        std::size_t nTruthMatchedTracks,
                        std::size_t nFakeTracks) {
  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();

  m_histograms.at("nRecoTracks_vs_pT")
      .fill({t_pT, static_cast<double>(nTruthMatchedTracks + nFakeTracks)});
  m_histograms.at("nTruthMatchedTracks_vs_pT")
      .fill({t_pT, static_cast<double>(nTruthMatchedTracks)});
  m_histograms.at("nFakeTracks_vs_pT")
      .fill({t_pT, static_cast<double>(nFakeTracks)});

  m_histograms.at("nRecoTracks_vs_eta")
      .fill({t_eta, static_cast<double>(nTruthMatchedTracks + nFakeTracks)});
  m_histograms.at("nTruthMatchedTracks_vs_eta")
      .fill({t_eta, static_cast<double>(nTruthMatchedTracks)});
  m_histograms.at("nFakeTracks_vs_eta")
      .fill({t_eta, static_cast<double>(nFakeTracks)});
}

}  // namespace ActsExamples
