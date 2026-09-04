// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/FakePlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"

#include <format>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace {

std::string capitalize(std::string_view s) {
  if (s.empty()) {
    return std::string{};
  }
  std::string r{s};
  r[0] = std::toupper(r[0]);
  return r;
}

}  // namespace

namespace ActsExamples {

FakePlotTool::FakePlotTool(const FakePlotTool::Config& cfg,
                           Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("FakePlotTool", lvl)) {
  ACTS_DEBUG("Initialize the histograms for fake ratio plots");

  std::string lt = capitalize(m_cfg.label);
  std::string lp = std::format("{}s", m_cfg.label);
  std::string rt = std::format("nReco{}s", lt);
  std::string tm = std::format("nTruthMatched{}s", lt);
  std::string ft = std::format("nFake{}s", lt);

  auto ins = [&](const std::string& stem, const std::string& nameSuffix,
                 const std::string& axisKey, const std::string& title) {
    std::string name = std::format("{}_{}", stem, nameSuffix);
    m_histograms.emplace(name,
                         Histogram2(name, title,
                                    std::array{m_cfg.varBinning.at(axisKey),
                                               m_cfg.varBinning.at("Num")}));
  };

  std::string recT = std::format("Number of reconstructed {} candidates", lp);
  std::string truT = std::format("Number of truth-matched {} candidates", lp);
  std::string fakT = std::format("Number of fake {} candidates", lp);

  ins(rt, "vs_pT", "Pt", recT);
  ins(tm, "vs_pT", "Pt", truT);
  ins(ft, "vs_pT", "Pt", fakT);
  ins(rt, "vs_eta", "Eta", recT);
  ins(tm, "vs_eta", "Eta", truT);
  ins(ft, "vs_eta", "Eta", fakT);

  std::string fakeRatioTitle = std::format("{} fake ratio", m_cfg.label);

  m_efficiencies.insert(
      {"fakeRatio_vs_pT",
       Efficiency1("fakeRatio_vs_pT", fakeRatioTitle,
                   std::array{m_cfg.recoVarBinning.at("Pt")})});
  m_efficiencies.insert(
      {"fakeRatio_vs_eta",
       Efficiency1("fakeRatio_vs_eta", fakeRatioTitle,
                   std::array{m_cfg.recoVarBinning.at("Eta")})});
  m_efficiencies.insert(
      {"fakeRatio_vs_phi",
       Efficiency1("fakeRatio_vs_phi", fakeRatioTitle,
                   std::array{m_cfg.recoVarBinning.at("Phi")})});
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
  std::string labelTitle = capitalize(m_cfg.label);

  auto key = [&](const std::string& stem,
                 const std::string& var) -> std::string {
    return std::format("{}_{}", stem, var);
  };

  std::string recoPrefix = std::format("nReco{}s", labelTitle);
  std::string truthPrefix = std::format("nTruthMatched{}s", labelTitle);
  std::string fakePrefix = std::format("nFake{}s", labelTitle);

  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();

  m_histograms.at(key(recoPrefix, "vs_pT"))
      .fill({t_pT, static_cast<double>(nTruthMatchedTracks + nFakeTracks)});
  m_histograms.at(key(truthPrefix, "vs_pT"))
      .fill({t_pT, static_cast<double>(nTruthMatchedTracks)});
  m_histograms.at(key(fakePrefix, "vs_pT"))
      .fill({t_pT, static_cast<double>(nFakeTracks)});

  m_histograms.at(key(recoPrefix, "vs_eta"))
      .fill({t_eta, static_cast<double>(nTruthMatchedTracks + nFakeTracks)});
  m_histograms.at(key(truthPrefix, "vs_eta"))
      .fill({t_eta, static_cast<double>(nTruthMatchedTracks)});
  m_histograms.at(key(fakePrefix, "vs_eta"))
      .fill({t_eta, static_cast<double>(nFakeTracks)});
}

}  // namespace ActsExamples
