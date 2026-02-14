// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/DuplicationPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

using namespace Acts::Experimental;

namespace ActsExamples {

namespace {

ProfileHistogram1 makeProfile(const DuplicationPlotTool::Config& cfg,
                              const std::string& name, const std::string& title,
                              const AxisVariant& ax) {
  const auto& yAxis = cfg.varBinning.at("Num");
  Acts::Range1D<double> yRange{yAxis.bin(0).lower(),
                               yAxis.bin(yAxis.size() - 1).upper()};
  return ProfileHistogram1(name, title, {ax}, yAxis.metadata(), yRange);
}

}  // namespace

DuplicationPlotTool::DuplicationPlotTool(const DuplicationPlotTool::Config& cfg,
                                         Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("DuplicationPlotTool", lvl)) {
  ACTS_DEBUG("Initialize the histograms for duplication ratio plots");

  m_profiles.insert(
      {"nDuplicated_vs_pT", makeProfile(m_cfg, "nDuplicated_vs_pT",
                                        "Number of duplicated track candidates",
                                        m_cfg.varBinning.at("Pt"))});
  m_profiles.insert({"nDuplicated_vs_eta",
                     makeProfile(m_cfg, "nDuplicated_vs_eta",
                                 "Number of duplicated track candidates",
                                 m_cfg.varBinning.at("Eta"))});
  m_profiles.insert({"nDuplicated_vs_phi",
                     makeProfile(m_cfg, "nDuplicated_vs_phi",
                                 "Number of duplicated track candidates",
                                 m_cfg.varBinning.at("Phi"))});

  m_efficiencies.insert(
      {"duplicationRatio_vs_pT",
       Efficiency1("duplicationRatio_vs_pT", "Duplication ratio",
                   std::array{m_cfg.varBinning.at("Pt")})});
  m_efficiencies.insert(
      {"duplicationRatio_vs_eta",
       Efficiency1("duplicationRatio_vs_eta", "Duplication ratio",
                   std::array{m_cfg.varBinning.at("Eta")})});
  m_efficiencies.insert(
      {"duplicationRatio_vs_phi",
       Efficiency1("duplicationRatio_vs_phi", "Duplication ratio",
                   std::array{m_cfg.varBinning.at("Phi")})});
}

void DuplicationPlotTool::fill(
    const Acts::BoundTrackParameters& fittedParameters, bool status) {
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  m_efficiencies.at("duplicationRatio_vs_pT").fill({fit_pT}, status);
  m_efficiencies.at("duplicationRatio_vs_eta").fill({fit_eta}, status);
  m_efficiencies.at("duplicationRatio_vs_phi").fill({fit_phi}, status);
}

void DuplicationPlotTool::fill(const SimParticleState& truthParticle,
                               std::size_t nMatchedTracks) {
  const auto t_phi = phi(truthParticle.direction());
  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();

  const auto nDuplicatedTracks = nMatchedTracks == 0 ? 0 : nMatchedTracks - 1;

  m_profiles.at("nDuplicated_vs_pT")
      .fill({t_pT}, static_cast<double>(nDuplicatedTracks));
  m_profiles.at("nDuplicated_vs_eta")
      .fill({t_eta}, static_cast<double>(nDuplicatedTracks));
  m_profiles.at("nDuplicated_vs_phi")
      .fill({t_phi}, static_cast<double>(nDuplicatedTracks));
}

}  // namespace ActsExamples
