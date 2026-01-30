// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/EffPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <format>
#include <limits>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace ActsExamples {

EffPlotTool::EffPlotTool(const EffPlotTool::Config& cfg,
                         Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("EffPlotTool", lvl)) {
  ACTS_DEBUG("Initialize the histograms for efficiency plots");

  using Efficiency1 = Acts::Experimental::Efficiency1;
  using Efficiency2 = Acts::Experimental::Efficiency2;

  // 1D efficiencies
  m_efficiencies1D.insert(
      {"trackeff_vs_eta",
       Efficiency1("trackeff_vs_eta",
                   std::format("Tracking efficiency with pT > {} GeV/c",
                               m_cfg.minTruthPt / Acts::UnitConstants::GeV),
                   std::array{m_cfg.varBinning.at("Eta")})});
  m_efficiencies1D.insert(
      {"trackeff_vs_phi",
       Efficiency1("trackeff_vs_phi",
                   std::format("Tracking efficiency with pT > {} GeV/c",
                               m_cfg.minTruthPt / Acts::UnitConstants::GeV),
                   std::array{m_cfg.varBinning.at("Phi")})});
  m_efficiencies1D.insert(
      {"trackeff_vs_pT", Efficiency1("trackeff_vs_pT", "Tracking efficiency",
                                     std::array{m_cfg.varBinning.at("Pt")})});
  m_efficiencies1D.insert(
      {"trackeff_vs_LogPt",
       Efficiency1("trackeff_vs_LogPt", "Tracking efficiency",
                   std::array{m_cfg.varBinning.at("LogPt")})});
  m_efficiencies1D.insert(
      {"trackeff_vs_LowPt",
       Efficiency1("trackeff_vs_LowPt", "Tracking efficiency",
                   std::array{m_cfg.varBinning.at("LowPt")})});
  m_efficiencies1D.insert(
      {"trackeff_vs_d0",
       Efficiency1("trackeff_vs_d0",
                   std::format("Tracking efficiency with pT > {} GeV/c",
                               m_cfg.minTruthPt / Acts::UnitConstants::GeV),
                   std::array{m_cfg.varBinning.at("D0")})});
  m_efficiencies1D.insert(
      {"trackeff_vs_z0",
       Efficiency1("trackeff_vs_z0",
                   std::format("Tracking efficiency with pT > {} GeV/c",
                               m_cfg.minTruthPt / Acts::UnitConstants::GeV),
                   std::array{m_cfg.varBinning.at("Z0")})});
  m_efficiencies1D.insert(
      {"trackeff_vs_DeltaR",
       Efficiency1("trackeff_vs_DeltaR",
                   std::format("Tracking efficiency with pT > {} GeV/c",
                               m_cfg.minTruthPt / Acts::UnitConstants::GeV),
                   std::array{m_cfg.varBinning.at("DeltaR")})});
  m_efficiencies1D.insert(
      {"trackeff_vs_prodR",
       Efficiency1("trackeff_vs_prodR",
                   std::format("Tracking efficiency with pT > {} GeV/c",
                               m_cfg.minTruthPt / Acts::UnitConstants::GeV),
                   std::array{m_cfg.varBinning.at("prodR")})});

  // 2D efficiencies
  m_efficiencies2D.insert(
      {"trackeff_vs_eta_phi",
       Efficiency2("trackeff_vs_eta_phi",
                   std::format("Tracking efficiency with pT > {} GeV/c",
                               m_cfg.minTruthPt / Acts::UnitConstants::GeV),
                   std::array{m_cfg.varBinning.at("Eta"),
                              m_cfg.varBinning.at("Phi")})});
  m_efficiencies2D.insert(
      {"trackeff_vs_eta_pt",
       Efficiency2(
           "trackeff_vs_eta_pt", "Tracking efficiency",
           std::array{m_cfg.varBinning.at("Eta"), m_cfg.varBinning.at("Pt")})});

  const auto& etaAxis = m_cfg.varBinning.at("Eta");
  const auto& ptAxis = m_cfg.varBinning.at("Pt");

  // efficiency vs eta in different pT ranges
  for (const auto& [i, ptRange] : Acts::enumerate(m_cfg.truthPtRangesForEta)) {
    const std::string name = std::format("trackeff_vs_eta_ptRange_{}", i);
    const std::string title =
        std::format("Tracking efficiency with pT in [{}, {}] GeV/c",
                    ptRange.first / Acts::UnitConstants::GeV,
                    ptRange.second / Acts::UnitConstants::GeV);
    m_trackEffVsEtaInPtRanges.emplace_back(name, title, std::array{etaAxis});
  }

  // efficiency vs pT in different abs(eta) ranges
  for (const auto& [i, absEtaRange] :
       Acts::enumerate(m_cfg.truthAbsEtaRangesForPt)) {
    const std::string name = std::format("trackeff_vs_pT_absEtaRange_{}", i);
    const std::string title =
        std::format("Tracking efficiency with |#eta| in [{}, {}]",
                    absEtaRange.first, absEtaRange.second);
    m_trackEffVsPtInAbsEtaRanges.emplace_back(name, title, std::array{ptAxis});
  }
}

void EffPlotTool::fill(const Acts::GeometryContext& gctx,
                       const SimParticleState& truthParticle,
                       const double deltaR, const bool status) {
  constexpr double nan = std::numeric_limits<double>::quiet_NaN();

  const auto intersection =
      m_cfg.beamline
          ->intersect(gctx, truthParticle.position(), truthParticle.direction())
          .closest();
  Acts::Vector2 d0z0{nan, nan};
  if (intersection.isValid()) {
    auto localRes = m_cfg.beamline->globalToLocal(gctx, intersection.position(),
                                                  truthParticle.direction());
    if (localRes.ok()) {
      d0z0 = localRes.value();
    }
  }

  const double t_phi = phi(truthParticle.direction());
  const double t_eta = eta(truthParticle.direction());
  const double t_absEta = std::abs(t_eta);
  const double t_pT = truthParticle.transverseMomentum();
  const double t_d0 = d0z0.x();
  const double t_z0 = d0z0.y();
  const double t_deltaR = deltaR;
  const double t_prodR = perp(truthParticle.position());

  // cut on truth pT with the global range for the relevant plots
  if (t_pT >= m_cfg.minTruthPt) {
    m_efficiencies1D.at("trackeff_vs_eta").fill({t_eta}, status);
    m_efficiencies1D.at("trackeff_vs_phi").fill({t_phi}, status);
    m_efficiencies1D.at("trackeff_vs_d0").fill({t_d0}, status);
    m_efficiencies1D.at("trackeff_vs_z0").fill({t_z0}, status);
    m_efficiencies1D.at("trackeff_vs_DeltaR").fill({t_deltaR}, status);
    m_efficiencies1D.at("trackeff_vs_prodR").fill({t_prodR}, status);

    m_efficiencies2D.at("trackeff_vs_eta_phi").fill({t_eta, t_phi}, status);
  }

  // do not cut on truth pT as it is a variable on the plot
  m_efficiencies1D.at("trackeff_vs_pT").fill({t_pT}, status);
  m_efficiencies1D.at("trackeff_vs_LogPt").fill({t_pT}, status);
  m_efficiencies1D.at("trackeff_vs_LowPt").fill({t_pT}, status);
  m_efficiencies2D.at("trackeff_vs_eta_pt").fill({t_eta, t_pT}, status);

  // fill the efficiency vs eta in different pT ranges
  for (auto&& [ptRange, eff] :
       Acts::zip(m_cfg.truthPtRangesForEta, m_trackEffVsEtaInPtRanges)) {
    if (t_pT >= ptRange.first && t_pT < ptRange.second) {
      eff.fill({t_eta}, status);
    }
  }

  // fill the efficiency vs pT in different eta ranges
  for (auto&& [absEtaRange, eff] :
       Acts::zip(m_cfg.truthAbsEtaRangesForPt, m_trackEffVsPtInAbsEtaRanges)) {
    if (t_absEta >= absEtaRange.first && t_absEta < absEtaRange.second) {
      eff.fill({t_pT}, status);
    }
  }
}

}  // namespace ActsExamples
