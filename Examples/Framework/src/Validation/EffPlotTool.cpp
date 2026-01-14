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
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("EffPlotTool", lvl)) {}

void EffPlotTool::book(Cache& cache) const {
  const auto& etaAxis = m_cfg.varBinning.at("Eta");
  const auto& phiAxis = m_cfg.varBinning.at("Phi");
  const auto& ptAxis = m_cfg.varBinning.at("Pt");
  const auto& logPtAxis = m_cfg.varBinning.at("LogPt");
  const auto& lowPtAxis = m_cfg.varBinning.at("LowPt");
  const auto& d0Axis = m_cfg.varBinning.at("D0");
  const auto& z0Axis = m_cfg.varBinning.at("Z0");
  const auto& deltaRAxis = m_cfg.varBinning.at("DeltaR");
  const auto& prodRAxis = m_cfg.varBinning.at("prodR");

  ACTS_DEBUG("Initialize the histograms for efficiency plots");

  const std::string ptCutStr =
      std::format("pT > {} GeV/c", m_cfg.minTruthPt / Acts::UnitConstants::GeV);

  // efficiency vs eta
  cache.trackEff_vs_eta = Acts::Experimental::Efficiency1(
      "trackeff_vs_eta",
      std::format("Tracking efficiency with {};Truth #eta;Efficiency", ptCutStr),
      std::array{etaAxis});

  // efficiency vs phi
  cache.trackEff_vs_phi = Acts::Experimental::Efficiency1(
      "trackeff_vs_phi",
      std::format("Tracking efficiency with {};Truth #phi;Efficiency", ptCutStr),
      std::array{phiAxis});

  // efficiency vs pT
  cache.trackEff_vs_pT = Acts::Experimental::Efficiency1(
      "trackeff_vs_pT", "Tracking efficiency;Truth pT [GeV/c];Efficiency",
      std::array{ptAxis});

  // efficiency vs log pT
  cache.trackEff_vs_LogPt = Acts::Experimental::Efficiency1(
      "trackeff_vs_LogPt", "Tracking efficiency;Truth pT [GeV/c];Efficiency",
      std::array{logPtAxis});

  // efficiency vs low pT
  cache.trackEff_vs_LowPt = Acts::Experimental::Efficiency1(
      "trackeff_vs_LowPt", "Tracking efficiency;Truth pT [GeV/c];Efficiency",
      std::array{lowPtAxis});

  // efficiency vs d0
  cache.trackEff_vs_d0 = Acts::Experimental::Efficiency1(
      "trackeff_vs_d0",
      std::format("Tracking efficiency with {};Truth d_0 [mm];Efficiency",
                  ptCutStr),
      std::array{d0Axis});

  // efficiency vs z0
  cache.trackEff_vs_z0 = Acts::Experimental::Efficiency1(
      "trackeff_vs_z0",
      std::format("Tracking efficiency with {};Truth z_0 [mm];Efficiency",
                  ptCutStr),
      std::array{z0Axis});

  // efficiency vs distance to the closest truth particle
  cache.trackEff_vs_DeltaR = Acts::Experimental::Efficiency1(
      "trackeff_vs_DeltaR",
      std::format(
          "Tracking efficiency with {};Closest track #Delta R;Efficiency",
          ptCutStr),
      std::array{deltaRAxis});

  // efficiency vs production radius
  cache.trackEff_vs_prodR = Acts::Experimental::Efficiency1(
      "trackeff_vs_prodR",
      std::format(
          "Tracking efficiency with {};Production radius [mm];Efficiency",
          ptCutStr),
      std::array{prodRAxis});

  // efficiency vs eta and phi
  cache.trackEff_vs_eta_phi = Acts::Experimental::Efficiency2(
      "trackeff_vs_eta_phi",
      std::format(
          "Tracking efficiency with {};Truth #eta;Truth #phi;Efficiency",
          ptCutStr),
      std::array{etaAxis, phiAxis});

  // efficiency vs eta and pT
  cache.trackEff_vs_eta_pt = Acts::Experimental::Efficiency2(
      "trackeff_vs_eta_pt",
      "Tracking efficiency;Truth #eta;Truth pT [GeV/c];Efficiency",
      std::array{etaAxis, ptAxis});

  // efficiency vs eta in different pT ranges
  for (const auto& [i, ptRange] : Acts::enumerate(m_cfg.truthPtRangesForEta)) {
    const std::string name = std::format("trackeff_vs_eta_ptRange_{}", i);
    const std::string title = std::format(
        "Tracking efficiency with pT in [{}, {}] GeV/c;Truth #eta;Efficiency",
        ptRange.first / Acts::UnitConstants::GeV,
        ptRange.second / Acts::UnitConstants::GeV);
    cache.trackEff_vs_eta_inPtRanges.emplace_back(name, title,
        std::array{etaAxis});
  }

  // efficiency vs pT in different abs(eta) ranges
  for (const auto& [i, absEtaRange] :
       Acts::enumerate(m_cfg.truthAbsEtaRangesForPt)) {
    const std::string name = std::format("trackeff_vs_pT_absEtaRange_{}", i);
    const std::string title = std::format(
        "Tracking efficiency with |#eta| in [{}, {}];Truth pT "
        "[GeV/c];Efficiency",
        absEtaRange.first, absEtaRange.second);
    cache.trackEff_vs_pT_inAbsEtaRanges.emplace_back(name, title,
        std::array{ptAxis});
  }
}

void EffPlotTool::fill(const Acts::GeometryContext& gctx, Cache& cache,
                       const SimParticleState& truthParticle,
                       const double deltaR, const bool status) const {
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
    cache.trackEff_vs_eta.fill({t_eta}, status);
    cache.trackEff_vs_phi.fill({t_phi}, status);
    cache.trackEff_vs_d0.fill({t_d0}, status);
    cache.trackEff_vs_z0.fill({t_z0}, status);
    cache.trackEff_vs_DeltaR.fill({t_deltaR}, status);
    cache.trackEff_vs_prodR.fill({t_prodR}, status);

    cache.trackEff_vs_eta_phi.fill({t_eta, t_phi}, status);
  }

  // do not cut on truth pT as it is a variable on the plot
  cache.trackEff_vs_pT.fill({t_pT}, status);
  cache.trackEff_vs_LogPt.fill({t_pT}, status);
  cache.trackEff_vs_LowPt.fill({t_pT}, status);
  cache.trackEff_vs_eta_pt.fill({t_eta, t_pT}, status);

  // fill the efficiency vs eta in different pT ranges
  for (const auto& [ptRange, eff] :
       Acts::zip(m_cfg.truthPtRangesForEta, cache.trackEff_vs_eta_inPtRanges)) {
    if (t_pT >= ptRange.first && t_pT < ptRange.second) {
      eff.fill({t_eta}, status);
    }
  }

  // fill the efficiency vs pT in different eta ranges
  for (const auto& [absEtaRange, eff] : Acts::zip(
           m_cfg.truthAbsEtaRangesForPt, cache.trackEff_vs_pT_inAbsEtaRanges)) {
    if (t_absEta >= absEtaRange.first && t_absEta < absEtaRange.second) {
      eff.fill({t_pT}, status);
    }
  }
}

}  // namespace ActsExamples
