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

#include <TEfficiency.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace ActsExamples {

struct EffPlotTool::Impl {
  /// Tracking efficiency vs eta
  std::unique_ptr<TEfficiency> trackEff_vs_eta{nullptr};
  /// Tracking efficiency vs phi
  std::unique_ptr<TEfficiency> trackEff_vs_phi{nullptr};
  /// Tracking efficiency vs pT
  std::unique_ptr<TEfficiency> trackEff_vs_pT{nullptr};
  /// Tracking efficiency vs log pT
  std::unique_ptr<TEfficiency> trackEff_vs_LogPt{nullptr};
  /// Tracking efficiency vs low pT
  std::unique_ptr<TEfficiency> trackEff_vs_LowPt{nullptr};
  /// Tracking efficiency vs d0
  std::unique_ptr<TEfficiency> trackEff_vs_d0{nullptr};
  /// Tracking efficiency vs z0
  std::unique_ptr<TEfficiency> trackEff_vs_z0{nullptr};
  /// Tracking efficiency vs distance to the closest truth particle
  std::unique_ptr<TEfficiency> trackEff_vs_DeltaR{nullptr};
  /// Tracking efficiency vs production radius
  std::unique_ptr<TEfficiency> trackEff_vs_prodR{nullptr};

  /// Tracking efficiency vs eta and phi
  std::unique_ptr<TEfficiency> trackEff_vs_eta_phi{nullptr};
  /// Tracking efficiency vs eta and pT
  std::unique_ptr<TEfficiency> trackEff_vs_eta_pt{nullptr};

  /// Tracking efficiency vs eta in different pT ranges
  std::vector<std::unique_ptr<TEfficiency>> trackEff_vs_eta_inPtRanges;
  /// Tracking efficiency vs pT in different abs(eta) ranges
  std::vector<std::unique_ptr<TEfficiency>> trackEff_vs_pT_inAbsEtaRanges;
};

EffPlotTool::EffPlotTool(const EffPlotTool::Config& cfg,
                         Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("EffPlotTool", lvl)),
      m_impl(std::make_unique<Impl>()) {}

EffPlotTool::~EffPlotTool() = default;

void EffPlotTool::book() {
  Impl& cache = *m_impl;
  const PlotHelpers::Binning& bEta = m_cfg.varBinning.at("Eta");
  const PlotHelpers::Binning& bPhi = m_cfg.varBinning.at("Phi");
  const PlotHelpers::Binning& bPt = m_cfg.varBinning.at("Pt");
  const PlotHelpers::Binning& bLogPt = m_cfg.varBinning.at("LogPt");
  const PlotHelpers::Binning& bLowPt = m_cfg.varBinning.at("LowPt");
  const PlotHelpers::Binning& bD0 = m_cfg.varBinning.at("D0");
  const PlotHelpers::Binning& bZ0 = m_cfg.varBinning.at("Z0");
  const PlotHelpers::Binning& bDeltaR = m_cfg.varBinning.at("DeltaR");
  const PlotHelpers::Binning& bProdR = m_cfg.varBinning.at("prodR");

  ACTS_DEBUG("Initialize the histograms for efficiency plots");

  const std::string ptCutStr =
      std::format("pT > {} GeV/c", m_cfg.minTruthPt / Acts::UnitConstants::GeV);

  // efficiency vs eta
  cache.trackEff_vs_eta = PlotHelpers::bookEff(
      "trackeff_vs_eta",
      std::format("Tracking efficiency with {};Truth #eta;Efficiency",
                  ptCutStr),
      bEta);
  // efficiency vs phi
  cache.trackEff_vs_phi = PlotHelpers::bookEff(
      "trackeff_vs_phi",
      std::format("Tracking efficiency with {};Truth #phi;Efficiency",
                  ptCutStr),
      bPhi);
  // efficiency vs pT
  cache.trackEff_vs_pT = PlotHelpers::bookEff(
      "trackeff_vs_pT", "Tracking efficiency;Truth pT [GeV/c];Efficiency", bPt);
  // efficiency vs log pT
  cache.trackEff_vs_LogPt = PlotHelpers::bookEff(
      "trackeff_vs_LogPt", "Tracking efficiency;Truth pT [GeV/c];Efficiency",
      bLogPt);
  // efficiency vs low pT
  cache.trackEff_vs_LowPt = PlotHelpers::bookEff(
      "trackeff_vs_LowPt", "Tracking efficiency;Truth pT [GeV/c];Efficiency",
      bLowPt);
  // efficiency vs d0
  cache.trackEff_vs_d0 = PlotHelpers::bookEff(
      "trackeff_vs_d0",
      std::format("Tracking efficiency with {};Truth d_0 [mm];Efficiency",
                  ptCutStr),
      bD0);
  // efficiency vs z0
  cache.trackEff_vs_z0 = PlotHelpers::bookEff(
      "trackeff_vs_z0",
      std::format("Tracking efficiency with {};Truth z_0 [mm];Efficiency",
                  ptCutStr),
      bZ0);
  // efficiancy vs distance to the closest truth particle
  cache.trackEff_vs_DeltaR = PlotHelpers::bookEff(
      "trackeff_vs_DeltaR",
      std::format(
          "Tracking efficiency with {};Closest track #Delta R;Efficiency",
          ptCutStr),
      bDeltaR);
  // efficiency vs production radius
  cache.trackEff_vs_prodR = PlotHelpers::bookEff(
      "trackeff_vs_prodR",
      std::format(
          "Tracking efficiency with {};Production radius [mm];Efficiency",
          ptCutStr),
      bProdR);

  // efficiency vs eta and phi
  cache.trackEff_vs_eta_phi = PlotHelpers::bookEff(
      "trackeff_vs_eta_phi",
      std::format(
          "Tracking efficiency with {};Truth #eta;Truth #phi;Efficiency",
          ptCutStr),
      bEta, bPhi);
  // efficiency vs eta and pT
  cache.trackEff_vs_eta_pt = PlotHelpers::bookEff(
      "trackeff_vs_eta_pt",
      "Tracking efficiency;Truth #eta;Truth pT [GeV/c];Efficiency", bEta, bPt);

  // efficiency vs eta in different pT ranges
  for (const auto& [i, ptRange] : Acts::enumerate(m_cfg.truthPtRangesForEta)) {
    const std::string name = std::format("trackeff_vs_eta_ptRange_{}", i);
    const std::string title = std::format(
        "Tracking efficiency with pT in [{}, {}] GeV/c;Truth #eta;Efficiency",
        ptRange.first / Acts::UnitConstants::GeV,
        ptRange.second / Acts::UnitConstants::GeV);
    cache.trackEff_vs_eta_inPtRanges.push_back(
        PlotHelpers::bookEff(name, title, bEta));
  }
  // efficiency vs pT in different abs(eta) ranges
  for (const auto& [i, absEtaRange] :
       Acts::enumerate(m_cfg.truthAbsEtaRangesForPt)) {
    const std::string name = std::format("trackeff_vs_pT_absEtaRange_{}", i);
    const std::string title = std::format(
        "Tracking efficiency with |#eta| in [{}, {}];Truth pT "
        "[GeV/c];Efficiency",
        absEtaRange.first, absEtaRange.second);
    cache.trackEff_vs_pT_inAbsEtaRanges.push_back(
        PlotHelpers::bookEff(name, title, bPt));
  }
}

void EffPlotTool::write() {
  Impl& cache = *m_impl;
  ACTS_DEBUG("Write the plots to output file.");

  cache.trackEff_vs_eta->Write();
  for (const auto& eff : cache.trackEff_vs_eta_inPtRanges) {
    eff->Write();
  }
  cache.trackEff_vs_eta_phi->Write();
  cache.trackEff_vs_eta_pt->Write();
  cache.trackEff_vs_phi->Write();
  cache.trackEff_vs_pT->Write();
  for (const auto& eff : cache.trackEff_vs_pT_inAbsEtaRanges) {
    eff->Write();
  }
  cache.trackEff_vs_LogPt->Write();
  cache.trackEff_vs_LowPt->Write();
  cache.trackEff_vs_d0->Write();
  cache.trackEff_vs_z0->Write();
  cache.trackEff_vs_DeltaR->Write();
  cache.trackEff_vs_prodR->Write();
}

void EffPlotTool::fill(const Acts::GeometryContext& gctx,
                       const SimParticleState& truthParticle,
                       const double deltaR, const bool status) {
  Impl& cache = *m_impl;
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
    PlotHelpers::fillEff(*cache.trackEff_vs_eta, t_eta, status);
    PlotHelpers::fillEff(*cache.trackEff_vs_phi, t_phi, status);
    PlotHelpers::fillEff(*cache.trackEff_vs_d0, t_d0, status);
    PlotHelpers::fillEff(*cache.trackEff_vs_z0, t_z0, status);
    PlotHelpers::fillEff(*cache.trackEff_vs_DeltaR, t_deltaR, status);
    PlotHelpers::fillEff(*cache.trackEff_vs_prodR, t_prodR, status);

    PlotHelpers::fillEff(*cache.trackEff_vs_eta_phi, t_eta, t_phi, status);
  }

  // do not cut on truth pT as it is a variable on the plot
  PlotHelpers::fillEff(*cache.trackEff_vs_pT, t_pT, status);
  PlotHelpers::fillEff(*cache.trackEff_vs_LogPt, t_pT, status);
  PlotHelpers::fillEff(*cache.trackEff_vs_LowPt, t_pT, status);
  PlotHelpers::fillEff(*cache.trackEff_vs_eta_pt, t_eta, t_pT, status);

  // fill the efficiency vs eta in different pT ranges
  for (const auto& [ptRange, eff] :
       Acts::zip(m_cfg.truthPtRangesForEta, cache.trackEff_vs_eta_inPtRanges)) {
    if (t_pT >= ptRange.first && t_pT < ptRange.second) {
      PlotHelpers::fillEff(*eff, t_eta, status);
    }
  }

  // fill the efficiency vs pT in different eta ranges
  for (const auto& [absEtaRange, eff] : Acts::zip(
           m_cfg.truthAbsEtaRangesForPt, cache.trackEff_vs_pT_inAbsEtaRanges)) {
    if (t_absEta >= absEtaRange.first && t_absEta < absEtaRange.second) {
      PlotHelpers::fillEff(*eff, t_pT, status);
    }
  }
}

}  // namespace ActsExamples
