// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/TrackQualityPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

using namespace Acts::Experimental;

namespace ActsExamples {

namespace {

ProfileHistogram1 makeProfile(const TrackQualityPlotTool::Config& cfg,
                              const std::string& name, const std::string& title,
                              const AxisVariant& ax,
                              const std::string& sampleTitle) {
  const auto& yAxis = cfg.varBinning.at("Num");
  Acts::Range1D<double> yRange{yAxis.bin(0).lower(),
                               yAxis.bin(yAxis.size() - 1).upper()};
  return ProfileHistogram1(name, title, {ax}, sampleTitle, yRange);
}

}  // namespace

TrackQualityPlotTool::TrackQualityPlotTool(const Config& cfg,
                                           Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("TrackQualityPlotTool", lvl)) {
  ACTS_DEBUG("Initialize the histograms for track quality plots");

  m_profiles.insert({"completeness_vs_pT",
                     makeProfile(m_cfg, "completeness_vs_pT", "Completeness",
                                 m_cfg.varBinning.at("Pt"), "Completeness")});
  m_profiles.insert({"completeness_vs_eta",
                     makeProfile(m_cfg, "completeness_vs_eta", "Completeness",
                                 m_cfg.varBinning.at("Eta"), "Completeness")});
  m_profiles.insert({"completeness_vs_phi",
                     makeProfile(m_cfg, "completeness_vs_phi", "Completeness",
                                 m_cfg.varBinning.at("Phi"), "Completeness")});
  m_profiles.insert(
      {"purity_vs_pT", makeProfile(m_cfg, "purity_vs_pT", "Purity",
                                   m_cfg.varBinning.at("Pt"), "Purity")});
  m_profiles.insert(
      {"purity_vs_eta", makeProfile(m_cfg, "purity_vs_eta", "Purity",
                                    m_cfg.varBinning.at("Eta"), "Purity")});
  m_profiles.insert(
      {"purity_vs_phi", makeProfile(m_cfg, "purity_vs_phi", "Purity",
                                    m_cfg.varBinning.at("Phi"), "Purity")});
}

void TrackQualityPlotTool::fill(
    const Acts::BoundTrackParameters& fittedParameters, double completeness,
    double purity) {
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  m_profiles.at("completeness_vs_pT").fill({fit_pT}, completeness);
  m_profiles.at("completeness_vs_eta").fill({fit_eta}, completeness);
  m_profiles.at("completeness_vs_phi").fill({fit_phi}, completeness);

  m_profiles.at("purity_vs_pT").fill({fit_pT}, purity);
  m_profiles.at("purity_vs_eta").fill({fit_eta}, purity);
  m_profiles.at("purity_vs_phi").fill({fit_phi}, purity);
}

}  // namespace ActsExamples
