// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/TrackQualityPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"

#include <TEfficiency.h>
#include <TProfile.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

namespace ActsExamples {

TrackQualityPlotTool::TrackQualityPlotTool(const Config& cfg,
                                           Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("TrackCompletenessPlotTool", lvl)) {}

void TrackQualityPlotTool::book(Cache& cache) const {
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
  PlotHelpers::Binning bNum = m_cfg.varBinning.at("Num");
  ACTS_DEBUG("Initialize the histograms for completeness plots");

  // completeness vs pT
  cache.completeness_vs_pT = PlotHelpers::bookProf(
      "completeness_vs_pT", "Completeness;pT [GeV/c];Completeness", bPt, bNum);
  // completeness vs eta
  cache.completeness_vs_eta = PlotHelpers::bookProf(
      "completeness_vs_eta", "Completeness;#eta;Completeness", bEta, bNum);
  // completeness vs phi
  cache.completeness_vs_phi = PlotHelpers::bookProf(
      "completeness_vs_phi", "Completeness;#phi;Completeness", bPhi, bNum);

  // purity vs pT
  cache.purity_vs_pT = PlotHelpers::bookProf(
      "purity_vs_pT", "Purity;pT [GeV/c];Purity", bPt, bNum);
  // purity vs eta
  cache.purity_vs_eta =
      PlotHelpers::bookProf("purity_vs_eta", "Purity;#eta;Purity", bEta, bNum);
  // purity vs phi
  cache.purity_vs_phi =
      PlotHelpers::bookProf("purity_vs_phi", "Purity;#phi;Purity", bPhi, bNum);
}

void TrackQualityPlotTool::clear(Cache& cache) const {
  delete cache.completeness_vs_pT;
  delete cache.completeness_vs_eta;
  delete cache.completeness_vs_phi;
  delete cache.purity_vs_pT;
  delete cache.purity_vs_eta;
  delete cache.purity_vs_phi;
}

void TrackQualityPlotTool::write(const Cache& cache) const {
  ACTS_DEBUG("Write the plots to output file.");
  cache.completeness_vs_pT->Write();
  cache.completeness_vs_eta->Write();
  cache.completeness_vs_phi->Write();
  cache.purity_vs_pT->Write();
  cache.purity_vs_eta->Write();
  cache.purity_vs_phi->Write();
}

void TrackQualityPlotTool::fill(
    Cache& cache, const Acts::BoundTrackParameters& fittedParameters,
    double completeness, double purity) const {
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  PlotHelpers::fillProf(cache.completeness_vs_pT, fit_pT, completeness);
  PlotHelpers::fillProf(cache.completeness_vs_eta, fit_eta, completeness);
  PlotHelpers::fillProf(cache.completeness_vs_phi, fit_phi, completeness);

  PlotHelpers::fillProf(cache.purity_vs_pT, fit_pT, purity);
  PlotHelpers::fillProf(cache.purity_vs_eta, fit_eta, purity);
  PlotHelpers::fillProf(cache.purity_vs_phi, fit_phi, purity);
}

}  // namespace ActsExamples
