// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"

#include <TProfile.h>

namespace ActsExamples {

TrackSummaryPlotTool::TrackSummaryPlotTool(
    const TrackSummaryPlotTool::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("TrackSummaryPlotTool", lvl)) {}

void TrackSummaryPlotTool::book(Cache& cache, const std::string& prefix) const {
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  PlotHelpers::Binning bNum = m_cfg.varBinning.at("Num");
  PlotHelpers::Binning bpmNum = m_cfg.varBinning.at("pmNum");
  ACTS_DEBUG("Initialize the histograms for track info plots, use prefix '"
             << prefix << "'");
  auto addPrefix = [&](const std::string& name) {
    return prefix.empty() ? name : prefix + "_" + name;
  };
  // number of track states versus eta
  cache.nStates_vs_eta =
      PlotHelpers::bookProf(addPrefix("nStates_vs_eta").c_str(),
                            "Number of total states vs. #eta", bEta, bNum);
  // number of measurements versus eta
  cache.nMeasurements_vs_eta =
      PlotHelpers::bookProf(addPrefix("nMeasurements_vs_eta").c_str(),
                            "Number of measurements vs. #eta", bEta, bNum);
  // number of holes versus eta
  cache.nHoles_vs_eta =
      PlotHelpers::bookProf(addPrefix("nHoles_vs_eta").c_str(),
                            "Number of holes vs. #eta", bEta, bNum);
  // number of outliers versus eta
  cache.nOutliers_vs_eta =
      PlotHelpers::bookProf(addPrefix("nOutliers_vs_eta").c_str(),
                            "Number of outliers vs. #eta", bEta, bNum);
  // number of Shared Hits versus eta
  cache.nSharedHits_vs_eta =
      PlotHelpers::bookProf(addPrefix("nSharedHits_vs_eta").c_str(),
                            "Number of Shared Hits vs. #eta", bEta, bNum);
  // number of track states versus pt
  cache.nStates_vs_pt =
      PlotHelpers::bookProf(addPrefix("nStates_vs_pT").c_str(),
                            "Number of total states vs. pT", bPt, bNum);
  // number of measurements versus pt
  cache.nMeasurements_vs_pt =
      PlotHelpers::bookProf(addPrefix("nMeasurements_vs_pT").c_str(),
                            "Number of measurements vs. pT", bPt, bNum);
  // number of holes versus pt
  cache.nHoles_vs_pt = PlotHelpers::bookProf(
      addPrefix("nHoles_vs_pT").c_str(), "Number of holes vs. pT", bPt, bNum);
  // number of outliers versus pt
  cache.nOutliers_vs_pt =
      PlotHelpers::bookProf(addPrefix("nOutliers_vs_pT").c_str(),
                            "Number of outliers vs. pT", bPt, bNum);
  // number of Shared Hits versus pt
  cache.nSharedHits_vs_pt =
      PlotHelpers::bookProf(addPrefix("nSharedHits_vs_pT").c_str(),
                            "Number of Shared Hits vs. pT", bPt, bNum);
       // 2D histograms for holes vs eta
  cache.nHoles_vs_eta_2D = PlotHelpers::bookHisto(
      addPrefix("nHoles_vs_eta_2D").c_str(),
      "Number of holes vs. #eta;#eta;Number of holes", bEta, bNum);
  // 2D histograms for measurements vs eta
  cache.nMeasurements_vs_eta_2D = PlotHelpers::bookHisto(
      addPrefix("nMeasurements_vs_eta_2D").c_str(),
      "Number of measurements vs. #eta;#eta;Number of measurements", bEta, bNum);
  // 2D histograms for outliers vs eta
  cache.nOutliers_vs_eta_2D = PlotHelpers::bookHisto(
      addPrefix("nOutliers_vs_eta_2D").c_str(),
      "Number of outliers vs. #eta;#eta;Number of outliers", bEta, bNum);
  // 2D histograms for Changed Measurements vs eta
  cache.nChangedMeasurements_vs_eta_2D = PlotHelpers::bookHisto(
      addPrefix("nChangedMeasurements_vs_eta_2D").c_str(),
      "Number of Changed Measurements vs. #eta;#eta;Number of Changed Measurements",
      bEta, bpmNum);
  // 2D histograms for nHoles vs pt
  cache.nHoles_vs_pt_2D = PlotHelpers::bookHisto(
      addPrefix("nHoles_vs_pT_2D").c_str(),
      "Number of holes vs. pT;pT;Number of holes", bPt, bNum);
  // 2D histograms for nMeasurements vs pt
  cache.nMeasurements_vs_pt_2D = PlotHelpers::bookHisto(
      addPrefix("nMeasurements_vs_pT_2D").c_str(),
      "Number of measurements vs. pT;pT;Number of measurements", bPt, bNum);
  // 2D histograms for nOutliers vs pt
  cache.nOutliers_vs_pt_2D = PlotHelpers::bookHisto(
      addPrefix("nOutliers_vs_pT_2D").c_str(),
      "Number of outliers vs. pT;pT;Number of outliers", bPt, bNum);
  // 2D histograms for Changed Measurements vs pt
  cache.nChangedMeasurements_vs_pt_2D = PlotHelpers::bookHisto(
      addPrefix("nChangedMeasurements_vs_pT_2D").c_str(),
      "Number of Changed Measurements vs. pT;pT;Number of Changed Measurements",
      bPt, bpmNum);
  // number of max measurements versus eta
  cache.maxMeasurements_vs_eta =
      PlotHelpers::bookProf(addPrefix("nProtoTrackMeasurements_vs_eta").c_str(),
                              "Number of ProtoTrack measurements vs. #eta", bEta, bNum);
  // number of max measurements versus pt
  cache.maxMeasurements_vs_pt =
      PlotHelpers::bookProf(addPrefix("nProtoTrackMeasurements_vs_pT").c_str(),
                              "Number of ProtoTrack measurements vs. pT", bPt, bNum);
  // 2D histograms for max measurements vs eta
  cache.maxMeasurements_vs_eta_2D = PlotHelpers::bookHisto(
      addPrefix("nProtoTrackMeasurements_vs_eta_2D").c_str(),
      "Number of ProtoTrack measurements vs. #eta;#eta;Number of ProtoTrack measurements",
      bEta, bNum);
  // 2D histograms for max measurements vs pt
  cache.maxMeasurements_vs_pt_2D = PlotHelpers::bookHisto(
      addPrefix("nProtoTrackMeasurements_vs_pT_2D").c_str(),
      "Number of ProtoTrack measurements vs. pT;pT;Number of ProtoTrack measurements",
      bPt, bNum);
  cache.nChangedMeasurements_vs_matchingprob = PlotHelpers::bookEff(
      addPrefix("nChangedMeasurements_vs_matchingprob").c_str(),
      "Number of Changed Measurements vs. Matching Probability;Matching Probability;Number of Changed Measurements",
      bpmNum);
  cache.nHoles_vs_matchingprob = PlotHelpers::bookEff(
      addPrefix("nHoles_vs_matchingprob").c_str(),
      "Number of Holes vs. Matching Probability;Matching Probability;Number of Holes",
      bNum);
  cache.nMeasurements_vs_matchingprob = PlotHelpers::bookEff(
      addPrefix("nMeasurements_vs_matchingprob").c_str(),
      "Number of Measurements vs. Matching Probability;Matching Probability;Number of Measurements",
      bNum);
  cache.nOutliers_vs_matchingprob = PlotHelpers::bookEff(
      addPrefix("nOutliers_vs_matchingprob").c_str(),
      "Number of Outliers vs. Matching Probability;Matching Probability;Number of Outliers",
      bNum);
  cache.maxMeasurements_vs_matchingprob = PlotHelpers::bookEff(
      addPrefix("nProtoTrackMeasurements_vs_matchingprob").c_str(),
      "Number of ProtoTrack measurements vs. Matching Probability;Matching Probability;Number of ProtoTrack measurements",
      bNum);
  cache.nChangedMeasurements_vs_eta_2D_matchingprob = PlotHelpers::bookEff(
      addPrefix("nChangedMeasurements_vs_eta_2D_matchingprob").c_str(),
      "Number of Changed Measurements vs. #eta with Matching Probability;#eta;Number of Changed Measurements",
      bEta, bpmNum);
  cache.nChangedMeasurements_vs_pt_2D_matchingprob = PlotHelpers::bookEff(
      addPrefix("nChangedMeasurements_vs_pT_2D_matchingprob").c_str(),
      "Number of Changed Measurements vs. pT with Matching Probability;pT;Number of Changed Measurements",
      bPt, bpmNum);
  cache.maxMeasurements_vs_eta_2D_matchingprob = PlotHelpers::bookEff(
      addPrefix("nProtoTrackMeasurements_vs_eta_2D_matchingprob").c_str(),
      "Number of ProtoTrack measurements vs. #eta with Matching Probability;#eta;Number of ProtoTrack measurements",
      bEta, bNum);
  cache.maxMeasurements_vs_pt_2D_matchingprob = PlotHelpers::bookEff(
      addPrefix("nProtoTrackMeasurements_vs_pT_2D_matchingprob").c_str(),
      "Number of ProtoTrack measurements vs. pT with Matching Probability;pT;Number of ProtoTrack measurements",
      bPt, bNum);

}

void TrackSummaryPlotTool::clear(Cache& cache) const {
  delete cache.nStates_vs_eta;
  delete cache.nMeasurements_vs_eta;
  delete cache.nOutliers_vs_eta;
  delete cache.nHoles_vs_eta;
  delete cache.nSharedHits_vs_eta;
  delete cache.nStates_vs_pt;
  delete cache.nMeasurements_vs_pt;
  delete cache.nOutliers_vs_pt;
  delete cache.nHoles_vs_pt;
  delete cache.nSharedHits_vs_pt;
  delete cache.nHoles_vs_eta_2D;
  delete cache.nMeasurements_vs_eta_2D;
  delete cache.nOutliers_vs_eta_2D;
  delete cache.nChangedMeasurements_vs_eta_2D;
  delete cache.nHoles_vs_pt_2D;
  delete cache.nMeasurements_vs_pt_2D;
  delete cache.nOutliers_vs_pt_2D;
  delete cache.nChangedMeasurements_vs_pt_2D;
  delete cache.maxMeasurements_vs_eta;
  delete cache.maxMeasurements_vs_pt;
  delete cache.maxMeasurements_vs_eta_2D;
  delete cache.maxMeasurements_vs_pt_2D;
  delete cache.nChangedMeasurements_vs_matchingprob;
  delete cache.nHoles_vs_matchingprob;
  delete cache.nMeasurements_vs_matchingprob;
  delete cache.nOutliers_vs_matchingprob;
  delete cache.maxMeasurements_vs_matchingprob;
  delete cache.nChangedMeasurements_vs_eta_2D_matchingprob;
  delete cache.nChangedMeasurements_vs_pt_2D_matchingprob;
  delete cache.maxMeasurements_vs_eta_2D_matchingprob;
  delete cache.maxMeasurements_vs_pt_2D_matchingprob;
}

void TrackSummaryPlotTool::write(const Cache& cache) const {
  ACTS_DEBUG("Write the plots to output file.");
  cache.nStates_vs_eta->Write();
  cache.nMeasurements_vs_eta->Write();
  cache.nOutliers_vs_eta->Write();
  cache.nHoles_vs_eta->Write();
  cache.nSharedHits_vs_eta->Write();
  cache.nStates_vs_pt->Write();
  cache.nMeasurements_vs_pt->Write();
  cache.nOutliers_vs_pt->Write();
  cache.nHoles_vs_pt->Write();
  cache.nSharedHits_vs_pt->Write();
  cache.nHoles_vs_eta_2D->Write();
  cache.nMeasurements_vs_eta_2D->Write();
  cache.nOutliers_vs_eta_2D->Write();
  cache.nChangedMeasurements_vs_eta_2D->Write();
  cache.nHoles_vs_pt_2D->Write();
  cache.nMeasurements_vs_pt_2D->Write();
  cache.nOutliers_vs_pt_2D->Write();
  cache.nChangedMeasurements_vs_pt_2D->Write();
  cache.maxMeasurements_vs_eta->Write();
  cache.maxMeasurements_vs_pt->Write();
  cache.maxMeasurements_vs_eta_2D->Write();
  cache.maxMeasurements_vs_pt_2D->Write();
  cache.nChangedMeasurements_vs_matchingprob->Write();
  cache.nHoles_vs_matchingprob->Write();
  cache.nMeasurements_vs_matchingprob->Write();
  cache.nOutliers_vs_matchingprob->Write();
  cache.maxMeasurements_vs_matchingprob->Write();
  cache.nChangedMeasurements_vs_eta_2D_matchingprob->Write();
  cache.nChangedMeasurements_vs_pt_2D_matchingprob->Write();
  cache.maxMeasurements_vs_eta_2D_matchingprob->Write();
  cache.maxMeasurements_vs_pt_2D_matchingprob->Write();
}

void TrackSummaryPlotTool::fill(
    Cache& cache, const Acts::BoundTrackParameters& fittedParameters,
    std::size_t nStates, std::size_t nMeasurements, std::size_t nOutliers,
    std::size_t nHoles, std::size_t nSharedHits, std::ptrdiff_t nChangedMeasurements, std::size_t maxMeasurements, bool matched) const {
  using Acts::VectorHelpers::eta;
  using Acts::VectorHelpers::perp;
  const auto momentum = fittedParameters.momentum();
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  PlotHelpers::fillProf(cache.nStates_vs_eta, fit_eta, nStates);
  PlotHelpers::fillProf(cache.nMeasurements_vs_eta, fit_eta, nMeasurements);
  PlotHelpers::fillProf(cache.nOutliers_vs_eta, fit_eta, nOutliers);
  PlotHelpers::fillProf(cache.nHoles_vs_eta, fit_eta, nHoles);
  PlotHelpers::fillProf(cache.nSharedHits_vs_eta, fit_eta, nSharedHits);

  PlotHelpers::fillProf(cache.nStates_vs_pt, fit_pT, nStates);
  PlotHelpers::fillProf(cache.nMeasurements_vs_pt, fit_pT, nMeasurements);
  PlotHelpers::fillProf(cache.nOutliers_vs_pt, fit_pT, nOutliers);
  PlotHelpers::fillProf(cache.nHoles_vs_pt, fit_pT, nHoles);
  PlotHelpers::fillProf(cache.nSharedHits_vs_pt, fit_pT, nSharedHits);

  PlotHelpers::fillHisto(cache.nHoles_vs_eta_2D, fit_eta,
                         nHoles);
  PlotHelpers::fillHisto(cache.nMeasurements_vs_eta_2D,
                         fit_eta, nMeasurements);
  PlotHelpers::fillHisto(cache.nOutliers_vs_eta_2D, fit_eta,
                         nOutliers);
  PlotHelpers::fillHisto(cache.nChangedMeasurements_vs_eta_2D,
                         fit_eta, nChangedMeasurements);
  PlotHelpers::fillHisto(cache.nHoles_vs_pt_2D, fit_pT, nHoles);
  PlotHelpers::fillHisto(cache.nMeasurements_vs_pt_2D,
                         fit_pT, nMeasurements);
  PlotHelpers::fillHisto(cache.nOutliers_vs_pt_2D, fit_pT,
                         nOutliers);
  PlotHelpers::fillHisto(cache.nChangedMeasurements_vs_pt_2D,
                         fit_pT, nChangedMeasurements);
  PlotHelpers::fillProf(cache.maxMeasurements_vs_eta, fit_eta, maxMeasurements);
  PlotHelpers::fillProf(cache.maxMeasurements_vs_pt, fit_pT, maxMeasurements);
  PlotHelpers::fillHisto(cache.maxMeasurements_vs_eta_2D, fit_eta,
                            maxMeasurements);
  PlotHelpers::fillHisto(cache.maxMeasurements_vs_pt_2D, fit_pT,
                            maxMeasurements);
  PlotHelpers::fillEff(cache.nChangedMeasurements_vs_matchingprob,
                        nChangedMeasurements, matched);
  PlotHelpers::fillEff(cache.nHoles_vs_matchingprob,
                       nHoles, matched);
  PlotHelpers::fillEff(cache.nMeasurements_vs_matchingprob,
                       nMeasurements, matched);
  PlotHelpers::fillEff(cache.nOutliers_vs_matchingprob,
                       nOutliers, matched);
  PlotHelpers::fillEff(cache.maxMeasurements_vs_matchingprob,
                       maxMeasurements, matched);
  PlotHelpers::fillEff(cache.nChangedMeasurements_vs_eta_2D_matchingprob,
                         fit_eta, nChangedMeasurements, matched);
  PlotHelpers::fillEff(cache.nChangedMeasurements_vs_pt_2D_matchingprob,
                       fit_pT, nChangedMeasurements, matched);
  PlotHelpers::fillEff(cache.maxMeasurements_vs_eta_2D_matchingprob,
                       fit_eta, maxMeasurements, matched);
  PlotHelpers::fillEff(cache.maxMeasurements_vs_pt_2D_matchingprob,
                       fit_pT, maxMeasurements, matched);

}

}  // namespace ActsExamples
