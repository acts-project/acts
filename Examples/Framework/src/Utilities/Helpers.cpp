// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/Helpers.hpp"

#include <cassert>

#include <TAxis.h>
#include <TEfficiency.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

namespace ActsExamples {

std::unique_ptr<TH1F> PlotHelpers::bookHisto(const std::string& histName,
                                             const std::string& histTitle,
                                             const Binning& varBinning) {
  auto hist = std::make_unique<TH1F>(histName.c_str(), histTitle.c_str(),
                                     varBinning.nBins(), varBinning.binEdges());
  hist->GetXaxis()->SetTitle(varBinning.title().c_str());
  hist->GetYaxis()->SetTitle("Entries");
  hist->Sumw2();
  return hist;
}

std::unique_ptr<TH2F> PlotHelpers::bookHisto(const std::string& histName,
                                             const std::string& histTitle,
                                             const Binning& varXBinning,
                                             const Binning& varYBinning) {
  auto hist = std::make_unique<TH2F>(
      histName.c_str(), histTitle.c_str(), varXBinning.nBins(),
      varXBinning.binEdges(), varYBinning.nBins(), varYBinning.binEdges());
  hist->GetXaxis()->SetTitle(varXBinning.title().c_str());
  hist->GetYaxis()->SetTitle(varYBinning.title().c_str());
  hist->Sumw2();
  return hist;
}

void PlotHelpers::fillHisto(TH1F& hist, float value, float weight) {
  hist.Fill(value, weight);
}

void PlotHelpers::fillHisto(TH2F& hist, float xValue, float yValue,
                            float weight) {
  hist.Fill(xValue, yValue, weight);
}

void PlotHelpers::anaHisto(TH1D* inputHist, int j, TH1F* meanHist,
                           TH1F* widthHist) {
  // evaluate mean and width via the Gauss fit
  assert(inputHist != nullptr);
  if (inputHist->GetEntries() > 0) {
    TFitResultPtr r = inputHist->Fit("gaus", "QS0");
    if ((r.Get() != nullptr) && ((r->Status() % 1000) == 0)) {
      // fill the mean and width into 'j'th bin of the meanHist and widthHist,
      // respectively
      meanHist->SetBinContent(j, r->Parameter(1));
      meanHist->SetBinError(j, r->ParError(1));
      widthHist->SetBinContent(j, r->Parameter(2));
      widthHist->SetBinError(j, r->ParError(2));
    }
  }
}

std::unique_ptr<TEfficiency> PlotHelpers::bookEff(const std::string& effName,
                                                  const std::string& effTitle,
                                                  const Binning& varBinning) {
  return std::make_unique<TEfficiency>(effName.c_str(), effTitle.c_str(),
                                       varBinning.nBins(),
                                       varBinning.binEdges());
}

std::unique_ptr<TEfficiency> PlotHelpers::bookEff(const std::string& effName,
                                                  const std::string& effTitle,
                                                  const Binning& varXBinning,
                                                  const Binning& varYBinning) {
  return std::make_unique<TEfficiency>(
      effName.c_str(), effTitle.c_str(), varXBinning.nBins(),
      varXBinning.binEdges(), varYBinning.nBins(), varYBinning.binEdges());
}

void PlotHelpers::fillEff(TEfficiency& efficiency, float value, bool status) {
  efficiency.Fill(status, value);
}

void PlotHelpers::fillEff(TEfficiency& efficiency, float xValue, float yValue,
                          bool status) {
  efficiency.Fill(status, xValue, yValue);
}

std::unique_ptr<TProfile> PlotHelpers::bookProf(const std::string& profName,
                                                const std::string& profTitle,
                                                const Binning& varXBinning,
                                                const Binning& varYBinning) {
  auto prof = std::make_unique<TProfile>(
      profName.c_str(), profTitle.c_str(), varXBinning.nBins(),
      varXBinning.binEdges(), varYBinning.low(), varYBinning.high());
  prof->GetXaxis()->SetTitle(varXBinning.title().c_str());
  prof->GetYaxis()->SetTitle(varYBinning.title().c_str());
  return prof;
}

void PlotHelpers::fillProf(TProfile& profile, float xValue, float yValue,
                           float weight) {
  profile.Fill(xValue, yValue, weight);
}

}  // namespace ActsExamples
