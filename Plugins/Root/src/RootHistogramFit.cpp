// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/RootHistogramFit.hpp"

#include "Acts/Utilities/detail/IterativeGaussianFit.hpp"

#include <cmath>
#include <format>

#include <TFitResult.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

using Acts::Experimental::GaussianFitResult;

namespace ActsPlugins {

std::optional<GaussianFitResult> RootHistogramFit::fit(TH1& hist) const {
  const TFitResultPtr result = hist.Fit("gaus", "SQ0", nullptr);
  if (result.Get() == nullptr || result->Status() % 1000 != 0) {
    return std::nullopt;
  }

  return GaussianFitResult{result->Parameter(1), result->Parameter(2),
                           result->ParError(1), result->ParError(2)};
}

std::optional<GaussianFitResult> RootHistogramFit::fit(TH1& hist, double xMin,
                                                       double xMax) const {
  const TFitResultPtr result = hist.Fit("gaus", "SRQ0", nullptr, xMin, xMax);
  if (result.Get() == nullptr || result->Status() % 1000 != 0) {
    return std::nullopt;
  }

  return GaussianFitResult{result->Parameter(1), result->Parameter(2),
                           result->ParError(1), result->ParError(2)};
}

std::optional<GaussianFitResult> RootHistogramFit::iterativeFit(
    TH1& hist, double sigmaRange, int iterations,
    const Acts::Logger& logger) const {
  return Acts::detail::iterativeGaussianFit(
      [&](double xMin, double xMax) {
        return (std::isinf(xMin) || std::isinf(xMax)) ? fit(hist)
                                                      : fit(hist, xMin, xMax);
      },
      sigmaRange, iterations, logger, hist.GetName());
}

std::tuple<std::unique_ptr<TH1F>, std::unique_ptr<TH1F>, double>
RootHistogramFit::extractMeanWidthProfiles(const TH2F& hist2d,
                                           const std::string& meanName,
                                           const std::string& widthName,
                                           const int minEntriesForFit,
                                           const double sigmaRange,
                                           const int iterations,
                                           const Acts::Logger& logger) const {
  const int nBinsX = hist2d.GetNbinsX();

  // Create mean and width histograms with same X binning as the 2D histogram
  auto meanHist = std::make_unique<TH1F>(
      meanName.c_str(), (std::string(hist2d.GetTitle()) + " mean").c_str(),
      nBinsX, hist2d.GetXaxis()->GetXmin(), hist2d.GetXaxis()->GetXmax());
  auto widthHist = std::make_unique<TH1F>(
      widthName.c_str(), (std::string(hist2d.GetTitle()) + " width").c_str(),
      nBinsX, hist2d.GetXaxis()->GetXmin(), hist2d.GetXaxis()->GetXmax());

  // Copy X axis bin edges for variable binning
  if (hist2d.GetXaxis()->GetXbins()->GetSize() > 0) {
    meanHist->SetBins(nBinsX, hist2d.GetXaxis()->GetXbins()->GetArray());
    widthHist->SetBins(nBinsX, hist2d.GetXaxis()->GetXbins()->GetArray());
  }

  // Project each X bin and extract mean/width via Gaussian fit
  const auto trySlice =
      [&](int flat) -> std::optional<std::optional<GaussianFitResult>> {
    const int i = flat + 1;
    const auto proj = std::unique_ptr<TH1D>(hist2d.ProjectionY(
        std::format("{}_projy_bin_{}", hist2d.GetName(), i).c_str(), i, i));

    if (proj->GetEntries() < minEntriesForFit) {
      return std::nullopt;
    }

    return iterativeFit(*proj, sigmaRange, iterations, logger);
  };

  const auto storeResult = [&](int flat, const GaussianFitResult& result) {
    const int i = flat + 1;
    meanHist->SetBinContent(i, result.mean);
    meanHist->SetBinError(i, result.meanError);
    widthHist->SetBinContent(i, result.sigma);
    widthHist->SetBinError(i, result.sigmaError);
  };

  const double fitFailureFraction =
      Acts::detail::extractMeanWidthProfilesImpl(nBinsX, trySlice, storeResult);

  meanHist->GetXaxis()->SetTitle(hist2d.GetXaxis()->GetTitle());
  meanHist->GetYaxis()->SetTitle(hist2d.GetYaxis()->GetTitle());

  widthHist->GetXaxis()->SetTitle(hist2d.GetXaxis()->GetTitle());
  widthHist->GetYaxis()->SetTitle(hist2d.GetYaxis()->GetTitle());

  return {std::move(meanHist), std::move(widthHist), fitFailureFraction};
}

std::tuple<std::unique_ptr<TH2F>, std::unique_ptr<TH2F>, double>
RootHistogramFit::extractMeanWidthProfiles(const TH3F& hist3d,
                                           const std::string& meanName,
                                           const std::string& widthName,
                                           const int minEntriesForFit,
                                           const double sigmaRange,
                                           const int iterations,
                                           const Acts::Logger& logger) const {
  const int nBinsX = hist3d.GetNbinsX();
  const int nBinsY = hist3d.GetNbinsY();

  // Create output histograms with same XY binning as input
  auto meanHist = std::make_unique<TH2F>(
      meanName.c_str(), (std::string(hist3d.GetTitle()) + " mean").c_str(),
      nBinsX, hist3d.GetXaxis()->GetXmin(), hist3d.GetXaxis()->GetXmax(),
      nBinsY, hist3d.GetYaxis()->GetXmin(), hist3d.GetYaxis()->GetXmax());

  auto widthHist = std::make_unique<TH2F>(
      widthName.c_str(), (std::string(hist3d.GetTitle()) + " width").c_str(),
      nBinsX, hist3d.GetXaxis()->GetXmin(), hist3d.GetXaxis()->GetXmax(),
      nBinsY, hist3d.GetYaxis()->GetXmin(), hist3d.GetYaxis()->GetXmax());

  // Copy X and Y axis bin edges for variable binning
  if (hist3d.GetXaxis()->GetXbins()->GetSize() > 0 ||
      hist3d.GetYaxis()->GetXbins()->GetSize() > 0) {
    meanHist->SetBins(nBinsX, hist3d.GetXaxis()->GetXbins()->GetArray(), nBinsY,
                      hist3d.GetYaxis()->GetXbins()->GetArray());
    widthHist->SetBins(nBinsX, hist3d.GetXaxis()->GetXbins()->GetArray(),
                       nBinsY, hist3d.GetYaxis()->GetXbins()->GetArray());
  }

  // Loop over all (X, Y) bins, Y fastest, matching the nested-loop order of
  // the original per-dimension implementation
  const auto trySlice =
      [&](int flat) -> std::optional<std::optional<GaussianFitResult>> {
    const int i = flat / nBinsY + 1;
    const int j = flat % nBinsY + 1;

    const auto proj = std::unique_ptr<TH1D>(hist3d.ProjectionZ(
        std::format("{}_projz_bin_{}_{}", hist3d.GetName(), i, j).c_str(), i, i,
        j, j));

    if (proj->GetEntries() < minEntriesForFit) {
      return std::nullopt;
    }

    return iterativeFit(*proj, sigmaRange, iterations, logger);
  };

  const auto storeResult = [&](int flat, const GaussianFitResult& result) {
    const int i = flat / nBinsY + 1;
    const int j = flat % nBinsY + 1;

    meanHist->SetBinContent(i, j, result.mean);
    meanHist->SetBinError(i, j, result.meanError);
    widthHist->SetBinContent(i, j, result.sigma);
    widthHist->SetBinError(i, j, result.sigmaError);
  };

  const double fitFailureFraction = Acts::detail::extractMeanWidthProfilesImpl(
      nBinsX * nBinsY, trySlice, storeResult);

  meanHist->GetXaxis()->SetTitle(hist3d.GetXaxis()->GetTitle());
  meanHist->GetYaxis()->SetTitle(hist3d.GetYaxis()->GetTitle());
  meanHist->GetZaxis()->SetTitle(hist3d.GetZaxis()->GetTitle());

  widthHist->GetXaxis()->SetTitle(hist3d.GetXaxis()->GetTitle());
  widthHist->GetYaxis()->SetTitle(hist3d.GetYaxis()->GetTitle());
  widthHist->GetZaxis()->SetTitle(hist3d.GetZaxis()->GetTitle());

  return {std::move(meanHist), std::move(widthHist), fitFailureFraction};
}

}  // namespace ActsPlugins
