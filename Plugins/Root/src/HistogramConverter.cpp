// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/HistogramConverter.hpp"

#include <vector>

#include <TEfficiency.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <boost/histogram/accumulators/weighted_sum.hpp>

using namespace Acts::Experimental;

std::unique_ptr<TH1F> ActsPlugins::toRoot(const Histogram1& boostHist) {
  const auto& bh = boostHist.histogram();
  const auto& axis = bh.axis(0);

  // Extract bin edges from boost histogram axis
  std::vector<double> edges = extractBinEdges(axis);

  // Create ROOT histogram with variable binning
  auto rootHist = std::make_unique<TH1F>(boostHist.name().c_str(),
                                         boostHist.title().c_str(), axis.size(),
                                         edges.data());

  // Copy bin contents from boost to ROOT
  for (auto&& x : boost::histogram::indexed(bh)) {
    // Dereference to get bin content
    double content = *x;

    // ROOT bin numbering starts at 1 (bin 0 is underflow)
    int rootBinIndex = x.index(0) + 1;
    rootHist->SetBinContent(rootBinIndex, content);
  }

  // Set axis titles from axis metadata
  rootHist->GetXaxis()->SetTitle(axis.metadata().c_str());

  return rootHist;
}

std::unique_ptr<TH2F> ActsPlugins::toRoot(const Histogram2& boostHist) {
  const auto& bh = boostHist.histogram();
  const auto& xAxis = bh.axis(0);
  const auto& yAxis = bh.axis(1);

  // Extract bin edges from X axis
  std::vector<double> xEdges = extractBinEdges(xAxis);

  // Extract bin edges from Y axis
  std::vector<double> yEdges = extractBinEdges(yAxis);

  // Create ROOT histogram with 2D variable binning
  auto rootHist = std::make_unique<TH2F>(
      boostHist.name().c_str(), boostHist.title().c_str(), xAxis.size(),
      xEdges.data(), yAxis.size(), yEdges.data());

  // Copy bin contents from boost to ROOT
  for (auto&& x : boost::histogram::indexed(bh)) {
    // Dereference to get bin content
    double content = *x;

    // ROOT bin numbering starts at 1 (bin 0 is underflow)
    // indexed() gives us 0-based bin indices for each axis
    int rootXBin = x.index(0) + 1;
    int rootYBin = x.index(1) + 1;
    rootHist->SetBinContent(rootXBin, rootYBin, content);
  }

  // Set axis titles from axis metadata
  rootHist->GetXaxis()->SetTitle(xAxis.metadata().c_str());
  rootHist->GetYaxis()->SetTitle(yAxis.metadata().c_str());

  return rootHist;
}

std::unique_ptr<TH3F> ActsPlugins::toRoot(const Histogram3& boostHist) {
  const auto& bh = boostHist.histogram();
  const auto& xAxis = bh.axis(0);
  const auto& yAxis = bh.axis(1);
  const auto& zAxis = bh.axis(2);

  // Extract bin edges from X axis
  std::vector<double> xEdges = extractBinEdges(xAxis);

  // Extract bin edges from Y axis
  std::vector<double> yEdges = extractBinEdges(yAxis);

  // Extract bin edges from Z axis
  std::vector<double> zEdges = extractBinEdges(zAxis);

  // Create ROOT histogram with 3D variable binning
  auto rootHist = std::make_unique<TH3F>(
      boostHist.name().c_str(), boostHist.title().c_str(), xAxis.size(),
      xEdges.data(), yAxis.size(), yEdges.data(), zAxis.size(), zEdges.data());

  // Copy bin contents from boost to ROOT
  for (auto&& x : boost::histogram::indexed(bh)) {
    // Dereference to get bin content
    double content = *x;

    // ROOT bin numbering starts at 1 (bin 0 is underflow)
    // indexed() gives us 0-based bin indices for each axis
    int rootXBin = x.index(0) + 1;
    int rootYBin = x.index(1) + 1;
    int rootZBin = x.index(2) + 1;
    rootHist->SetBinContent(rootXBin, rootYBin, rootZBin, content);
  }

  // Set axis titles from axis metadata
  rootHist->GetXaxis()->SetTitle(xAxis.metadata().c_str());
  rootHist->GetYaxis()->SetTitle(yAxis.metadata().c_str());
  rootHist->GetZaxis()->SetTitle(zAxis.metadata().c_str());

  return rootHist;
}

std::unique_ptr<TProfile> ActsPlugins::toRoot(
    const ProfileHistogram1& boostProfile) {
  const auto& bh = boostProfile.histogram();
  const auto& axis = bh.axis(0);

  // Extract bin edges from boost histogram axis
  std::vector<double> edges = extractBinEdges(axis);

  // Create ROOT TProfile with variable binning
  auto rootProfile = std::make_unique<TProfile>(boostProfile.name().c_str(),
                                                boostProfile.title().c_str(),
                                                axis.size(), edges.data());

  // Enable sum of weights squared storage for proper error calculation
  rootProfile->Sumw2();

  // Copy data from boost profile to ROOT profile
  // - count(): number of fills
  // - value(): mean of y-values
  // - variance(): sample variance = sum((y - mean)^2) / (n - 1)

  using Accumulator = boost::histogram::accumulators::mean<double>;

  for (auto&& x : boost::histogram::indexed(bh)) {
    const Accumulator& acc = *x;

    // ROOT bin numbering starts at 1 (bin 0 is underflow)
    int rootBinIndex = x.index(0) + 1;

    double count = acc.count();

    if (count == 0) {
      continue;
    }
    double mean = acc.value();
    double variance = acc.variance();  // Sample variance from boost

    // ROOT TProfile internally stores the following data (all weights=1):
    // - fArray[bin]      = sum of (w*x)    = count * mean
    // - fBinEntries[bin] = sum of w        = count
    // - fSumw2[bin]      = sum of (w*x)^2  = sum of x^2
    // - fBinSumw2[bin]   = sum of w^2      = count

    double sum = mean * count;

    // To compute sum of y^2 from sample variance s^2:
    // s^2 = Sum((x - μ)^2) / (n-1) = (Sum(x^2) - n*µ^2) / (n-1)
    // Sum(x^2) = (n-1) * s^2 + n*µ^2
    double sumOfSquares = (count - 1.0) * variance + count * mean * mean;

    // Set bin content (sum) and entries via public interface
    rootProfile->SetBinContent(rootBinIndex, sum);
    rootProfile->SetBinEntries(rootBinIndex, count);

    // Access internal arrays to set sum of squares for error calculation
    // This is necessary because ROOT has no public API to set these arrays
    TArrayD* sumw2 = rootProfile->GetSumw2();        // fSumw2 array
    TArrayD* binSumw2 = rootProfile->GetBinSumw2();  // fBinSumw2 array

    assert(sumw2 && "Sumw2 is null");
    assert(binSumw2 && "BinSumw2 is null");

    // Set sum of (weight * value)^2 = sum of y^2 for unweighted data
    sumw2->fArray[rootBinIndex] = sumOfSquares;
    // Set sum of weights^2 = count for unweighted data (all weights = 1)
    binSumw2->fArray[rootBinIndex] = count;
  }

  // Set X axis title from axis metadata
  rootProfile->GetXaxis()->SetTitle(axis.metadata().c_str());
  // Set Y axis title from sampleAxisTitle member
  rootProfile->GetYaxis()->SetTitle(boostProfile.sampleAxisTitle().c_str());

  return rootProfile;
}

std::unique_ptr<TEfficiency> ActsPlugins::toRoot(const Efficiency1& boostEff) {
  const auto& accepted = boostEff.acceptedHistogram();
  const auto& total = boostEff.totalHistogram();
  const auto& axis = accepted.axis(0);

  // Extract bin edges
  std::vector<double> edges = extractBinEdges(axis);

  // Create accepted and total TH1F histograms
  auto acceptedHist =
      std::make_unique<TH1F>((boostEff.name() + "_accepted").c_str(), "Passed",
                             axis.size(), edges.data());

  auto totalHist = std::make_unique<TH1F>((boostEff.name() + "_total").c_str(),
                                          "Total", axis.size(), edges.data());

  // Fill histograms with counts
  for (int i = 0; i < axis.size(); ++i) {
    auto acceptedCount = static_cast<double>(accepted.at(i));
    auto totalCount = static_cast<double>(total.at(i));

    acceptedHist->SetBinContent(i + 1, acceptedCount);
    totalHist->SetBinContent(i + 1, totalCount);
  }

  // Set X axis title from axis metadata
  acceptedHist->GetXaxis()->SetTitle(axis.metadata().c_str());
  totalHist->GetXaxis()->SetTitle(axis.metadata().c_str());

  // Set Y axis titles; use "Efficiency" for total histogram since TEfficiency
  // will inherit this title for the efficiency graph
  totalHist->GetYaxis()->SetTitle("Efficiency");

  // Create TEfficiency from the two histograms
  // TEfficiency takes ownership of the histograms
  auto rootEff = std::make_unique<TEfficiency>(*acceptedHist, *totalHist);
  rootEff->SetName(boostEff.name().c_str());
  rootEff->SetTitle(boostEff.title().c_str());

  return rootEff;
}

std::unique_ptr<TEfficiency> ActsPlugins::toRoot(const Efficiency2& boostEff) {
  const auto& accepted = boostEff.acceptedHistogram();
  const auto& total = boostEff.totalHistogram();
  const auto& xAxis = accepted.axis(0);
  const auto& yAxis = accepted.axis(1);

  // Extract X bin edges
  std::vector<double> xEdges = extractBinEdges(xAxis);

  // Extract Y bin edges
  std::vector<double> yEdges = extractBinEdges(yAxis);

  // Create accepted and total TH2F histograms
  auto acceptedHist = std::make_unique<TH2F>(
      (boostEff.name() + "_accepted").c_str(), "Accepted", xAxis.size(),
      xEdges.data(), yAxis.size(), yEdges.data());

  auto totalHist = std::make_unique<TH2F>((boostEff.name() + "_total").c_str(),
                                          "Total", xAxis.size(), xEdges.data(),
                                          yAxis.size(), yEdges.data());

  // Fill histograms with counts
  for (int i = 0; i < xAxis.size(); ++i) {
    for (int j = 0; j < yAxis.size(); ++j) {
      auto acceptedCount = static_cast<double>(accepted.at(i, j));
      auto totalCount = total.at(i, j);

      acceptedHist->SetBinContent(i + 1, j + 1, acceptedCount);
      totalHist->SetBinContent(i + 1, j + 1, totalCount);
    }
  }

  // Set X axis title from axis metadata
  acceptedHist->GetXaxis()->SetTitle(xAxis.metadata().c_str());
  totalHist->GetXaxis()->SetTitle(xAxis.metadata().c_str());

  // Set Y axis title from axis metadata
  acceptedHist->GetYaxis()->SetTitle(yAxis.metadata().c_str());
  totalHist->GetYaxis()->SetTitle(yAxis.metadata().c_str());

  // Set Z axis titles; use "Efficiency" for total histogram since TEfficiency
  // will inherit this title for the efficiency graph
  totalHist->GetZaxis()->SetTitle("Efficiency");

  // Create TEfficiency from the two histograms
  // TEfficiency takes ownership of the histograms
  auto rootEff = std::make_unique<TEfficiency>(*acceptedHist, *totalHist);
  rootEff->SetName(boostEff.name().c_str());
  rootEff->SetTitle(boostEff.title().c_str());

  return rootEff;
}

std::pair<std::unique_ptr<TH1F>, std::unique_ptr<TH1F>>
ActsPlugins::extractMeanWidth1DProfiles(const TH2F& hist2d,
                                        const std::string& meanName,
                                        const std::string& widthName) {
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
  for (int i = 1; i <= nBinsX; ++i) {
    const auto proj = std::unique_ptr<TH1D>(
        hist2d.ProjectionY(Form("%s_projy_bin%d", hist2d.GetName(), i), i, i));
    if (proj->GetEntries() <= 0) {
      continue;
    }

    TFitResultPtr r = proj->Fit("gaus", "QS0");
    if ((r.Get() == nullptr) || ((r->Status() % 1000) != 0)) {
      continue;
    }

    // Fill mean
    meanHist->SetBinContent(i, r->Parameter(1));
    meanHist->SetBinError(i, r->ParError(1));

    // Fill width (sigma)
    widthHist->SetBinContent(i, r->Parameter(2));
    widthHist->SetBinError(i, r->ParError(2));
  }

  meanHist->GetXaxis()->SetTitle(hist2d.GetXaxis()->GetTitle());
  meanHist->GetYaxis()->SetTitle(hist2d.GetYaxis()->GetTitle());

  widthHist->GetXaxis()->SetTitle(hist2d.GetXaxis()->GetTitle());
  widthHist->GetYaxis()->SetTitle(hist2d.GetYaxis()->GetTitle());

  return {std::move(meanHist), std::move(widthHist)};
}

std::pair<std::unique_ptr<TH2F>, std::unique_ptr<TH2F>>
ActsPlugins::extractMeanWidth2DProfiles(const TH3F& hist3d,
                                        const std::string& meanName,
                                        const std::string& widthName) {
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

  // Loop over all (X,Y) bins
  for (int i = 1; i <= nBinsX; ++i) {
    for (int j = 1; j <= nBinsY; ++j) {
      const auto proj = std::unique_ptr<TH1D>(hist3d.ProjectionZ(
          Form("%s_projz_bin%d_%d", hist3d.GetName(), i, j), i, i, j, j));

      if (proj->GetEntries() <= 0) {
        continue;
      }

      TFitResultPtr r = proj->Fit("gaus", "QS0");
      if ((r.Get() == nullptr) || ((r->Status() % 1000) != 0)) {
        continue;
      }

      // Fill mean
      meanHist->SetBinContent(i, j, r->Parameter(1));
      meanHist->SetBinError(i, j, r->ParError(1));

      // Fill width (sigma)
      widthHist->SetBinContent(i, j, r->Parameter(2));
      widthHist->SetBinError(i, j, r->ParError(2));
    }
  }

  meanHist->GetXaxis()->SetTitle(hist3d.GetXaxis()->GetTitle());
  meanHist->GetYaxis()->SetTitle(hist3d.GetYaxis()->GetTitle());
  meanHist->GetZaxis()->SetTitle(hist3d.GetZaxis()->GetTitle());

  widthHist->GetXaxis()->SetTitle(hist3d.GetXaxis()->GetTitle());
  widthHist->GetYaxis()->SetTitle(hist3d.GetYaxis()->GetTitle());
  widthHist->GetZaxis()->SetTitle(hist3d.GetZaxis()->GetTitle());

  return {std::move(meanHist), std::move(widthHist)};
}
