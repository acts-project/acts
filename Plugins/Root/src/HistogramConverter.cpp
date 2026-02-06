// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/HistogramConverter.hpp"

#include <cmath>
#include <vector>

#include <TEfficiency.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <boost/histogram/accumulators/weighted_sum.hpp>

using namespace Acts::Experimental;

namespace ActsPlugins {

std::unique_ptr<TH1F> toRoot(const Histogram1& boostHist) {
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

std::unique_ptr<TH2F> toRoot(const Histogram2& boostHist) {
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

std::unique_ptr<TProfile> toRoot(const ProfileHistogram1& boostProfile) {
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

std::unique_ptr<TEfficiency> toRoot(const Efficiency1& boostEff) {
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

  // Create TEfficiency from the two histograms
  // TEfficiency takes ownership of the histograms
  auto rootEff = std::make_unique<TEfficiency>(*acceptedHist, *totalHist);
  rootEff->SetName(boostEff.name().c_str());
  rootEff->SetTitle(boostEff.title().c_str());

  return rootEff;
}

std::unique_ptr<TEfficiency> toRoot(const Efficiency2& boostEff) {
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

  // Create TEfficiency from the two histograms
  // TEfficiency takes ownership of the histograms
  auto rootEff = std::make_unique<TEfficiency>(*acceptedHist, *totalHist);
  rootEff->SetName(boostEff.name().c_str());
  rootEff->SetTitle(boostEff.title().c_str());

  return rootEff;
}

}  // namespace ActsPlugins
