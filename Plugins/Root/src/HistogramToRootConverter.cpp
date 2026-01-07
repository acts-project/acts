// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/HistogramToRootConverter.hpp"

#include <cmath>
#include <vector>

#include <TEfficiency.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <boost/histogram/accumulators/weighted_sum.hpp>

using namespace Acts;

namespace ActsPlugins {

TH1F* toRoot(const Histogram1D& boostHist) {
  const auto& bh = boostHist.histogram();
  const auto& axis = bh.axis(0);

  // Extract bin edges from boost histogram axis
  std::vector<double> edges;
  edges.reserve(axis.size() + 1);
  for (int i = 0; i < axis.size(); ++i) {
    edges.push_back(axis.bin(i).lower());
  }
  edges.push_back(axis.bin(axis.size() - 1).upper());

  // Create ROOT histogram with variable binning
  TH1F* rootHist = new TH1F(boostHist.name().c_str(), boostHist.title().c_str(),
                            static_cast<int>(axis.size()), edges.data());

  // Copy bin contents from boost to ROOT
  for (auto&& x : boost::histogram::indexed(bh)) {
    // Dereference to get bin content
    double content = *x;

    // ROOT bin numbering starts at 1 (bin 0 is underflow)
    int rootBinIndex = static_cast<int>(x.index(0)) + 1;
    rootHist->SetBinContent(rootBinIndex, content);
  }

  // Set axis titles
  rootHist->GetXaxis()->SetTitle(boostHist.axisTitle().c_str());

  return rootHist;
}

TH2F* toRoot(const Histogram2D& boostHist) {
  const auto& bh = boostHist.histogram();
  const auto& xAxis = bh.axis(0);
  const auto& yAxis = bh.axis(1);

  // Extract bin edges from X axis
  std::vector<double> xEdges;
  xEdges.reserve(xAxis.size() + 1);
  for (int i = 0; i < xAxis.size(); ++i) {
    xEdges.push_back(xAxis.bin(i).lower());
  }
  xEdges.push_back(xAxis.bin(xAxis.size() - 1).upper());

  // Extract bin edges from Y axis
  std::vector<double> yEdges;
  yEdges.reserve(yAxis.size() + 1);
  for (int i = 0; i < yAxis.size(); ++i) {
    yEdges.push_back(yAxis.bin(i).lower());
  }
  yEdges.push_back(yAxis.bin(yAxis.size() - 1).upper());

  // Create ROOT histogram with 2D variable binning
  TH2F* rootHist = new TH2F(boostHist.name().c_str(), boostHist.title().c_str(),
                            static_cast<int>(xAxis.size()), xEdges.data(),
                            static_cast<int>(yAxis.size()), yEdges.data());

  // Copy bin contents from boost to ROOT
  for (auto&& x : boost::histogram::indexed(bh)) {
    // Dereference to get bin content
    double content = *x;

    // ROOT bin numbering starts at 1 (bin 0 is underflow)
    // indexed() gives us 0-based bin indices for each axis
    int rootXBin = static_cast<int>(x.index(0)) + 1;
    int rootYBin = static_cast<int>(x.index(1)) + 1;
    rootHist->SetBinContent(rootXBin, rootYBin, content);
  }

  // Set axis titles
  rootHist->GetXaxis()->SetTitle(boostHist.xAxisTitle().c_str());
  rootHist->GetYaxis()->SetTitle(boostHist.yAxisTitle().c_str());

  return rootHist;
}

TProfile* toRoot(const ProfileHistogram& boostProfile) {
  const auto& bh = boostProfile.histogram();
  const auto& axis = bh.axis(0);

  // Extract bin edges from boost histogram axis
  std::vector<double> edges;
  edges.reserve(axis.size() + 1);
  for (int i = 0; i < axis.size(); ++i) {
    edges.push_back(axis.bin(i).lower());
  }
  edges.push_back(axis.bin(axis.size() - 1).upper());

  // Create ROOT TProfile with variable binning
  TProfile* rootProfile =
      new TProfile(boostProfile.name().c_str(), boostProfile.title().c_str(),
                   static_cast<int>(axis.size()), edges.data());

  // Enable sum of weights squared storage for proper error calculation
  rootProfile->Sumw2();

  // Copy data from boost profile to ROOT profile
  //
  // The boost histogram uses a weighted_mean accumulator that tracks:
  // - count(): number of fills
  // - value(): mean of y-values
  // - variance(): sample variance = sum((y - mean)^2) / (n - 1)
  //
  // ROOT TProfile needs different information for correct error calculation.
  // We must reconstruct ROOT's internal arrays from the boost accumulator data.
  using Accumulator = boost::histogram::accumulators::weighted_mean<double>;

  for (auto&& x : boost::histogram::indexed(bh)) {
    const Accumulator& acc = *x;

    // ROOT bin numbering starts at 1 (bin 0 is underflow)
    int rootBinIndex = static_cast<int>(x.index(0)) + 1;

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
    double sum_of_squares = (count - 1.0) * variance + count * mean * mean;

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
    sumw2->fArray[rootBinIndex] = sum_of_squares;
    // Set sum of weights^2 = count for unweighted data (all weights = 1)
    binSumw2->fArray[rootBinIndex] = count;
  }

  // Set axis titles
  rootProfile->GetXaxis()->SetTitle(boostProfile.xAxisTitle().c_str());
  rootProfile->GetYaxis()->SetTitle(boostProfile.yAxisTitle().c_str());

  return rootProfile;
}

TEfficiency* toRoot(const Efficiency1D& boostEff) {
  const auto& passed = boostEff.passedHistogram();
  const auto& total = boostEff.totalHistogram();
  const auto& axis = passed.axis(0);

  // Extract bin edges
  std::vector<double> edges;
  edges.reserve(axis.size() + 1);
  for (int i = 0; i < axis.size(); ++i) {
    edges.push_back(axis.bin(i).lower());
  }
  edges.push_back(axis.bin(axis.size() - 1).upper());

  // Create passed and total TH1F histograms
  TH1F* passedHist = new TH1F((boostEff.name() + "_passed").c_str(), "Passed",
                              static_cast<int>(axis.size()), edges.data());

  TH1F* totalHist = new TH1F((boostEff.name() + "_total").c_str(), "Total",
                             static_cast<int>(axis.size()), edges.data());

  // Fill histograms with counts
  for (int i = 0; i < axis.size(); ++i) {
    double passedCount = static_cast<double>(passed.at(i));
    double totalCount = static_cast<double>(total.at(i));

    passedHist->SetBinContent(i + 1, passedCount);
    totalHist->SetBinContent(i + 1, totalCount);
  }

  // Create TEfficiency from the two histograms
  // TEfficiency takes ownership of the histograms
  TEfficiency* rootEff = new TEfficiency(*passedHist, *totalHist);
  rootEff->SetName(boostEff.name().c_str());
  rootEff->SetTitle(boostEff.title().c_str());

  // Clean up temporary histograms (TEfficiency made copies)
  delete passedHist;
  delete totalHist;

  return rootEff;
}

TEfficiency* toRoot(const Efficiency2D& boostEff) {
  const auto& passed = boostEff.passedHistogram();
  const auto& total = boostEff.totalHistogram();
  const auto& xAxis = passed.axis(0);
  const auto& yAxis = passed.axis(1);

  // Extract X bin edges
  std::vector<double> xEdges;
  xEdges.reserve(xAxis.size() + 1);
  for (int i = 0; i < xAxis.size(); ++i) {
    xEdges.push_back(xAxis.bin(i).lower());
  }
  xEdges.push_back(xAxis.bin(xAxis.size() - 1).upper());

  // Extract Y bin edges
  std::vector<double> yEdges;
  yEdges.reserve(yAxis.size() + 1);
  for (int i = 0; i < yAxis.size(); ++i) {
    yEdges.push_back(yAxis.bin(i).lower());
  }
  yEdges.push_back(yAxis.bin(yAxis.size() - 1).upper());

  // Create passed and total TH2F histograms
  TH2F* passedHist = new TH2F((boostEff.name() + "_passed").c_str(), "Passed",
                              static_cast<int>(xAxis.size()), xEdges.data(),
                              static_cast<int>(yAxis.size()), yEdges.data());

  TH2F* totalHist = new TH2F((boostEff.name() + "_total").c_str(), "Total",
                             static_cast<int>(xAxis.size()), xEdges.data(),
                             static_cast<int>(yAxis.size()), yEdges.data());

  // Fill histograms with counts
  for (int i = 0; i < xAxis.size(); ++i) {
    for (int j = 0; j < yAxis.size(); ++j) {
      double passedCount = static_cast<double>(passed.at(i, j));
      double totalCount = static_cast<double>(total.at(i, j));

      passedHist->SetBinContent(i + 1, j + 1, passedCount);
      totalHist->SetBinContent(i + 1, j + 1, totalCount);
    }
  }

  // Create TEfficiency from the two histograms
  // TEfficiency takes ownership of the histograms
  TEfficiency* rootEff = new TEfficiency(*passedHist, *totalHist);
  rootEff->SetName(boostEff.name().c_str());
  rootEff->SetTitle(boostEff.title().c_str());

  // Clean up temporary histograms (TEfficiency made copies)
  delete passedHist;
  delete totalHist;

  return rootEff;
}

}  // namespace ActsPlugins
