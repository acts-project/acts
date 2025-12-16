// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/BoostHistogramToRootConverter.hpp"

#include <TEfficiency.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

#include <boost/histogram/accumulators/weighted_sum.hpp>

#include <cmath>
#include <vector>

namespace ActsExamples::BoostHistogramToRoot {

TH1F* toTH1F(const BoostHistogram1D& boostHist) {
  const auto& bh = boostHist.histogram();
  const auto& axis = bh.axis(0);

  // Extract bin edges from boost histogram axis
  std::vector<double> edges;
  edges.reserve(axis.size() + 1);
  for (std::size_t i = 0; i < axis.size(); ++i) {
    edges.push_back(axis.bin(i).lower());
  }
  edges.push_back(axis.bin(axis.size() - 1).upper());

  // Create ROOT histogram with variable binning
  TH1F* rootHist = new TH1F(boostHist.name().c_str(),
                            boostHist.title().c_str(),
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

TH2F* toTH2F(const BoostHistogram2D& boostHist) {
  const auto& bh = boostHist.histogram();
  const auto& xAxis = bh.axis(0);
  const auto& yAxis = bh.axis(1);

  // Extract bin edges from X axis
  std::vector<double> xEdges;
  xEdges.reserve(xAxis.size() + 1);
  for (std::size_t i = 0; i < xAxis.size(); ++i) {
    xEdges.push_back(xAxis.bin(i).lower());
  }
  xEdges.push_back(xAxis.bin(xAxis.size() - 1).upper());

  // Extract bin edges from Y axis
  std::vector<double> yEdges;
  yEdges.reserve(yAxis.size() + 1);
  for (std::size_t i = 0; i < yAxis.size(); ++i) {
    yEdges.push_back(yAxis.bin(i).lower());
  }
  yEdges.push_back(yAxis.bin(yAxis.size() - 1).upper());

  // Create ROOT histogram with 2D variable binning
  TH2F* rootHist = new TH2F(boostHist.name().c_str(),
                            boostHist.title().c_str(),
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

TProfile* toTProfile(const BoostProfileHistogram& boostProfile) {
  const auto& bh = boostProfile.histogram();
  const auto& axis = bh.axis(0);

  // Extract bin edges from boost histogram axis
  std::vector<double> edges;
  edges.reserve(axis.size() + 1);
  for (std::size_t i = 0; i < axis.size(); ++i) {
    edges.push_back(axis.bin(i).lower());
  }
  edges.push_back(axis.bin(axis.size() - 1).upper());

  // Create ROOT TProfile with variable binning
  // TProfile constructor: (name, title, nbins, edges)
  TProfile* rootProfile = new TProfile(
      boostProfile.name().c_str(), boostProfile.title().c_str(),
      static_cast<int>(axis.size()), edges.data());

  // Copy mean values from boost profile to ROOT profile
  for (auto&& x : boost::histogram::indexed(bh)) {
    const auto& acc = *x;  // Get the accumulator (weighted_mean)

    // ROOT bin numbering starts at 1
    int rootBinIndex = static_cast<int>(x.index(0)) + 1;

    // Get mean value and count from accumulator
    double mean = acc.value();
    double count = acc.count();

    // TProfile stores sum and entries, then computes mean = sum / entries
    // We need to reconstruct the sum from the mean
    if (count > 0) {
      double sum = mean * count;
      rootProfile->SetBinContent(rootBinIndex, sum);
      rootProfile->SetBinEntries(rootBinIndex, count);
    }
  }

  // Set axis titles
  rootProfile->GetXaxis()->SetTitle(boostProfile.xAxisTitle().c_str());
  rootProfile->GetYaxis()->SetTitle(boostProfile.yAxisTitle().c_str());

  return rootProfile;
}

TEfficiency* toTEfficiency(const BoostEfficiency1D& boostEff) {
  const auto& passed = boostEff.passedHistogram();
  const auto& total = boostEff.totalHistogram();
  const auto& axis = passed.axis(0);

  // Extract bin edges
  std::vector<double> edges;
  edges.reserve(axis.size() + 1);
  for (std::size_t i = 0; i < axis.size(); ++i) {
    edges.push_back(axis.bin(i).lower());
  }
  edges.push_back(axis.bin(axis.size() - 1).upper());

  // Create passed and total TH1F histograms
  TH1F* passedHist = new TH1F(
      (boostEff.name() + "_passed").c_str(), "Passed",
      static_cast<int>(axis.size()), edges.data());

  TH1F* totalHist = new TH1F(
      (boostEff.name() + "_total").c_str(), "Total",
      static_cast<int>(axis.size()), edges.data());

  // Fill histograms with counts
  for (std::size_t i = 0; i < axis.size(); ++i) {
    double passedCount = static_cast<double>(passed.at(i));
    double totalCount = static_cast<double>(total.at(i));

    passedHist->SetBinContent(static_cast<int>(i + 1), passedCount);
    totalHist->SetBinContent(static_cast<int>(i + 1), totalCount);
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

TEfficiency* toTEfficiency(const BoostEfficiency2D& boostEff) {
  const auto& passed = boostEff.passedHistogram();
  const auto& total = boostEff.totalHistogram();
  const auto& xAxis = passed.axis(0);
  const auto& yAxis = passed.axis(1);

  // Extract X bin edges
  std::vector<double> xEdges;
  xEdges.reserve(xAxis.size() + 1);
  for (std::size_t i = 0; i < xAxis.size(); ++i) {
    xEdges.push_back(xAxis.bin(i).lower());
  }
  xEdges.push_back(xAxis.bin(xAxis.size() - 1).upper());

  // Extract Y bin edges
  std::vector<double> yEdges;
  yEdges.reserve(yAxis.size() + 1);
  for (std::size_t i = 0; i < yAxis.size(); ++i) {
    yEdges.push_back(yAxis.bin(i).lower());
  }
  yEdges.push_back(yAxis.bin(yAxis.size() - 1).upper());

  // Create passed and total TH2F histograms
  TH2F* passedHist = new TH2F(
      (boostEff.name() + "_passed").c_str(), "Passed",
      static_cast<int>(xAxis.size()), xEdges.data(),
      static_cast<int>(yAxis.size()), yEdges.data());

  TH2F* totalHist = new TH2F(
      (boostEff.name() + "_total").c_str(), "Total",
      static_cast<int>(xAxis.size()), xEdges.data(),
      static_cast<int>(yAxis.size()), yEdges.data());

  // Fill histograms with counts
  for (std::size_t i = 0; i < xAxis.size(); ++i) {
    for (std::size_t j = 0; j < yAxis.size(); ++j) {
      double passedCount = static_cast<double>(passed.at(i, j));
      double totalCount = static_cast<double>(total.at(i, j));

      passedHist->SetBinContent(static_cast<int>(i + 1),
                                static_cast<int>(j + 1), passedCount);
      totalHist->SetBinContent(static_cast<int>(i + 1),
                               static_cast<int>(j + 1), totalCount);
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

}  // namespace ActsExamples::BoostHistogramToRoot
