// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/BoostHistogramToRootConverter.hpp"

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

}  // namespace ActsExamples::BoostHistogramToRoot
