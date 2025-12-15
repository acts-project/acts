// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/BoostProfileEfficiency.hpp"

namespace ActsExamples {

BoostProfileHistogram::BoostProfileHistogram(
    std::string name, std::string title,
    const PlotHelpers::Binning& xBinning, std::string yAxisTitle)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_xAxisTitle(xBinning.title()),
      m_yAxisTitle(std::move(yAxisTitle)) {
  // Convert PlotHelpers::Binning to std::vector<double> for axis
  std::vector<double> edges(xBinning.nBins() + 1);
  for (std::size_t i = 0; i <= xBinning.nBins(); ++i) {
    edges[i] = xBinning.binEdges()[i];
  }

  // Create boost histogram with variable axis and weighted_mean accumulator
  m_hist = boost::histogram::make_histogram_with(
      boost::histogram::dense_storage<
          boost::histogram::accumulators::weighted_mean<>>(),
      boost::histogram::axis::variable<double>(edges));
}

void BoostProfileHistogram::fill(double xValue, double yValue) {
  m_hist(xValue, boost::histogram::sample(yValue));
}

BoostEfficiency1D::BoostEfficiency1D(std::string name, std::string title,
                                     const PlotHelpers::Binning& binning)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_axisTitle(binning.title()) {
  // Convert binning to edge vector
  std::vector<double> edges(binning.nBins() + 1);
  for (std::size_t i = 0; i <= binning.nBins(); ++i) {
    edges[i] = binning.binEdges()[i];
  }

  // Create two histograms: one for passed, one for total
  auto axis = boost::histogram::axis::variable<double>(edges);
  m_passed = boost::histogram::make_histogram(axis);
  m_total = boost::histogram::make_histogram(axis);
}

void BoostEfficiency1D::fill(double value, bool passed) {
  m_total(value);
  if (passed) {
    m_passed(value);
  }
}

BoostEfficiency2D::BoostEfficiency2D(std::string name, std::string title,
                                     const PlotHelpers::Binning& xBinning,
                                     const PlotHelpers::Binning& yBinning)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_xAxisTitle(xBinning.title()),
      m_yAxisTitle(yBinning.title()) {
  // Convert binning to edge vectors
  std::vector<double> xEdges(xBinning.nBins() + 1);
  for (std::size_t i = 0; i <= xBinning.nBins(); ++i) {
    xEdges[i] = xBinning.binEdges()[i];
  }

  std::vector<double> yEdges(yBinning.nBins() + 1);
  for (std::size_t i = 0; i <= yBinning.nBins(); ++i) {
    yEdges[i] = yBinning.binEdges()[i];
  }

  // Create two 2D histograms: one for passed, one for total
  auto xAxis = boost::histogram::axis::variable<double>(xEdges);
  auto yAxis = boost::histogram::axis::variable<double>(yEdges);

  m_passed = boost::histogram::make_histogram(xAxis, yAxis);
  m_total = boost::histogram::make_histogram(xAxis, yAxis);
}

void BoostEfficiency2D::fill(double xValue, double yValue, bool passed) {
  m_total(xValue, yValue);
  if (passed) {
    m_passed(xValue, yValue);
  }
}

}  // namespace ActsExamples
