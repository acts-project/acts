// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/BoostHistogramWrappers.hpp"

#include <vector>

namespace ActsExamples {

BoostHistogram1D::BoostHistogram1D(std::string name, std::string title,
                                   const PlotHelpers::Binning& binning)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_axisTitle(binning.title()) {
  // Create variable axis from bin edges
  std::vector<double> edges(binning.binEdges(),
                            binning.binEdges() + binning.nBins() + 1);

  // Construct boost histogram with variable axis and unlimited storage
  m_hist = boost::histogram::make_histogram(
      boost::histogram::axis::variable<double>(edges));
}

void BoostHistogram1D::fill(double value) {
  m_hist(value);
}

BoostHistogram2D::BoostHistogram2D(std::string name, std::string title,
                                   const PlotHelpers::Binning& xBinning,
                                   const PlotHelpers::Binning& yBinning)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_xAxisTitle(xBinning.title()),
      m_yAxisTitle(yBinning.title()) {
  // Create variable axes from bin edges
  std::vector<double> xEdges(xBinning.binEdges(),
                             xBinning.binEdges() + xBinning.nBins() + 1);
  std::vector<double> yEdges(yBinning.binEdges(),
                             yBinning.binEdges() + yBinning.nBins() + 1);

  // Construct boost histogram with two variable axes and unlimited storage
  m_hist = boost::histogram::make_histogram(
      boost::histogram::axis::variable<double>(xEdges),
      boost::histogram::axis::variable<double>(yEdges));
}

void BoostHistogram2D::fill(double xValue, double yValue) {
  m_hist(xValue, yValue);
}

BoostProfileHistogram::BoostProfileHistogram(
    std::string name, std::string title,
    const PlotHelpers::Binning& xBinning, std::string yAxisTitle)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_xAxisTitle(xBinning.title()),
      m_yAxisTitle(std::move(yAxisTitle)) {
  // Create variable axis from bin edges
  std::vector<double> edges(xBinning.binEdges(),
                            xBinning.binEdges() + xBinning.nBins() + 1);

  // Create profile histogram using make_profile
  m_hist = boost::histogram::make_profile(
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
  // Create variable axis from bin edges
  std::vector<double> edges(binning.binEdges(),
                            binning.binEdges() + binning.nBins() + 1);

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
  // Create variable axes from bin edges
  std::vector<double> xEdges(xBinning.binEdges(),
                             xBinning.binEdges() + xBinning.nBins() + 1);
  std::vector<double> yEdges(yBinning.binEdges(),
                             yBinning.binEdges() + yBinning.nBins() + 1);

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
