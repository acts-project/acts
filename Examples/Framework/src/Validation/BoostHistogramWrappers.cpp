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

}  // namespace ActsExamples
