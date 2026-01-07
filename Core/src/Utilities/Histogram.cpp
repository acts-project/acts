// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Histogram.hpp"

#include <cmath>
#include <vector>

namespace Acts {

HistBinning HistBinning::Uniform(std::string title, std::size_t bins,
                                 double bMin, double bMax) {
  std::vector<double> binEdges(bins + 1);
  const double step = (bMax - bMin) / bins;
  std::generate(binEdges.begin(), binEdges.end(), [&, v = bMin]() mutable {
    const double r = v;
    v += step;
    return r;
  });
  return HistBinning(std::move(title), std::move(binEdges));
}

HistBinning HistBinning::Variable(std::string title,
                                  std::vector<double> binEdges) {
  return HistBinning(std::move(title), std::move(binEdges));
}

HistBinning HistBinning::Logarithmic(std::string title, std::size_t bins,
                                     double bMin, double bMax) {
  std::vector<double> binEdges(bins + 1);
  const double logMin = std::log10(bMin);
  const double logMax = std::log10(bMax);
  const double step = (logMax - logMin) / bins;
  for (std::size_t i = 0; i <= bins; ++i) {
    binEdges[i] = std::pow(10, logMin + i * step);
  }
  return HistBinning(std::move(title), std::move(binEdges));
}

Histogram1D::Histogram1D(std::string name, std::string title,
                         const HistBinning& binning)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_axisTitle(binning.title()),
      m_hist(boost::histogram::make_histogram(
          detail::BoostVariableAxis(binning.binEdges()))) {}

Histogram1D::Histogram1D(std::string name, std::string title,
                         std::string axisTitle, detail::BoostHist1D hist)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_axisTitle(std::move(axisTitle)),
      m_hist(std::move(hist)) {}

void Histogram1D::fill(double value) {
  m_hist(value);
}

Histogram2D::Histogram2D(std::string name, std::string title,
                         const HistBinning& xBinning,
                         const HistBinning& yBinning)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_xAxisTitle(xBinning.title()),
      m_yAxisTitle(yBinning.title()),
      m_hist(boost::histogram::make_histogram(
          detail::BoostVariableAxis(xBinning.binEdges()),
          detail::BoostVariableAxis(yBinning.binEdges()))) {}

void Histogram2D::fill(double xValue, double yValue) {
  m_hist(xValue, yValue);
}

Histogram1D Histogram2D::projectionX() const {
  auto projectedHist = boost::histogram::algorithm::project(
      m_hist, std::integral_constant<unsigned, 0>{});
  return Histogram1D(m_name + "_projX", m_title + " projection X", m_xAxisTitle,
                     std::move(projectedHist));
}

Histogram1D Histogram2D::projectionY() const {
  auto projectedHist = boost::histogram::algorithm::project(
      m_hist, std::integral_constant<unsigned, 1>{});
  return Histogram1D(m_name + "_projY", m_title + " projection Y", m_yAxisTitle,
                     std::move(projectedHist));
}

ProfileHistogram::ProfileHistogram(std::string name, std::string title,
                                   const HistBinning& xBinning,
                                   std::string yAxisTitle)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_xAxisTitle(xBinning.title()),
      m_yAxisTitle(std::move(yAxisTitle)),
      m_hist(boost::histogram::make_weighted_profile(
          detail::BoostVariableAxis(xBinning.binEdges()))) {}

void ProfileHistogram::fill(double xValue, double yValue) {
  // Fill with weight=1.0. This way, the weighted_mean accumulator behaves
  // identically to the regular mean accumulator for unweighted data, but
  // additionally provides sum_of_weights_squared() which is needed for the
  // error calculation in ROOT.
  m_hist(xValue, boost::histogram::weight(1.0),
         boost::histogram::sample(yValue));
}

Efficiency1D::Efficiency1D(std::string name, std::string title,
                           const HistBinning& binning)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_axisTitle(binning.title()),
      m_passed(boost::histogram::make_histogram(
          detail::BoostVariableAxis(binning.binEdges()))),
      m_total(boost::histogram::make_histogram(
          detail::BoostVariableAxis(binning.binEdges()))) {}

void Efficiency1D::fill(double value, bool passed) {
  m_total(value);
  if (passed) {
    m_passed(value);
  }
}

Efficiency2D::Efficiency2D(std::string name, std::string title,
                           const HistBinning& xBinning,
                           const HistBinning& yBinning)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_xAxisTitle(xBinning.title()),
      m_yAxisTitle(yBinning.title()),
      m_passed(boost::histogram::make_histogram(
          detail::BoostVariableAxis(xBinning.binEdges()),
          detail::BoostVariableAxis(yBinning.binEdges()))),
      m_total(boost::histogram::make_histogram(
          detail::BoostVariableAxis(xBinning.binEdges()),
          detail::BoostVariableAxis(yBinning.binEdges()))) {}

void Efficiency2D::fill(double xValue, double yValue, bool passed) {
  m_total(xValue, yValue);
  if (passed) {
    m_passed(xValue, yValue);
  }
}

}  // namespace Acts
