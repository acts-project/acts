// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Histogram.hpp"

namespace Acts {

Histogram1D::Histogram1D(std::string name, std::string title,
                         BoostVariableAxis axis)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_hist(boost::histogram::make_histogram(std::move(axis))) {}

Histogram1D::Histogram1D(std::string name, std::string title, BoostHist1D hist)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_hist(std::move(hist)) {}

void Histogram1D::fill(double value) {
  m_hist(value);
}

Histogram2D::Histogram2D(std::string name, std::string title,
                         BoostVariableAxis xAxis, BoostVariableAxis yAxis)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_hist(boost::histogram::make_histogram(std::move(xAxis),
                                              std::move(yAxis))) {}

void Histogram2D::fill(double xValue, double yValue) {
  m_hist(xValue, yValue);
}

Histogram1D Histogram2D::projectionX() const {
  auto projectedHist = boost::histogram::algorithm::project(
      m_hist, std::integral_constant<unsigned, 0>{});
  return Histogram1D(m_name + "_projX", m_title + " projection X",
                     std::move(projectedHist));
}

Histogram1D Histogram2D::projectionY() const {
  auto projectedHist = boost::histogram::algorithm::project(
      m_hist, std::integral_constant<unsigned, 1>{});
  return Histogram1D(m_name + "_projY", m_title + " projection Y",
                     std::move(projectedHist));
}

ProfileHistogram::ProfileHistogram(std::string name, std::string title,
                                   BoostVariableAxis xAxis,
                                   std::string yAxisTitle)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_yAxisTitle(std::move(yAxisTitle)),
      m_hist(boost::histogram::make_profile(std::move(xAxis))) {}

void ProfileHistogram::fill(double xValue, double yValue) {
  m_hist(xValue, boost::histogram::sample(yValue));
}

Efficiency1D::Efficiency1D(std::string name, std::string title,
                           BoostVariableAxis axis)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_accepted(boost::histogram::make_histogram(axis)),
      m_total(boost::histogram::make_histogram(std::move(axis))) {}

void Efficiency1D::fill(double value, bool accepted) {
  m_total(value);
  if (accepted) {
    m_accepted(value);
  }
}

Efficiency2D::Efficiency2D(std::string name, std::string title,
                           BoostVariableAxis xAxis, BoostVariableAxis yAxis)
    : m_name(std::move(name)),
      m_title(std::move(title)),
      m_accepted(boost::histogram::make_histogram(xAxis, yAxis)),
      m_total(boost::histogram::make_histogram(std::move(xAxis),
                                               std::move(yAxis))) {}

void Efficiency2D::fill(double xValue, double yValue, bool accepted) {
  m_total(xValue, yValue);
  if (accepted) {
    m_accepted(xValue, yValue);
  }
}

}  // namespace Acts
