// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Utilities/Helpers.hpp"

#include <boost/histogram.hpp>
#include <string>

namespace ActsExamples {

/// @brief Profile histogram using boost::histogram with weighted_mean accumulator
///
/// This class wraps boost::histogram to provide a ROOT-independent profile
/// histogram implementation. For each X bin, it tracks the mean (and variance)
/// of Y values.
class BoostProfileHistogram {
 public:
  /// Construct profile histogram from X binning specification
  ///
  /// @param name Histogram name
  /// @param title Histogram title
  /// @param xBinning X-axis binning specification
  /// @param yAxisTitle Y-axis title
  BoostProfileHistogram(std::string name, std::string title,
                        const PlotHelpers::Binning& xBinning,
                        std::string yAxisTitle);

  /// Fill profile with (x, y) pair
  ///
  /// @param xValue X value (bin coordinate)
  /// @param yValue Y value (profiled quantity)
  void fill(double xValue, double yValue);

  /// Get histogram name
  const std::string& name() const { return m_name; }

  /// Get histogram title
  const std::string& title() const { return m_title; }

  /// Get X-axis title
  const std::string& xAxisTitle() const { return m_xAxisTitle; }

  /// Get Y-axis title
  const std::string& yAxisTitle() const { return m_yAxisTitle; }

  /// Direct access to boost::histogram (for converters and tests)
  const auto& histogram() const { return m_hist; }

 private:
  std::string m_name;
  std::string m_title;
  std::string m_xAxisTitle;
  std::string m_yAxisTitle;

  /// Boost profile histogram with weighted_mean accumulator
  boost::histogram::histogram<
      std::tuple<boost::histogram::axis::variable<double>>,
      boost::histogram::dense_storage<
          boost::histogram::accumulators::weighted_mean<>>>
      m_hist;
};

/// @brief 1D efficiency histogram using boost::histogram
///
/// This class tracks pass/total counts for efficiency calculation.
/// It internally uses two simple histograms: one for passed events,
/// one for total events.
class BoostEfficiency1D {
 public:
  /// Construct 1D efficiency histogram
  ///
  /// @param name Histogram name
  /// @param title Histogram title
  /// @param binning Binning specification
  BoostEfficiency1D(std::string name, std::string title,
                    const PlotHelpers::Binning& binning);

  /// Fill efficiency histogram
  ///
  /// @param value Value to fill
  /// @param passed Whether the event passed selection
  void fill(double value, bool passed);

  /// Get histogram name
  const std::string& name() const { return m_name; }

  /// Get histogram title
  const std::string& title() const { return m_title; }

  /// Get axis title
  const std::string& axisTitle() const { return m_axisTitle; }

  /// Access to passed histogram (for converters and tests)
  const auto& passedHistogram() const { return m_passed; }

  /// Access to total histogram (for converters and tests)
  const auto& totalHistogram() const { return m_total; }

 private:
  std::string m_name;
  std::string m_title;
  std::string m_axisTitle;

  /// Histogram for passed events
  boost::histogram::histogram<
      std::tuple<boost::histogram::axis::variable<double>>,
      boost::histogram::unlimited_storage<>>
      m_passed;

  /// Histogram for total events
  boost::histogram::histogram<
      std::tuple<boost::histogram::axis::variable<double>>,
      boost::histogram::unlimited_storage<>>
      m_total;
};

/// @brief 2D efficiency histogram using boost::histogram
///
/// This class tracks pass/total counts for 2D efficiency calculation.
class BoostEfficiency2D {
 public:
  /// Construct 2D efficiency histogram
  ///
  /// @param name Histogram name
  /// @param title Histogram title
  /// @param xBinning X-axis binning specification
  /// @param yBinning Y-axis binning specification
  BoostEfficiency2D(std::string name, std::string title,
                    const PlotHelpers::Binning& xBinning,
                    const PlotHelpers::Binning& yBinning);

  /// Fill efficiency histogram
  ///
  /// @param xValue X value
  /// @param yValue Y value
  /// @param passed Whether the event passed selection
  void fill(double xValue, double yValue, bool passed);

  /// Get histogram name
  const std::string& name() const { return m_name; }

  /// Get histogram title
  const std::string& title() const { return m_title; }

  /// Get X-axis title
  const std::string& xAxisTitle() const { return m_xAxisTitle; }

  /// Get Y-axis title
  const std::string& yAxisTitle() const { return m_yAxisTitle; }

  /// Access to passed histogram (for converters and tests)
  const auto& passedHistogram() const { return m_passed; }

  /// Access to total histogram (for converters and tests)
  const auto& totalHistogram() const { return m_total; }

 private:
  std::string m_name;
  std::string m_title;
  std::string m_xAxisTitle;
  std::string m_yAxisTitle;

  /// Histogram for passed events
  boost::histogram::histogram<
      std::tuple<boost::histogram::axis::variable<double>,
                 boost::histogram::axis::variable<double>>,
      boost::histogram::unlimited_storage<>>
      m_passed;

  /// Histogram for total events
  boost::histogram::histogram<
      std::tuple<boost::histogram::axis::variable<double>,
                 boost::histogram::axis::variable<double>>,
      boost::histogram::unlimited_storage<>>
      m_total;
};

}  // namespace ActsExamples
