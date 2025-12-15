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

/// @brief 1D histogram wrapper using boost::histogram for data collection
///
/// This class wraps boost::histogram to provide a ROOT-independent histogram
/// implementation. It supports variable binning (uniform, logarithmic, and
/// custom bin edges) and weighted fills with proper error tracking.
///
/// The class provides minimal API - converters access the boost histogram
/// directly via histogram() for iteration and bin content extraction.
class BoostHistogram1D {
 public:
  /// Construct 1D histogram from binning specification
  ///
  /// @param name Histogram name (for identification and output)
  /// @param title Histogram title (for plotting)
  /// @param binning Binning specification (uniform, variable, or logarithmic)
  BoostHistogram1D(std::string name, std::string title,
                   const PlotHelpers::Binning& binning);

  /// Fill histogram with value
  ///
  /// @param value Value to fill
  void fill(double value);

  /// Get histogram name
  const std::string& name() const { return m_name; }

  /// Get histogram title
  const std::string& title() const { return m_title; }

  /// Get axis title
  const std::string& axisTitle() const { return m_axisTitle; }

  /// Direct access to boost::histogram (for converters and tests)
  const auto& histogram() const { return m_hist; }

 private:
  std::string m_name;
  std::string m_title;
  std::string m_axisTitle;

  /// Boost histogram with variable axis and unlimited storage
  ///
  /// - axis::variable<double>: supports non-uniform binning
  /// - unlimited_storage<>: tracks sum-of-weights and sum-of-weights-squared
  ///   (equivalent to TH1F::Sumw2() behavior)
  boost::histogram::histogram<
      std::tuple<boost::histogram::axis::variable<double>>,
      boost::histogram::unlimited_storage<>> m_hist;
};

/// @brief 2D histogram wrapper using boost::histogram for data collection
///
/// This class wraps boost::histogram to provide a ROOT-independent 2D histogram
/// implementation. It supports variable binning on both axes and weighted fills
/// with proper error tracking.
///
/// The class provides minimal API - converters access the boost histogram
/// directly via histogram(). Projections are handled by ROOT converters.
class BoostHistogram2D {
 public:
  /// Construct 2D histogram from binning specifications
  ///
  /// @param name Histogram name (for identification and output)
  /// @param title Histogram title (for plotting)
  /// @param xBinning X-axis binning specification
  /// @param yBinning Y-axis binning specification
  BoostHistogram2D(std::string name, std::string title,
                   const PlotHelpers::Binning& xBinning,
                   const PlotHelpers::Binning& yBinning);

  /// Fill histogram with x, y values
  ///
  /// @param xValue X-axis value to fill
  /// @param yValue Y-axis value to fill
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

  /// Boost histogram with two variable axes and unlimited storage
  ///
  /// - Two axis::variable<double>: supports non-uniform binning on both axes
  /// - unlimited_storage<>: tracks sum-of-weights and sum-of-weights-squared
  boost::histogram::histogram<
      std::tuple<boost::histogram::axis::variable<double>,
                 boost::histogram::axis::variable<double>>,
      boost::histogram::unlimited_storage<>> m_hist;
};

}  // namespace ActsExamples
