// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <string>
#include <vector>

#include <boost/histogram.hpp>

namespace Acts {
namespace detail {
using BoostVariableAxis = boost::histogram::axis::variable<double>;
using BoostHist1D = decltype(boost::histogram::make_histogram(
    std::declval<BoostVariableAxis>()));
using BoostHist2D = decltype(boost::histogram::make_histogram(
    std::declval<BoostVariableAxis>(), std::declval<BoostVariableAxis>()));
using BoostProfileHist =
    decltype(boost::histogram::make_profile(std::declval<BoostVariableAxis>()));
}  // namespace detail

/// @brief Nested binning struct for booking plots
class HistBinning {
 public:
  static HistBinning Uniform(std::string title, std::size_t bins, double bMin,
                             double bMax);
  static HistBinning Variable(std::string title, std::vector<double> binEdges);
  static HistBinning Logarithmic(std::string title, std::size_t bins,
                                 double bMin, double bMax);

  HistBinning(std::string title, std::vector<double> binEdges)
      : m_title(std::move(title)), m_binEdges(std::move(binEdges)) {}

  const std::string& title() const { return m_title; }
  std::size_t nBins() const { return m_binEdges.size() - 1; }
  const std::vector<double>& binEdges() const { return m_binEdges; }
  double low() const { return m_binEdges.front(); }
  double high() const { return m_binEdges.back(); }

 private:
  std::string m_title;
  std::vector<double> m_binEdges;
};

/// @brief 1D histogram wrapper using boost::histogram for data collection
///
/// This class wraps boost::histogram to provide a ROOT-independent histogram
/// implementation. It supports variable binning (uniform, logarithmic, and
/// custom bin edges) and weighted fills with proper error tracking.
///
/// The class provides minimal API - converters access the boost histogram
/// directly via histogram() for iteration and bin content extraction.
class Histogram1D {
 public:
  /// Construct 1D histogram from binning specification
  ///
  /// @param name Histogram name (for identification and output)
  /// @param title Histogram title (for plotting)
  /// @param binning Binning specification (uniform, variable, or logarithmic)
  Histogram1D(std::string name, std::string title, const HistBinning& binning);

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

  detail::BoostHist1D m_hist;
};

/// @brief 2D histogram wrapper using boost::histogram for data collection
///
/// This class wraps boost::histogram to provide a ROOT-independent 2D histogram
/// implementation. It supports variable binning on both axes and weighted fills
/// with proper error tracking.
///
/// The class provides minimal API - converters access the boost histogram
/// directly via histogram(). Projections are handled by ROOT converters.
class Histogram2D {
 public:
  /// Construct 2D histogram from binning specifications
  ///
  /// @param name Histogram name (for identification and output)
  /// @param title Histogram title (for plotting)
  /// @param xBinning X-axis binning specification
  /// @param yBinning Y-axis binning specification
  Histogram2D(std::string name, std::string title, const HistBinning& xBinning,
              const HistBinning& yBinning);

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

  detail::BoostHist2D m_hist;
};

/// @brief Profile histogram using boost::histogram
///
/// This class wraps boost::histogram to provide a ROOT-independent profile
/// histogram implementation. For each X bin, it tracks the mean (and variance)
/// of Y values.
class ProfileHistogram {
 public:
  /// Construct profile histogram from X binning specification
  ///
  /// @param name Histogram name
  /// @param title Histogram title
  /// @param xBinning X-axis binning specification
  /// @param yAxisTitle Y-axis title
  ProfileHistogram(std::string name, std::string title,
                   const HistBinning& xBinning, std::string yAxisTitle);

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

  detail::BoostProfileHist m_hist;
};

/// @brief 1D efficiency histogram using boost::histogram
///
/// This class tracks pass/total counts for efficiency calculation.
/// It internally uses two simple histograms: one for passed events,
/// one for total events.
class Efficiency1D {
 public:
  /// Construct 1D efficiency histogram
  ///
  /// @param name Histogram name
  /// @param title Histogram title
  /// @param binning Binning specification
  Efficiency1D(std::string name, std::string title, const HistBinning& binning);

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

  detail::BoostHist1D m_passed, m_total;
};

/// @brief 2D efficiency histogram using boost::histogram
///
/// This class tracks pass/total counts for 2D efficiency calculation.
class Efficiency2D {
 public:
  /// Construct 2D efficiency histogram
  ///
  /// @param name Histogram name
  /// @param title Histogram title
  /// @param xBinning X-axis binning specification
  /// @param yBinning Y-axis binning specification
  Efficiency2D(std::string name, std::string title, const HistBinning& xBinning,
               const HistBinning& yBinning);

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

  detail::BoostHist2D m_passed, m_total;
};

}  // namespace Acts
