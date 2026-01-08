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

/// @brief Boost axis type to use for histograms with metadata support
using BoostVariableAxis = boost::histogram::axis::variable<double, std::string>;

/// @brief Underlying Boost type for Histogram1D
using BoostHist1D = decltype(boost::histogram::make_histogram(
    std::declval<BoostVariableAxis>()));

/// @brief Underlying Boost type for Histogram2D
using BoostHist2D = decltype(boost::histogram::make_histogram(
    std::declval<BoostVariableAxis>(), std::declval<BoostVariableAxis>()));

/// @brief Underlying Boost type for ProfileHistogram
using BoostProfileHist =
    decltype(boost::histogram::make_profile(std::declval<BoostVariableAxis>()));

/// @brief 1D histogram wrapper using boost::histogram for data collection
///
/// This class wraps boost::histogram to provide a ROOT-independent histogram
/// implementation.
class Histogram1D {
 public:
  /// Construct 1D histogram from axis
  ///
  /// @param name Histogram name (for identification and output)
  /// @param title Histogram title (for plotting)
  /// @param axis Axis with binning and metadata
  Histogram1D(std::string name, std::string title, BoostVariableAxis axis);

  /// Fill histogram with value
  ///
  /// @param value Value to fill
  void fill(double value);

  /// Get histogram name
  const std::string& name() const { return m_name; }

  /// Get histogram title
  const std::string& title() const { return m_title; }

  /// Get axis title from axis metadata
  const std::string& axisTitle() const { return m_hist.axis(0).metadata(); }

  /// Direct access to boost::histogram (for converters and tests)
  const BoostHist1D& histogram() const { return m_hist; }

 private:
  friend class Histogram2D;

  /// Construct 1D histogram from existing boost histogram
  ///
  /// @param name Histogram name (for identification and output)
  /// @param title Histogram title (for plotting)
  /// @param hist Boost histogram to wrap
  Histogram1D(std::string name, std::string title, BoostHist1D hist);

  std::string m_name;
  std::string m_title;

  BoostHist1D m_hist;
};

/// @brief 2D histogram wrapper using boost::histogram for data collection
///
/// This class wraps boost::histogram to provide a ROOT-independent 2D histogram
/// implementation.
class Histogram2D {
 public:
  /// Construct 2D histogram from axes
  ///
  /// @param name Histogram name (for identification and output)
  /// @param title Histogram title (for plotting)
  /// @param xAxis X-axis with binning and metadata
  /// @param yAxis Y-axis with binning and metadata
  Histogram2D(std::string name, std::string title, BoostVariableAxis xAxis,
              BoostVariableAxis yAxis);

  /// Fill histogram with x, y values
  ///
  /// @param xValue X-axis value to fill
  /// @param yValue Y-axis value to fill
  void fill(double xValue, double yValue);

  /// Get histogram name
  const std::string& name() const { return m_name; }

  /// Get histogram title
  const std::string& title() const { return m_title; }

  /// Get X-axis title from axis metadata
  const std::string& xAxisTitle() const { return m_hist.axis(0).metadata(); }

  /// Get Y-axis title from axis metadata
  const std::string& yAxisTitle() const { return m_hist.axis(1).metadata(); }

  /// Project the histogram onto x
  Histogram1D projectionX() const;

  /// Project the histogram onto y
  Histogram1D projectionY() const;

  /// Direct access to boost::histogram (for converters and tests)
  const BoostHist2D& histogram() const { return m_hist; }

 private:
  std::string m_name;
  std::string m_title;

  BoostHist2D m_hist;
};

/// @brief Profile histogram using boost::histogram
///
/// This class wraps boost::histogram to provide a ROOT-independent profile
/// histogram implementation. For each X bin, it tracks the mean and variance
/// of Y values.
class ProfileHistogram {
 public:
  /// Construct profile histogram from X axis
  ///
  /// @param name Histogram name
  /// @param title Histogram title
  /// @param xAxis X-axis with binning and metadata
  /// @param yAxisTitle Y-axis title
  ProfileHistogram(std::string name, std::string title, BoostVariableAxis xAxis,
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

  /// Get X-axis title from axis metadata
  const std::string& xAxisTitle() const { return m_hist.axis(0).metadata(); }

  /// Get Y-axis title
  const std::string& yAxisTitle() const { return m_yAxisTitle; }

  /// Direct access to boost::histogram (for converters and tests)
  const BoostProfileHist& histogram() const { return m_hist; }

 private:
  std::string m_name;
  std::string m_title;
  std::string m_yAxisTitle;

  BoostProfileHist m_hist;
};

/// @brief 1D efficiency histogram using boost::histogram
///
/// This class tracks pass/total counts for efficiency calculation.
/// It internally uses two 1D histograms: one for accepted events,
/// one for total events.
class Efficiency1D {
 public:
  /// Construct 1D efficiency histogram
  ///
  /// @param name Histogram name
  /// @param title Histogram title
  /// @param axis Axis with binning and metadata
  Efficiency1D(std::string name, std::string title, BoostVariableAxis axis);

  /// Fill efficiency histogram
  ///
  /// @param value Value to fill
  /// @param accepted Whether the event accepted selection
  void fill(double value, bool accepted);

  /// Get histogram name
  const std::string& name() const { return m_name; }

  /// Get histogram title
  const std::string& title() const { return m_title; }

  /// Get axis title from axis metadata
  const std::string& axisTitle() const { return m_accepted.axis(0).metadata(); }

  /// Access to accepted histogram (for converters and tests)
  const BoostHist1D& acceptedHistogram() const { return m_accepted; }

  /// Access to total histogram (for converters and tests)
  const BoostHist1D& totalHistogram() const { return m_total; }

 private:
  std::string m_name;
  std::string m_title;

  BoostHist1D m_accepted;
  BoostHist1D m_total;
};

/// @brief 2D efficiency histogram using boost::histogram
///
/// This class tracks pass/total counts for 2D efficiency calculation.
/// It internally uses two 1D histograms: one for accepted events,
/// one for total events.
class Efficiency2D {
 public:
  /// Construct 2D efficiency histogram
  ///
  /// @param name Histogram name
  /// @param title Histogram title
  /// @param xAxis X-axis with binning and metadata
  /// @param yAxis Y-axis with binning and metadata
  Efficiency2D(std::string name, std::string title, BoostVariableAxis xAxis,
               BoostVariableAxis yAxis);

  /// Fill efficiency histogram
  ///
  /// @param xValue X value
  /// @param yValue Y value
  /// @param accepted Whether the event accepted selection
  void fill(double xValue, double yValue, bool accepted);

  /// Get histogram name
  const std::string& name() const { return m_name; }

  /// Get histogram title
  const std::string& title() const { return m_title; }

  /// Get X-axis title from axis metadata
  const std::string& xAxisTitle() const {
    return m_accepted.axis(0).metadata();
  }

  /// Get Y-axis title from axis metadata
  const std::string& yAxisTitle() const {
    return m_accepted.axis(1).metadata();
  }

  /// Access to accepted histogram (for converters and tests)
  const BoostHist2D& acceptedHistogram() const { return m_accepted; }

  /// Access to total histogram (for converters and tests)
  const BoostHist2D& totalHistogram() const { return m_total; }

 private:
  std::string m_name;
  std::string m_title;

  BoostHist2D m_accepted;
  BoostHist2D m_total;
};

}  // namespace Acts
