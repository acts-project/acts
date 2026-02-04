// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/RangeXD.hpp"

#include <array>
#include <string>
#include <tuple>

#include <boost/histogram.hpp>

namespace Acts::Experimental {

/// Variable-width histogram axis with string metadata
using BoostVariableAxis = boost::histogram::axis::variable<double, std::string>;
/// Regular-width histogram axis with string metadata
using BoostRegularAxis =
    boost::histogram::axis::regular<double, boost::histogram::use_default,
                                    std::string>;
/// Logarithmic-scale histogram axis with string metadata
using BoostLogAxis = boost::histogram::axis::regular<
    double, boost::histogram::axis::transform::log, std::string>;

/// @brief Boost axis variant supporting variable, regular, and log-scale axes with metadata
/// @note It seems not to be possible to combine compile-time fixed number of
///       axes with `boost::histogram::axis::variant`. Therefore we use
///       `std::vector<AxisVariant>` internally.
using AxisVariant =
    boost::histogram::axis::variant<BoostVariableAxis, BoostRegularAxis,
                                    BoostLogAxis>;

/// @brief Underlying Boost type for histograms
using BoostHist = decltype(boost::histogram::make_histogram(
    std::declval<std::vector<AxisVariant>>()));

/// @brief Underlying Boost type for ProfileHistogram
using BoostProfileHist = decltype(boost::histogram::make_profile(
    std::declval<std::vector<AxisVariant>>()));

/// @brief Multi-dimensional histogram wrapper using boost::histogram for data collection
///
/// This class wraps boost::histogram to provide a ROOT-independent histogram
/// implementation with compile-time dimensionality.
///
/// @tparam Dim Number of dimensions
template <std::size_t Dim>
class Histogram {
 public:
  /// Construct multi-dimensional histogram from axes
  ///
  /// @param name Histogram name (for identification and output)
  /// @param title Histogram title (for plotting)
  /// @param axes Array of axes with binning and metadata
  Histogram(std::string name, std::string title,
            const std::array<AxisVariant, Dim>& axes)
      : m_name(std::move(name)),
        m_title(std::move(title)),
        m_hist(boost::histogram::make_histogram(axes.begin(), axes.end())) {}

  /// Copy constructor
  /// @param other The other histogram to copy from
  Histogram(const Histogram& other) = default;

  /// Move constructor
  /// @param other The other histogram to move from
  Histogram(Histogram&& other) noexcept = default;

  /// Copy assignment operator
  /// @param other The other histogram to copy from
  /// @return The copied histogram
  Histogram& operator=(const Histogram& other) = default;

  /// Move assignment operator
  /// @param other The other histogram to move from
  /// @return The moved histogram
  Histogram& operator=(Histogram&& other) noexcept = default;

  /// Fill histogram with values
  ///
  /// @param values Values to fill (one per axis)
  void fill(const std::array<double, Dim>& values) {
    std::apply([this](auto... v) { m_hist(v...); }, std::tuple_cat(values));
  }

  /// Get histogram name
  /// @return The histogram name
  const std::string& name() const { return m_name; }

  /// Get histogram title
  /// @return The histogram title
  const std::string& title() const { return m_title; }

  /// Get number of dimensions (compile-time constant)
  /// @return The number of dimensions
  static constexpr std::size_t rank() { return Dim; }

  /// Direct access to boost::histogram (for converters and tests)
  /// @return The underlying boost histogram
  const BoostHist& histogram() const { return m_hist; }

 private:
  std::string m_name;
  std::string m_title;

  BoostHist m_hist;
};

/// Type aliases for common dimensions
using Histogram1 = Histogram<1>;
/// 2D histogram
using Histogram2 = Histogram<2>;

/// @brief Multi-dimensional profile histogram using boost::histogram
///
/// This class wraps boost::histogram to provide a ROOT-independent profile
/// histogram implementation with compile-time dimensionality. For each bin,
/// it tracks the mean and variance of sample values.
///
/// @tparam Dim Number of dimensions
template <std::size_t Dim>
class ProfileHistogram {
 public:
  /// Construct multi-dimensional profile histogram from axes
  ///
  /// @param name Histogram name (for identification and output)
  /// @param title Histogram title (for plotting)
  /// @param axes Array of axes with binning and metadata
  /// @param sampleAxisTitle Title for the sampled axis (profiled quantity)
  /// @param sampleRange Samples are discarded when outside range
  ProfileHistogram(std::string name, std::string title,
                   const std::array<AxisVariant, Dim>& axes,
                   std::string sampleAxisTitle,
                   Range1D<double> sampleRange = {})
      : m_name(std::move(name)),
        m_title(std::move(title)),
        m_sampleAxisTitle(std::move(sampleAxisTitle)),
        m_sampleRange(sampleRange),
        m_hist(boost::histogram::make_profile(axes.begin(), axes.end())) {}

  /// Fill profile with values and sample
  ///
  /// @param values Bin coordinate values (one per axis)
  /// @param sample Sample value (profiled quantity)
  void fill(const std::array<double, Dim>& values, double sample) {
    if (!m_sampleRange.contains(sample)) {
      return;
    }

    std::apply(
        [&](auto... v) { m_hist(v..., boost::histogram::sample(sample)); },
        std::tuple_cat(values));
  }

  /// Get histogram name
  /// @return The histogram name
  const std::string& name() const { return m_name; }

  /// Get histogram title
  /// @return The histogram title
  const std::string& title() const { return m_title; }

  /// Get number of dimensions (compile-time constant)
  /// @return The number of dimensions
  static constexpr std::size_t rank() { return Dim; }

  /// Get title of the sample axis
  /// @return The sample axis title
  const std::string& sampleAxisTitle() const { return m_sampleAxisTitle; }

  /// Direct access to boost::histogram (for converters and tests)
  /// @return The underlying boost profile histogram
  const BoostProfileHist& histogram() const { return m_hist; }

 private:
  std::string m_name;
  std::string m_title;
  std::string m_sampleAxisTitle;
  Range1D<double> m_sampleRange;

  BoostProfileHist m_hist;
};

/// Type aliases for common dimensions
using ProfileHistogram1 = ProfileHistogram<1>;

/// @brief Multi-dimensional efficiency histogram using boost::histogram
///
/// This class tracks pass/total counts for efficiency calculation.
/// It internally uses two multi-dimensional histograms: one for accepted
/// events, one for total events.
///
/// @tparam Dim Number of dimensions
template <std::size_t Dim>
class Efficiency {
 public:
  /// Construct multi-dimensional efficiency histogram
  ///
  /// @param name Histogram name
  /// @param title Histogram title
  /// @param axes Array of axes with binning and metadata
  Efficiency(std::string name, std::string title,
             const std::array<AxisVariant, Dim>& axes)
      : m_name(std::move(name)),
        m_title(std::move(title)),
        m_accepted(boost::histogram::make_histogram(axes.begin(), axes.end())),
        m_total(boost::histogram::make_histogram(axes.begin(), axes.end())) {}

  /// Fill efficiency histogram
  ///
  /// @param values Values to fill (one per axis)
  /// @param accepted Whether the event passed selection
  void fill(const std::array<double, Dim>& values, bool accepted) {
    std::apply(
        [&](auto... v) {
          m_total(v...);
          if (accepted) {
            m_accepted(v...);
          }
        },
        std::tuple_cat(values));
  }

  /// Get histogram name
  /// @return The histogram name
  const std::string& name() const { return m_name; }

  /// Get histogram title
  /// @return The histogram title
  const std::string& title() const { return m_title; }

  /// Get number of dimensions (compile-time constant)
  /// @return The histogram dimension
  static constexpr std::size_t rank() { return Dim; }

  /// Access to accepted histogram (for converters and tests)
  /// @return The accepted histogram
  const BoostHist& acceptedHistogram() const { return m_accepted; }

  /// Access to total histogram (for converters and tests)
  /// @return The total histogram
  const BoostHist& totalHistogram() const { return m_total; }

 private:
  std::string m_name;
  std::string m_title;

  BoostHist m_accepted;
  BoostHist m_total;
};

/// Type aliases for common dimensions
using Efficiency1 = Efficiency<1>;
/// 2D efficiency histogram
using Efficiency2 = Efficiency<2>;

/// Project a 2D histogram onto the X axis (axis 0)
///
/// @param hist2d The 2D histogram to project
/// @return A 1D histogram containing the projection
Histogram1 projectionX(const Histogram2& hist2d);

/// Project a 2D histogram onto the Y axis (axis 1)
///
/// @param hist2d The 2D histogram to project
/// @return A 1D histogram containing the projection
Histogram1 projectionY(const Histogram2& hist2d);

/// Extract bin edges from an AxisVariant
///
/// Works with all axis types (regular, variable, log) in the variant by
/// accessing the generic axis interface.
///
/// @param axis The axis variant to extract edges from
/// @return Vector of bin edges (size = nBins + 1)
std::vector<double> extractBinEdges(const AxisVariant& axis);

}  // namespace Acts::Experimental
