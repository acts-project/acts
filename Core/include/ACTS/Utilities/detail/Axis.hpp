// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace Acts {

namespace detail {

  /// @brief calculate bin indices from a given binning structure
  ///
  /// This class provides some basic functionality for calculating bin indices
  /// for a given binning configuration. Both equidistant as well as variable
  /// binning structures are supported.
  ///
  /// Bin intervals are defined such that the lower bound is closed and the
  /// uper bound is open.
  ///
  /// @tparam equidistant flag whether binning is equidistant (@c true)
  ///                     or not (@c false)
  template <bool equidistant>
  class Axis;

  typedef Axis<true>  EquidistantAxis;
  typedef Axis<false> VariableAxis;

  /// @brief calculate bin indices for an equidistant binning
  ///
  /// This class provides some basic functionality for calculating bin indices
  /// for a given equidistant binning.
  template <>
  class Axis<true>
  {
  public:
    /// @brief default constructor
    ///
    /// @param [in] xmin lower boundary of axis range
    /// @param [in] xmax upper boundary of axis range
    /// @param [in] nBins number of bins to divide the axis range into
    ///
    /// Divide the range \f$[xmin,xmax)\f$ into \f$nBins\f$ equidistant bins.
    Axis(double xmin, double xmax, size_t nBins)
      : m_min(xmin), m_max(xmax), m_width((xmax - xmin) / nBins), m_bins(nBins)
    {
    }

    /// @brief get corresponding bin index for given coordinate
    ///
    /// @param  [in] x input coordinate
    /// @return index of bin containing the given value
    ///
    /// @note Bin intervals are defined with closed lower bounds and open upper
    ///       bounds, that is \f$l <= x < u\f$ if the value @c x lies within a
    ///       bin with lower bound @c l and upper bound @c u.
    /// @note Bin indices start at @c 1. The underflow bin has the index @c 0
    ///       while the index <tt>nBins + 1</tt> indicates the overflow bin .
    size_t
    getBin(double x) const
    {
      if (x < getMin()) return 0u;
      if (x >= getMax()) return getNBins() + 1;

      return std::floor((x - getMin()) / getBinWidth()) + 1;
    }

    /// @brief get bin width
    ///
    /// @return constant bin width for all bins
    double
    getBinWidth() const
    {
      return m_width;
    }

    /// @brief get maximum of binning range
    ///
    /// @return maximum of binning range
    double
    getMax() const
    {
      return m_max;
    }

    /// @brief get minimum of binning range
    ///
    /// @return minimum of binning range
    double
    getMin() const
    {
      return m_min;
    }

    /// @brief get total number of bins
    ///
    /// @return total number of bins (excluding under-/overflow bins)
    size_t
    getNBins() const
    {
      return m_bins;
    }

  private:
    /// minimum of binning range
    double m_min;
    /// maximum of binning range
    double m_max;
    /// constant bin width
    double m_width;
    /// number of bins (excluding under-/overflow bins)
    size_t m_bins;
  };

  /// @brief calculate bin indices for a variable binning
  ///
  /// This class provides some basic functionality for calculating bin indices
  /// for a given binning with variable bin sizes.
  template <>
  class Axis<false>
  {
  public:
    /// @brief default constructor
    ///
    /// @param [in] binEdges vector of bin edges
    /// @pre @c binEdges must be strictly sorted in ascending order.
    /// @pre @c binEdges must contain at least two entries.
    ///
    /// Create a binning structure with @c nBins variable-sized bins from the
    /// given bin boundaries. @c nBins is given by the number of bin edges
    /// reduced by one.
    Axis(std::vector<double> binEdges) : m_binEdges(std::move(binEdges)) {}

    /// @brief get corresponding bin index for given coordinate
    ///
    /// @param  [in] x input coordinate
    /// @return index of bin containing the given value
    ///
    /// @note Bin intervals are defined with closed lower bounds and open upper
    ///       bounds, that is \f$l <= x < u\f$ if the value @c x lies within a
    ///       bin with lower bound @c l and upper bound @c u.
    /// @note Bin indices start at @c 1. The underflow bin has the index @c 0
    ///       while the index <tt>nBins + 1</tt> indicates the overflow bin .
    size_t
    getBin(double x) const
    {
      if (x < getMin()) return 0u;
      if (x >= getMax()) return getNBins() + 1;

      const auto it
          = std::upper_bound(std::begin(m_binEdges), std::end(m_binEdges), x);
      return std::distance(std::begin(m_binEdges), it);
    }

    /// @brief get maximum of binning range
    ///
    /// @return maximum of binning range
    double
    getMax() const
    {
      return m_binEdges.back();
    }

    /// @brief get minimum of binning range
    ///
    /// @return minimum of binning range
    double
    getMin() const
    {
      return m_binEdges.front();
    }

    /// @brief get total number of bins
    ///
    /// @return total number of bins (excluding under-/overflow bins)
    size_t
    getNBins() const
    {
      return m_binEdges.size() - 1;
    }

  private:
    /// vector of bin edges (sorted in ascending order)
    std::vector<double> m_binEdges;
  };
}  // namespace detail

}  // namespace Acts
