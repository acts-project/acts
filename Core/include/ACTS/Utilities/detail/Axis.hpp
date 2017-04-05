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
    /// Divide the range \f$[\text{xmin},\text{xmax})\f$ into \f$\text{nBins}\f$
    /// equidistant bins.
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
    /// @return constant width for all bins
    double getBinWidth(size_t = 0) const { return m_width; }

    /// @brief get lower bound of bin
    ///
    /// @param  [in] bin index of bin
    /// @return lower bin boundary
    ///
    /// @pre @c bin must be a valid bin index (excluding under-/overflow bins),
    ///      i.e. \f$1 \le \text{bin} \le \text{nBins}\f$
    ///
    /// @note Bin intervals have a closed lower bound, i.e. the lower boundary
    ///       belongs to the bin with the given bin index.
    double
    getBinLowerBound(size_t bin) const
    {
      return getMin() + (bin - 1) * getBinWidth();
    }

    /// @brief get upper bound of bin
    ///
    /// @param  [in] bin index of bin
    /// @return upper bin boundary
    ///
    /// @pre @c bin must be a valid bin index (excluding under-/overflow bins),
    ///      i.e. \f$1 \le \text{bin} \le \text{nBins}\f$
    ///
    /// @note Bin intervals have an open upper bound, i.e. the upper boundary
    ///       does @b not belong to the bin with the given bin index.
    double
    getBinUpperBound(size_t bin) const
    {
      return getMin() + bin * getBinWidth();
    }

    /// @brief get bin center
    ///
    /// @param  [in] bin index of bin
    /// @return bin center position
    ///
    /// @pre @c bin must be a valid bin index (excluding under-/overflow bins),
    ///      i.e. \f$1 \le \text{bin} \le \text{nBins}\f$
    double
    getBinCenter(size_t bin) const
    {
      return getMin() + (bin - 0.5) * getBinWidth();
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

    /// @brief get bin width
    ///
    /// @param  [in] bin index of bin
    /// @return width of given bin
    ///
    /// @pre @c bin must be a valid bin index (excluding under-/overflow bins),
    ///      i.e. \f$1 \le \text{bin} \le \text{nBins}\f$
    double
    getBinWidth(size_t bin) const
    {
      return m_binEdges.at(bin) - m_binEdges.at(bin - 1);
    }

    /// @brief get lower bound of bin
    ///
    /// @param  [in] bin index of bin
    /// @return lower bin boundary
    ///
    /// @pre @c bin must be a valid bin index (excluding under-/overflow bins),
    ///      i.e. \f$1 \le \text{bin} \le \text{nBins}\f$
    ///
    /// @note Bin intervals have a closed lower bound, i.e. the lower boundary
    ///       belongs to the bin with the given bin index.
    double
    getBinLowerBound(size_t bin) const
    {
      return m_binEdges.at(bin - 1);
    }

    /// @brief get upper bound of bin
    ///
    /// @param  [in] bin index of bin
    /// @return upper bin boundary
    ///
    /// @pre @c bin must be a valid bin index (excluding under-/overflow bins),
    ///      i.e. \f$1 \le \text{bin} \le \text{nBins}\f$
    ///
    /// @note Bin intervals have an open upper bound, i.e. the upper boundary
    ///       does @b not belong to the bin with the given bin index.
    double
    getBinUpperBound(size_t bin) const
    {
      return m_binEdges.at(bin);
    }

    /// @brief get bin center
    ///
    /// @param  [in] bin index of bin
    /// @return bin center position
    ///
    /// @pre @c bin must be a valid bin index (excluding under-/overflow bins),
    ///      i.e. \f$1 \le \text{bin} \le \text{nBins}\f$
    double
    getBinCenter(size_t bin) const
    {
      return 0.5 * (getBinLowerBound(bin) + getBinUpperBound(bin));
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
