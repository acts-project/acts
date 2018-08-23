// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>
#include "Acts/Utilities/IAxis.hpp"

namespace Acts {

namespace detail {

  /// Enum which determines how the axis handle its outer boundaries
  /// possible values values
  /// - Open is the default behaviour: out of bounds
  /// positions are filled into the over or underflow bins
  /// - Bound: out-of-bounds positions resolve to first/last bin
  /// respectively
  /// - Closed: out-of-bounds positions resolve to the outermost
  /// bin on the oppsite side
  enum class AxisBoundaryType { Open, Bound, Closed };

  /// Enum which determines the binning type of the axis
  enum class AxisType { Equidistant, Variable };

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
  template <AxisType type, AxisBoundaryType bdt = AxisBoundaryType::Open>
  class Axis;

  using EquidistantAxis = Axis<AxisType::Equidistant>;
  using VariableAxis    = Axis<AxisType::Variable>;

  /// @brief calculate bin indices for an equidistant binning
  ///
  /// This class provides some basic functionality for calculating bin indices
  /// for a given equidistant binning.
  template <AxisBoundaryType bdt>
  class Axis<AxisType::Equidistant, bdt> final : public IAxis
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

    /// @brief returns whether the axis is equidistant
    ///
    /// @return bool is equidistant
    bool
    isEquidistant() const override
    {
      return true;
    }

    /// @brief returns whether the axis is variable
    ///
    /// @return bool is variable
    bool
    isVariable() const override
    {
      return false;
    }

    /// @brief returns the boundary type set in the template param
    ///
    /// @return @c AxisBoundaryType of this axis
    AxisBoundaryType
    getBoundaryType() const override
    {
      return bdt;
    }

    /// @brief Get #size bins which neighbor the one given
    ///
    /// Generic overload with symmetric size
    ///
    /// @param [in] idx requested bin index
    /// @param [in] sizes how many neighboring bins (up/down)
    /// @return std::set of neighboring bin indices (global)
    std::set<size_t>
    neighborHoodIndices(size_t idx, size_t size = 1) const
    {
      return neighborHoodIndices(idx, std::make_pair(size, size));
    }

    /// @brief Get #size bins which neighbor the one given
    ///
    /// This is the version for Open
    ///
    /// @param [in] idx requested bin index
    /// @param [in] sizes how many neighboring bins (up/down)
    /// @return std::set of neighboring bin indices (global)
    /// @note Open varies given bin and allows 0 and NBins+1 (underflow,
    /// overflow)
    ///       as neighbors
    template <AxisBoundaryType T = bdt,
              std::enable_if_t<T == AxisBoundaryType::Open, int> = 0>
    std::set<size_t>
    neighborHoodIndices(size_t idx,
                        std::pair<size_t, size_t> sizes = {1, 1}) const
    {
      std::set<size_t> result;
      constexpr int    min   = 0;
      const int        max   = getNBins() + 1;
      const int        itmin = std::max(min, int(idx - sizes.first));
      const int        itmax = std::min(max, int(idx + sizes.second));
      for (int i = itmin; i <= itmax; i++) {
        result.insert(i);
      }
      return result;
    }

    /// @brief Get #size bins which neighbor the one given
    ///
    /// This is the version for Bound
    ///
    /// @param [in] idx requested bin index
    /// @param [in] sizes how many neighboring bins (up/down)
    /// @return std::set of neighboring bin indices (global)
    /// @note Bound varies given bin and allows 1 and NBins (regular bins)
    ///       as neighbors
    template <AxisBoundaryType T = bdt,
              std::enable_if_t<T == AxisBoundaryType::Bound, int> = 0>
    std::set<size_t>
    neighborHoodIndices(size_t idx,
                        std::pair<size_t, size_t> sizes = {1, 1}) const
    {
      std::set<size_t> result;
      if (idx <= 0 || idx >= (getNBins() + 1)) {
        return result;
      }
      constexpr int min   = 1;
      const int     max   = getNBins();
      const int     itmin = std::max(min, int(idx - sizes.first));
      const int     itmax = std::min(max, int(idx + sizes.second));
      for (int i = itmin; i <= itmax; i++) {
        result.insert(i);
      }
      return result;
    }

    /// @brief Get #size bins which neighbor the one given
    ///
    /// This is the version for Closed
    ///
    /// @param [in] idx requested bin index
    /// @param [in] sizes how many neighboring bins (up/down)
    /// @return std::set of neighboring bin indices (global)
    /// @note Closed varies given bin and allows bins on the opposite
    ///       side of the axis as neighbors. (excludes underflow / overflow)
    template <AxisBoundaryType T = bdt,
              std::enable_if_t<T == AxisBoundaryType::Closed, int> = 0>
    std::set<size_t>
    neighborHoodIndices(size_t idx,
                        std::pair<size_t, size_t> sizes = {1, 1}) const
    {
      std::set<size_t> result;
      if (idx <= 0 || idx >= (getNBins() + 1)) {
        return result;
      }
      const int itmin = idx - sizes.first;
      const int itmax = idx + sizes.second;
      for (int i = itmin; i <= itmax; i++) {
        result.insert(wrapBin(i));
      }
      return result;
    }

    /// @brief Converts bin index into a valid one for this axis.
    ///
    /// @note Open: bin index is clamped to [0, nBins+1]
    ///
    /// @param [in] bin The bin to wrap
    /// @return valid bin index
    template <AxisBoundaryType T = bdt,
              std::enable_if_t<T == AxisBoundaryType::Open, int> = 0>
    size_t
    wrapBin(int bin) const
    {
      return std::max(std::min(bin, static_cast<int>(getNBins()) + 1), 0);
    }

    /// @brief Converts bin index into a valid one for this axis.
    ///
    /// @note Bound: bin index is clamped to [1, nBins]
    ///
    /// @param [in] bin The bin to wrap
    /// @return valid bin index
    template <AxisBoundaryType T = bdt,
              std::enable_if_t<T == AxisBoundaryType::Bound, int> = 0>
    size_t
    wrapBin(int bin) const
    {
      return std::max(std::min(bin, static_cast<int>(getNBins())), 1);
    }

    /// @brief Converts bin index into a valid one for this axis.
    ///
    /// @note Closed: bin index wraps around to other side
    ///
    /// @param [in] bin The bin to wrap
    /// @return valid bin index
    template <AxisBoundaryType T = bdt,
              std::enable_if_t<T == AxisBoundaryType::Closed, int> = 0>
    size_t
    wrapBin(int bin) const
    {
      const int w = getNBins();
      return 1 + (w + ((bin - 1) % w)) % w;
      // return int(bin<1)*w - int(bin>w)*w + bin;
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
      return wrapBin(std::floor((x - getMin()) / getBinWidth()) + 1);
    }

    /// @brief get bin width
    ///
    /// @return constant width for all bins
    double getBinWidth(size_t /*bin*/ = 0) const { return m_width; }

    /// @brief get lower bound of bin
    ///
    /// @param  [in] bin index of bin
    /// @return lower bin boundary
    ///
    /// @pre @c bin must be a valid bin index (excluding the underflow bin),
    ///      i.e. \f$1 \le \text{bin} \le \text{nBins} + 1\f$
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
    /// @pre @c bin must be a valid bin index (excluding the overflow bin),
    ///      i.e. \f$0 \le \text{bin} \le \text{nBins}\f$
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
    getMax() const override
    {
      return m_max;
    }

    /// @brief get minimum of binning range
    ///
    /// @return minimum of binning range
    double
    getMin() const override
    {
      return m_min;
    }

    /// @brief get total number of bins
    ///
    /// @return total number of bins (excluding under-/overflow bins)
    size_t
    getNBins() const override
    {
      return m_bins;
    }

    /// @brief check whether value is inside axis limits
    ///
    /// @return @c true if \f$\text{xmin} \le x < \text{xmax}\f$, otherwise
    ///         @c false
    ///
    /// @post If @c true is returned, the bin containing the given value is a
    ///       valid bin, i.e. it is neither the underflow nor the overflow bin.
    bool
    isInside(double x) const
    {
      return (m_min <= x) && (x < m_max);
    }

    /// @brief Return a vector of bin edges
    /// @return Vector which contains the bin edges
    std::vector<double>
    getBinEdges() const override
    {
      std::vector<double> binEdges;
      for (size_t i = 1; i <= m_bins; i++) {
        binEdges.push_back(getBinLowerBound(i));
      }
      binEdges.push_back(getBinUpperBound(m_bins));
      return binEdges;
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
  template <AxisBoundaryType bdt>
  class Axis<AxisType::Variable, bdt> final : public IAxis
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

    /// @brief returns whether the axis is equidistante
    ///
    /// @return bool is equidistant
    bool
    isEquidistant() const override
    {
      return false;
    }

    /// @brief returns whether the axis is variable
    ///
    /// @return bool is variable
    bool
    isVariable() const override
    {
      return true;
    }

    /// @brief returns the boundary type set in the template param
    ///
    /// @return @c AxisBoundaryType of this axis
    AxisBoundaryType
    getBoundaryType() const override
    {
      return bdt;
    }

    /// @brief Get #size bins which neighbor the one given
    ///
    /// Generic overload with symmetric size
    ///
    /// @param [in] idx requested bin index
    /// @param [in] size how many neighboring bins
    /// @return std::set of neighboring bin indices (global)
    std::set<size_t>
    neighborHoodIndices(size_t idx, size_t size = 1) const
    {
      return neighborHoodIndices(idx, std::make_pair(size, size));
    }

    /// @brief Get #size bins which neighbor the one given
    ///
    /// This is the version for Open
    ///
    /// @param [in] idx requested bin index
    /// @param [in] sizes how many neighboring bins (up/down)
    /// @return std::set of neighboring bin indices (global)
    /// @note Open varies given bin and allows 0 and NBins+1 (underflow,
    /// overflow)
    ///       as neighbors
    template <AxisBoundaryType T = bdt,
              std::enable_if_t<T == AxisBoundaryType::Open, int> = 0>
    std::set<size_t>
    neighborHoodIndices(size_t idx,
                        std::pair<size_t, size_t> sizes = {1, 1}) const
    {
      std::set<size_t> result;
      constexpr int    min   = 0;
      const int        max   = getNBins() + 1;
      const int        itmin = std::max(min, int(idx - sizes.first));
      const int        itmax = std::min(max, int(idx + sizes.second));
      for (int i = itmin; i <= itmax; i++) {
        result.insert(i);
      }
      return result;
    }

    /// @brief Get #size bins which neighbor the one given
    ///
    /// This is the version for Bound
    ///
    /// @param [in] idx requested bin index
    /// @param [in] sizes how many neighboring bins (up/down)
    /// @return std::set of neighboring bin indices (global)
    /// @note Bound varies given bin and allows 1 and NBins (regular bins)
    ///       as neighbors
    template <AxisBoundaryType T = bdt,
              std::enable_if_t<T == AxisBoundaryType::Bound, int> = 0>
    std::set<size_t>
    neighborHoodIndices(size_t idx,
                        std::pair<size_t, size_t> sizes = {1, 1}) const
    {
      std::set<size_t> result;
      if (idx <= 0 || idx >= (getNBins() + 1)) {
        return result;
      }
      constexpr int min   = 1;
      const int     max   = getNBins();
      const int     itmin = std::max(min, int(idx - sizes.first));
      const int     itmax = std::min(max, int(idx + sizes.second));
      for (int i = itmin; i <= itmax; i++) {
        result.insert(i);
      }
      return result;
    }

    /// @brief Get #size bins which neighbor the one given
    ///
    /// This is the version for Closed
    ///
    /// @param [in] idx requested bin index
    /// @param [in] sizes how many neighboring bins (up/down)
    /// @return std::set of neighboring bin indices (global)
    /// @note Closed varies given bin and allows bins on the opposite
    ///       side of the axis as neighbors. (excludes underflow / overflow)
    template <AxisBoundaryType T = bdt,
              std::enable_if_t<T == AxisBoundaryType::Closed, int> = 0>
    std::set<size_t>
    neighborHoodIndices(size_t idx,
                        std::pair<size_t, size_t> sizes = {1, 1}) const
    {
      std::set<size_t> result;
      if (idx <= 0 || idx >= (getNBins() + 1)) {
        return result;
      }
      const int itmin = idx - sizes.first;
      const int itmax = idx + sizes.second;
      for (int i = itmin; i <= itmax; i++) {
        result.insert(wrapBin(i));
      }
      return result;
    }

    /// @brief Converts bin index into a valid one for this axis.
    ///
    /// @note Open: bin index is clamped to [0, nBins+1]
    ///
    /// @param [in] bin The bin to wrap
    /// @return valid bin index
    template <AxisBoundaryType T = bdt,
              std::enable_if_t<T == AxisBoundaryType::Open, int> = 0>
    size_t
    wrapBin(int bin) const
    {
      return std::max(std::min(bin, static_cast<int>(getNBins()) + 1), 0);
    }

    /// @brief Converts bin index into a valid one for this axis.
    ///
    /// @note Bound: bin index is clamped to [1, nBins]
    ///
    /// @param [in] bin The bin to wrap
    /// @return valid bin index
    template <AxisBoundaryType T = bdt,
              std::enable_if_t<T == AxisBoundaryType::Bound, int> = 0>
    size_t
    wrapBin(int bin) const
    {
      return std::max(std::min(bin, static_cast<int>(getNBins())), 1);
    }

    /// @brief Converts bin index into a valid one for this axis.
    ///
    /// @note Closed: bin index wraps around to other side
    ///
    /// @param [in] bin The bin to wrap
    /// @return valid bin index
    template <AxisBoundaryType T = bdt,
              std::enable_if_t<T == AxisBoundaryType::Closed, int> = 0>
    size_t
    wrapBin(int bin) const
    {
      const int w = getNBins();
      return 1 + (w + ((bin - 1) % w)) % w;
      // return int(bin<1)*w - int(bin>w)*w + bin;
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
      const auto it
          = std::upper_bound(std::begin(m_binEdges), std::end(m_binEdges), x);
      return wrapBin(std::distance(std::begin(m_binEdges), it));
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
    /// @pre @c bin must be a valid bin index (excluding the underflow bin),
    ///      i.e. \f$1 \le \text{bin} \le \text{nBins} + 1\f$
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
    /// @pre @c bin must be a valid bin index (excluding the overflow bin),
    ///      i.e. \f$0 \le \text{bin} \le \text{nBins}\f$
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
    getMax() const override
    {
      return m_binEdges.back();
    }

    /// @brief get minimum of binning range
    ///
    /// @return minimum of binning range
    double
    getMin() const override
    {
      return m_binEdges.front();
    }

    /// @brief get total number of bins
    ///
    /// @return total number of bins (excluding under-/overflow bins)
    size_t
    getNBins() const override
    {
      return m_binEdges.size() - 1;
    }

    /// @brief check whether value is inside axis limits
    ///
    /// @return @c true if \f$\text{xmin} \le x < \text{xmax}\f$, otherwise
    ///         @c false
    ///
    /// @post If @c true is returned, the bin containing the given value is a
    ///       valid bin, i.e. it is neither the underflow nor the overflow bin.
    bool
    isInside(double x) const
    {
      return (m_binEdges.front() <= x) && (x < m_binEdges.back());
    }

    /// @brief Return a vector of bin edges
    /// @return Vector which contains the bin edges
    std::vector<double>
    getBinEdges() const override
    {
      return m_binEdges;
    }

  private:
    /// vector of bin edges (sorted in ascending order)
    std::vector<double> m_binEdges;
  };
}  // namespace detail

}  // namespace Acts
